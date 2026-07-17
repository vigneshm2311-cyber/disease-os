"""
core/sources/chembl.py

ChEMBL 37 drug-target mechanism loader for Disease-OS.

Sources:
  - ChEMBL REST API /mechanism endpoint (7,561 approved drug-target pairs)
  - data/raw/chembl_uniprot_mapping.txt (17,258 ChEMBL target → UniProt)

Adds:
  - TARGETS edges: Drug → Protein (UniProt) with action type
  - Action type stored as source_relationship_type
    (INHIBITOR / AGONIST / ANTAGONIST / BLOCKER / etc.)
  - PMID evidence on each edge from mechanism_refs
  - ChEMBL molecule ID stored on matched drug nodes

Edge confidence:
  - max_phase=4 (approved)      : 0.90
  - max_phase=3 (phase 3)       : 0.75
  - max_phase=2 (phase 2)       : 0.60
  - max_phase=1 (phase 1)       : 0.45
  - max_phase=0 / null          : 0.30

Run from project root:
    python3 core/sources/chembl.py
"""

import sqlite3, json, time, urllib.request
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
try:
    from core.config import DB_PATH
except ImportError:
    DB_PATH = None

UNIPROT_MAP = Path.home() / "disease-os/data/raw/chembl_uniprot_mapping.txt"
CHEMBL_API  = "https://www.ebi.ac.uk/chembl/api/data"
PAGE_SIZE   = 1000


def load_uniprot_mapping(path: Path) -> dict:
    """Parse chembl_uniprot_mapping.txt → {chembl_target_id: uniprot_acc}"""
    mapping = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                uniprot = parts[0].strip()
                chembl_target = parts[1].strip()
                # Only map SINGLE PROTEIN targets
                if len(parts) < 4 or "SINGLE PROTEIN" in parts[3]:
                    mapping[chembl_target] = uniprot
    return mapping


def fetch_mechanisms() -> list:
    """Fetch all drug mechanisms from ChEMBL API with pagination."""
    all_mechs = []
    offset    = 0
    total     = None

    while True:
        url  = (f"{CHEMBL_API}/mechanism?format=json"
                f"&limit={PAGE_SIZE}&offset={offset}")
        try:
            resp = urllib.request.urlopen(url, timeout=15)
            data = json.loads(resp.read())
        except Exception as e:
            print(f"  API error at offset {offset}: {e}")
            time.sleep(2)
            continue

        mechs = data.get("mechanisms", [])
        if total is None:
            total = data.get("page_meta", {}).get("total_count", 0)
            print(f"  Total mechanisms to fetch: {total:,}")

        all_mechs.extend(mechs)
        offset += PAGE_SIZE

        if offset % 3000 == 0:
            print(f"  Fetched {len(all_mechs):,}/{total:,}...")

        time.sleep(0.3)  # polite rate limiting

        if not data.get("page_meta", {}).get("next"):
            break
        if len(mechs) < PAGE_SIZE:
            break

    return all_mechs


def phase_to_confidence(max_phase) -> float:
    try:
        p = int(max_phase or 0)
    except Exception:
        p = 0
    return {4: 0.90, 3: 0.75, 2: 0.60, 1: 0.45}.get(p, 0.30)


def run(db_path=DB_PATH):
    now  = datetime.now(timezone.utc).isoformat()
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")

    # ── 1. Load UniProt mapping ───────────────────────────────
    print("[1/5] Loading ChEMBL→UniProt mapping...")
    target_to_uniprot = load_uniprot_mapping(UNIPROT_MAP)
    print(f"  {len(target_to_uniprot):,} ChEMBL target → UniProt mappings")

    # ── 2. Build drug name index (same prefix approach as DrugBank) ──
    print("[2/5] Building drug node index...")
    all_drugs = conn.execute("""
        SELECT id, primary_id, primary_system, label, synonyms, confidence
        FROM nodes
        WHERE entity_type IN ('Drug_clinical','Pharmacologic Substance')
    """).fetchall()

    from collections import defaultdict
    prefix_idx = defaultdict(list)
    for node_id, pid, psys, label, syns_json, conf in all_drugs:
        ll = label.lower().strip()
        fw = ll.split()[0] if ll else ""
        if fw:
            prefix_idx[fw].append((ll, pid, psys, node_id, conf))
        if syns_json:
            try:
                for syn in json.loads(syns_json):
                    sl = syn.lower().strip()
                    sfw = sl.split()[0] if sl else ""
                    if sfw:
                        prefix_idx[sfw].append((sl, pid, psys, node_id, conf))
            except Exception:
                pass

    def find_drug(name: str):
        nl = name.lower().strip()
        fw = nl.split()[0] if nl else ""
        for ll, pid, psys, nid, conf in prefix_idx.get(fw, []):
            if ll == nl:
                return pid, psys, nid
        best = None
        for ll, pid, psys, nid, conf in prefix_idx.get(fw, []):
            if ll.startswith(nl + " ") or ll.startswith(nl + ","):
                if best is None or conf > best[4]:
                    best = (ll, pid, psys, nid, conf)
        return (best[1], best[2], best[3]) if best else None

    print(f"  {sum(len(v) for v in prefix_idx.values()):,} drug entries indexed")

    # Also build ChEMBL molecule ID → graph node index via API lookup cache
    chembl_mol_cache = {}  # chembl_mol_id → (pid, psys, node_id)

    def resolve_chembl_mol(chembl_mol_id: str):
        if chembl_mol_id in chembl_mol_cache:
            return chembl_mol_cache[chembl_mol_id]
        try:
            url  = (f"{CHEMBL_API}/molecule/{chembl_mol_id}"
                    f"?format=json")
            resp = urllib.request.urlopen(url, timeout=8)
            data = json.loads(resp.read())
            pref = data.get("pref_name") or ""
            syns = [s.get("molecule_synonym","")
                    for s in data.get("molecule_synonyms", [])]
            # Try preferred name first, then synonyms
            for name in [pref] + syns:
                if name:
                    hit = find_drug(name)
                    if hit:
                        chembl_mol_cache[chembl_mol_id] = hit
                        return hit
        except Exception:
            pass
        chembl_mol_cache[chembl_mol_id] = None
        return None

    # ── 3. Fetch all mechanisms from API ─────────────────────
    print("[3/5] Fetching drug mechanisms from ChEMBL API...")
    mechanisms = fetch_mechanisms()
    print(f"  {len(mechanisms):,} mechanisms fetched")

    # ── 4. Build edges ────────────────────────────────────────
    print("[4/5] Resolving drug→UniProt edges...")

    edge_cols = [r[1] for r in
                 conn.execute("PRAGMA table_info(edges)").fetchall()
                 if r[1] != "id"]
    col = {c: i for i, c in enumerate(edge_cols)}
    ph  = ",".join(["?"] * len(edge_cols))
    cs  = ",".join(edge_cols)

    def make_edge(src_id, src_sys, tgt_id, tgt_sys,
                  rel, src_rel, conf, source, pmids):
        e = [None] * len(edge_cols)
        e[col["source_id"]]               = src_id
        e[col["source_system"]]           = src_sys
        e[col["target_id"]]               = tgt_id
        e[col["target_system"]]           = tgt_sys
        e[col["relationship_type"]]       = rel
        e[col["source_relationship_type"]]= src_rel
        e[col["confidence"]]              = conf
        e[col["primary_source"]]          = source
        e[col["imported_via"]]            = "ChEMBL_37_API"
        e[col["study_design"]]            = "clinical_review"
        e[col["source_version"]]          = "37"
        e[col["loaded_at"]]               = now
        e[col["feedback"]]                = 0
        return tuple(e)

    target_edges = []
    seen         = set()
    resolved     = 0
    no_uniprot   = 0
    no_drug      = 0
    api_calls    = 0

    for i, mech in enumerate(mechanisms):
        chembl_mol    = mech.get("molecule_chembl_id", "")
        chembl_target = mech.get("target_chembl_id", "")
        action_type   = mech.get("action_type", "") or "TARGETS"
        moa_text      = mech.get("mechanism_of_action", "") or ""
        max_phase     = mech.get("max_phase", 0)
        conf          = phase_to_confidence(max_phase)

        # Get UniProt for this target
        uniprot = target_to_uniprot.get(chembl_target)
        if not uniprot:
            no_uniprot += 1
            continue

        # Check UniProt node exists in graph
        up_row = conn.execute(
            "SELECT primary_id, primary_system FROM nodes "
            "WHERE primary_id=? AND primary_system='UniProt' LIMIT 1",
            (uniprot,)
        ).fetchone()
        if not up_row:
            no_uniprot += 1
            continue

        # Resolve drug node — first try cache, then API
        drug_node = chembl_mol_cache.get(chembl_mol)
        if drug_node is None and chembl_mol not in chembl_mol_cache:
            api_calls += 1
            drug_node = resolve_chembl_mol(chembl_mol)
            if api_calls % 100 == 0:
                time.sleep(0.5)  # rate limit

        if not drug_node:
            no_drug += 1
            continue

        drug_pid, drug_psys, drug_nid = drug_node
        pair = (drug_pid, uniprot)
        if pair in seen:
            continue
        seen.add(pair)

        # Extract PMIDs from mechanism_refs
        pmids = [r["ref_id"] for r in mech.get("mechanism_refs", [])
                 if r.get("ref_type") == "PubMed"]
        source = pmids[0] if pmids else chembl_mol

        # Relationship type from action_type
        rel = "TARGETS"

        # Store action type as source_relationship_type
        src_rel = f"chembl_{action_type.lower()}"

        target_edges.append(make_edge(
            drug_pid, drug_psys,
            uniprot, "UniProt",
            rel, src_rel, conf, source, pmids
        ))
        resolved += 1

        if (i + 1) % 500 == 0:
            print(f"  {i+1:,}/{len(mechanisms):,} processed | "
                  f"resolved={resolved} no_uniprot={no_uniprot} "
                  f"no_drug={no_drug} api_calls={api_calls}")

    print(f"\n  Done: {resolved} TARGETS edges | "
          f"no_uniprot={no_uniprot} | no_drug={no_drug}")

    # ── 5. Commit ─────────────────────────────────────────────
    print("[5/5] Committing to database...")
    inserted = 0
    conn.execute("BEGIN")
    for e in target_edges:
        try:
            conn.execute(
                f"INSERT OR IGNORE INTO edges ({cs}) VALUES ({ph})", e
            )
            inserted += 1
        except Exception:
            pass
    conn.execute("COMMIT")
    print(f"  {inserted:,} TARGETS edges inserted")

    # Validation
    total_targets = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='TARGETS'"
    ).fetchone()[0]
    print(f"\n  Total TARGETS edges in graph: {total_targets:,}")

    sample = conn.execute("""
        SELECT nd.label, np.label, e.source_relationship_type, e.confidence
        FROM edges e
        JOIN nodes nd ON nd.primary_id=e.source_id AND nd.primary_system=e.source_system
        JOIN nodes np ON np.primary_id=e.target_id AND np.primary_system=e.target_system
        WHERE e.relationship_type='TARGETS'
        LIMIT 10
    """).fetchall()
    print(f"\n  Sample TARGETS edges:")
    for drug, protein, action, conf in sample:
        print(f"    {drug[:35]:<37} --[{action}]--> {protein[:35]} ({conf})")

    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.close()
    print("\n✓ ChEMBL load complete")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    run()
