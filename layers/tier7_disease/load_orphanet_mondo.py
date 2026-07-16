"""
layers/tier7_disease/load_orphanet_mondo.py

Loads Orphanet and MONDO rare disease nodes into Disease-OS,
then rebuilds ClinVar edges for previously isolated pathogenic variants.

Run from project root:
    python3 layers/tier7_disease/load_orphanet_mondo.py
"""

import sys, csv, sqlite3, json
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, ".")
from core.graph_store import GraphStore
from core.node        import Node
from core.edge        import Edge
from core.config      import DB_PATH
from core.sources.orphanet import OrphanetSource
from core.sources.mondo    import MondoSource

CLINVAR = Path.home() / "disease-os/data/raw/clinvar/variant_summary.txt"
BATCH   = 3_000


def bulk_nodes(conn, nodes):
    if not nodes: return
    sample = Node("x","x","x",0,"Gene")
    cols   = list(sample.to_dict().keys())
    ph     = ",".join(["?"]*len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO nodes ({','.join(cols)}) VALUES ({ph})",
            [list(n.to_dict().values()) for n in nodes]
        )
        conn.execute("COMMIT")
    except Exception as e:
        conn.execute("ROLLBACK"); raise e


def bulk_edges(conn, edges):
    if not edges: return
    sample = Edge("x","x","y","y","ISA","x","x","x")
    cols   = list(sample.to_dict().keys())
    ph     = ",".join(["?"]*len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO edges ({','.join(cols)}) VALUES ({ph})",
            [list(e.to_dict().values()) for e in edges]
        )
        conn.execute("COMMIT")
    except Exception as e:
        conn.execute("ROLLBACK"); raise e


def rebuild_clinvar_edges(conn):
    """
    Rebuild ClinVar variant→disease edges for previously isolated
    pathogenic variants using Orphanet/MONDO nodes now in graph.
    """
    print(f"\n[3/3] Rebuilding ClinVar edges for isolated variants...")

    # Build lookup indices from newly loaded nodes
    # Orphanet ID → (primary_id, primary_system)
    orphanet_idx = {}
    for row in conn.execute("""
        SELECT primary_id, primary_system FROM nodes
        WHERE primary_system = 'Orphanet'
    """).fetchall():
        # Strip "Orphanet:" prefix to get bare code for matching
        code = row[0].replace("Orphanet:","")
        orphanet_idx[code]      = (row[0], row[1])
        orphanet_idx[row[0]]    = (row[0], row[1])  # full form too

    # MONDO ID → (primary_id, primary_system)
    mondo_idx = {}
    for row in conn.execute("""
        SELECT primary_id, primary_system FROM nodes
        WHERE primary_system = 'MONDO_disease'
    """).fetchall():
        mondo_idx[row[0]] = (row[0], row[1])

    # MedGen via xrefs
    medgen_idx = {}
    for row in conn.execute("""
        SELECT primary_id, primary_system, xrefs FROM nodes
        WHERE xrefs LIKE '%MedGen%'
          AND entity_type IN ('Disease_clinical','Disease_scientific')
    """).fetchall():
        try:
            xrefs = json.loads(row[2] or "{}")
            mg    = xrefs.get("MedGen")
            if mg:
                medgen_idx[mg] = (row[0], row[1])
        except: pass

    print(f"  Orphanet idx : {len(orphanet_idx):,}")
    print(f"  MONDO idx    : {len(mondo_idx):,}")
    print(f"  MedGen idx   : {len(medgen_idx):,}")

    # Get isolated rsID variants
    isolated = set(r[0] for r in conn.execute("""
        SELECT primary_id FROM nodes
        WHERE primary_system='dbSNP_rsID'
          AND NOT EXISTS (
            SELECT 1 FROM edges e WHERE e.source_id=nodes.primary_id
          )
    """).fetchall())
    print(f"  Isolated rsID variants: {len(isolated):,}")

    def find_disease(pheno_ids_str):
        if not pheno_ids_str:
            return None
        first = pheno_ids_str.split("|")[0]
        for part in first.split(","):
            part = part.strip()
            if ":" not in part:
                continue
            sys_raw, _, id_raw = part.partition(":")
            sys_raw = sys_raw.strip()
            id_raw  = id_raw.strip()
            if "MONDO" in sys_raw:
                # Handle MONDO:MONDO:0001234 → MONDO:0001234
                mondo_id = f"MONDO:{id_raw}" if not id_raw.startswith("MONDO:") \
                           else id_raw
                hit = mondo_idx.get(mondo_id)
                if hit: return hit
            elif "rphanet" in sys_raw:
                hit = orphanet_idx.get(id_raw) or orphanet_idx.get(f"Orphanet:{id_raw}")
                if hit: return hit
            elif sys_raw == "MedGen":
                hit = medgen_idx.get(id_raw)
                if hit: return hit
            elif sys_raw == "OMIM":
                hit = conn.execute(
                    "SELECT primary_id,primary_system FROM nodes "
                    "WHERE primary_id=? AND primary_system='OMIM' LIMIT 1",
                    (id_raw,)
                ).fetchone()
                if hit: return hit
        return None

    def clinsig_to_edge(sig):
        s = sig.lower()
        if "pathogenic/likely pathogenic" in s: return ("CAUSES",0.80)
        if "pathogenic" in s and "likely" not in s: return ("CAUSES",0.85)
        if "likely pathogenic" in s: return ("CONTRIBUTES_TO",0.70)
        if "risk factor" in s or "risk allele" in s: return ("INCREASES_RISK_OF",0.60)
        if "protective" in s: return ("PROTECTS_AGAINST",0.70)
        return None

    EDGE_COLS = [r[1] for r in conn.execute("PRAGMA table_info(edges)").fetchall()
                 if r[1] != "id"]
    PH      = ",".join(["?"]*len(EDGE_COLS))
    COL_STR = ",".join(EDGE_COLS)
    src_idx = EDGE_COLS.index("source_id")
    tgt_idx = EDGE_COLS.index("target_id")
    tsys_idx= EDGE_COLS.index("target_system")
    rel_idx = EDGE_COLS.index("relationship_type")
    srel_idx= EDGE_COLS.index("source_relationship_type")
    conf_idx= EDGE_COLS.index("confidence")
    ps_idx  = EDGE_COLS.index("primary_source")
    iv_idx  = EDGE_COLS.index("imported_via")
    sd_idx  = EDGE_COLS.index("study_design")
    sv_idx  = EDGE_COLS.index("source_version")
    la_idx  = EDGE_COLS.index("loaded_at")

    now      = datetime.now(timezone.utc).isoformat()
    rebuilt  = no_match = seen = 0
    batch    = []
    seen_set = set()

    with open(CLINVAR, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 16 or row[16] != "GRCh38":
                continue
            rs_raw = row[9].strip()
            if not rs_raw or rs_raw in ("-1","0"):
                continue
            rsid = f"rs{rs_raw}"
            if rsid not in isolated:
                continue

            clinsig   = row[6].strip()
            pheno_ids = row[12].strip()
            rcv       = row[11].strip().split(";")[0]

            edge_info = clinsig_to_edge(clinsig)
            if not edge_info:
                continue

            rel_type, confidence = edge_info
            disease = find_disease(pheno_ids)
            if not disease:
                no_match += 1
                continue

            pair = (rsid, disease[0], rel_type)
            if pair in seen_set:
                seen += 1
                continue
            seen_set.add(pair)

            e = [None] * len(EDGE_COLS)
            e[src_idx]  = rsid
            e[EDGE_COLS.index("source_system")] = "dbSNP_rsID"
            e[tgt_idx]  = disease[0]
            e[tsys_idx] = disease[1]
            e[rel_idx]  = rel_type
            e[srel_idx] = f"clinvar_{clinsig[:30].lower()}"
            e[conf_idx] = confidence
            e[ps_idx]   = rcv or f"ClinVar_{row[0]}"
            e[iv_idx]   = "ClinVar_variant_summary_2024-06"
            e[sd_idx]   = "clinical_review"
            e[sv_idx]   = "2024-06"
            e[la_idx]   = now
            batch.append(tuple(e))
            rebuilt += 1

            if len(batch) >= BATCH:
                conn.execute("BEGIN")
                conn.executemany(
                    f"INSERT OR IGNORE INTO edges ({COL_STR}) VALUES ({PH})",
                    batch
                )
                conn.execute("COMMIT")
                batch = []
                print(f"  {rebuilt:,} edges rebuilt | "
                      f"no_match={no_match:,} | dupes={seen:,}")

    if batch:
        conn.execute("BEGIN")
        conn.executemany(
            f"INSERT OR IGNORE INTO edges ({COL_STR}) VALUES ({PH})", batch
        )
        conn.execute("COMMIT")

    return rebuilt, no_match


if __name__ == "__main__":
    print(f"[Disease-OS] Loading Orphanet + MONDO + rebuilding ClinVar edges")
    print(f"  DB: {DB_PATH}")

    gs   = GraphStore(DB_PATH)
    conn = gs.conn
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    t0 = datetime.now()

    # ── 1. Orphanet ────────────────────────────────────────────────────────
    print(f"\n[1/3] Loading Orphanet...")
    src   = OrphanetSource()
    batch = []
    n_nodes = 0
    for node in src.nodes():
        batch.append(node)
        if len(batch) >= BATCH:
            bulk_nodes(conn, batch)
            n_nodes += len(batch)
            batch    = []
    if batch:
        bulk_nodes(conn, batch)
        n_nodes += len(batch)
    print(f"  {n_nodes:,} Orphanet nodes loaded")

    # ── 2. MONDO ───────────────────────────────────────────────────────────
    print(f"\n[2/3] Loading MONDO...")
    src2    = MondoSource()
    batch   = []
    n_nodes2 = 0
    n_edges2 = 0
    for node in src2.nodes():
        batch.append(node)
        if len(batch) >= BATCH:
            bulk_nodes(conn, batch)
            n_nodes2 += len(batch)
            batch     = []
    if batch:
        bulk_nodes(conn, batch)
        n_nodes2 += len(batch)
    print(f"  {n_nodes2:,} MONDO nodes loaded")

    edge_batch = []
    for edge in src2.edges():
        edge_batch.append(edge)
        if len(edge_batch) >= BATCH:
            bulk_edges(conn, edge_batch)
            n_edges2    += len(edge_batch)
            edge_batch   = []
    if edge_batch:
        bulk_edges(conn, edge_batch)
        n_edges2 += len(edge_batch)
    print(f"  {n_edges2:,} MONDO ISA edges loaded")

    # ── 3. Rebuild ClinVar edges ───────────────────────────────────────────
    rebuilt, no_match = rebuild_clinvar_edges(conn)

    conn.execute("PRAGMA synchronous=NORMAL;")
    elapsed = int((datetime.now()-t0).total_seconds())

    # ── Validation ─────────────────────────────────────────────────────────
    print(f"\n── Results ──────────────────────────────────────────────")
    print(f"  Orphanet nodes   : {n_nodes:,}")
    print(f"  MONDO nodes      : {n_nodes2:,}")
    print(f"  MONDO ISA edges  : {n_edges2:,}")
    print(f"  ClinVar edges rebuilt: {rebuilt:,}")
    print(f"  Still no match   : {no_match:,}")
    print(f"  Time             : {elapsed}s")

    total_nodes = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    total_edges = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    still_iso   = conn.execute("""
        SELECT COUNT(*) FROM nodes
        WHERE primary_system='dbSNP_rsID'
          AND NOT EXISTS (SELECT 1 FROM edges e WHERE e.source_id=primary_id)
    """).fetchone()[0]

    print(f"\n  Total nodes      : {total_nodes:,}")
    print(f"  Total edges      : {total_edges:,}")
    print(f"  rsID still isolated: {still_iso:,}")

    # Spot check
    foxred = conn.execute("""
        SELECT target_id, target_system, relationship_type, confidence
        FROM edges WHERE source_id='rs267606830' LIMIT 3
    """).fetchall()
    print(f"\n  rs267606830 (FOXRED1) edges: {len(foxred)}")
    for r in foxred:
        print(f"    → {r[0]} ({r[1]}) [{r[2]}] conf={r[3]}")

    print(f"\n✓  Orphanet + MONDO + ClinVar rebuild complete")
