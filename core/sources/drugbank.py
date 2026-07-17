"""
core/sources/drugbank.py

DrugBank 5.1.10 CSV loader for Disease-OS.
Input: data/raw/drugbank/drugbank_clean.csv (40MB)

Adds:
  - MOA annotations to 3,954 existing drug nodes
  - 1,190 TREATS edges from indication text matching
  - DrugBank ID xref on matched drug nodes

Notes:
  - Prefix-based name matching (ingredient -> formulation)
  - TREATS edges from keyword matching on indication text
  - Approved drugs only
  - Confidence 0.80 on TREATS edges
  - BE target IDs stored as properties for future UniProt resolution
    (requires full DrugBank XML for BE->UniProt mapping)

Run from project root:
    python3 core/sources/drugbank.py
"""

import csv, sqlite3, json
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
from core.config import DB_PATH

CSV_PATH = Path.home() / "disease-os/data/raw/drugbank/drugbank_clean.csv"

DISEASE_SYNONYMS = {
    "E11.9":  ["type 2 diabetes","T2D","T2DM","diabetes mellitus type 2","type II diabetes"],
    "E10.9":  ["type 1 diabetes","T1D","T1DM","juvenile diabetes"],
    "G30":    ["alzheimer's disease","alzheimer disease"],
    "I10":    ["essential hypertension","high blood pressure"],
    "F32.9":  ["major depressive disorder","major depression","MDD"],
    "I25.10": ["coronary artery disease","CAD","ischemic heart disease"],
    "J45":    ["asthma","bronchial asthma"],
    "N18.9":  ["chronic kidney disease","CKD","chronic renal failure"],
    "I50.9":  ["heart failure","congestive heart failure","CHF"],
    "J44.9":  ["COPD","chronic obstructive pulmonary disease"],
    "E78.5":  ["hyperlipidemia","dyslipidemia","hypercholesterolemia"],
    "M06.9":  ["rheumatoid arthritis","RA"],
    "G35":    ["multiple sclerosis","MS"],
    "L40.0":  ["psoriasis","plaque psoriasis"],
    "G20":    ["parkinson's disease","parkinson disease"],
}

DISEASE_KW = {
    "type 2 diabetes":    "E11.9",
    "type ii diabetes":   "E11.9",
    "diabetes mellitus":  "E11.9",
    "type 1 diabetes":    "E10.9",
    "type i diabetes":    "E10.9",
    "rheumatoid arthritis": "M06.9",
    "multiple sclerosis": "G35",
    "crohn":              "K50.90",
    "inflammatory bowel": "K51.90",
    "psoriasis":          "L40.0",
    "alzheimer":          "G30",
    "parkinson":          "G20",
    "hypertension":       "I10",
    "heart failure":      "I50.9",
    "coronary":           "I25.10",
    "angina":             "I25.10",
    "asthma":             "J45",
    "copd":               "J44.9",
    "depression":         "F32.9",
    "schizophrenia":      "F20",
    "epilepsy":           "G40.909",
    "hepatitis c":        "B18.2",
    "hepatitis b":        "B18.1",
    "hiv":                "B20",
    "breast cancer":      "C50.9",
    "prostate cancer":    "C61",
    "lung cancer":        "C34.10",
    "leukemia":           "C91.0",
    "lymphoma":           "C82.90",
    "anemia":             "D61.9",
    "anaemia":            "D61.9",
    "hemophilia":         "D66",
    "osteoporosis":       "M81.0",
    "chronic kidney":     "N18.9",
    "renal failure":      "N18.9",
    "cystic fibrosis":    "J98.09",
    "gaucher":            "E75.22",
    "acromegaly":         "E22.0",
    "infertility":        "N97.9",
    "hyperlipidemia":     "E78.5",
    "hypercholesterol":   "E78.0",
    "gout":               "M10.9",
    "transplant rejection": "T86.19",
    "neutropenia":        "D70.9",
    "thrombocytopenia":   "D69.6",
    "sepsis":             "A41.9",
    "osteoarthritis":     "M19.90",
    "hypothyroid":        "E03.9",
    "hyperthyroid":       "E05.90",
}


def run(db_path=DB_PATH, csv_path=CSV_PATH):
    now  = datetime.now(timezone.utc).isoformat()
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")

    # ── Add disease synonyms ──────────────────────────────────
    print("[1/4] Adding disease synonyms...")
    updated = 0
    for icd10, syns in DISEASE_SYNONYMS.items():
        row = conn.execute(
            "SELECT id, synonyms FROM nodes WHERE primary_id=? LIMIT 1",
            (icd10,)
        ).fetchone()
        if row:
            existing = json.loads(row[1] or "[]")
            el = {s.lower() for s in existing}
            added = [s for s in syns if s.lower() not in el]
            if added:
                conn.execute(
                    "UPDATE nodes SET synonyms=? WHERE id=?",
                    (json.dumps(existing + added), row[0])
                )
                updated += 1
    conn.commit()
    print(f"    {updated} disease nodes updated")

    # ── Build prefix index ────────────────────────────────────
    print("[2/4] Building drug prefix index...")
    all_drugs = conn.execute("""
        SELECT id, primary_id, primary_system, label, synonyms, confidence
        FROM nodes
        WHERE entity_type IN ('Drug_clinical','Pharmacologic Substance')
    """).fetchall()

    prefix_index = defaultdict(list)
    for node_id, pid, psys, label, syns_json, conf in all_drugs:
        ll = label.lower().strip()
        fw = ll.split()[0] if ll else ""
        if fw:
            prefix_index[fw].append((ll, pid, psys, node_id, conf))
        if syns_json:
            try:
                for syn in json.loads(syns_json):
                    sl = syn.lower().strip()
                    sfw = sl.split()[0] if sl else ""
                    if sfw:
                        prefix_index[sfw].append((sl, pid, psys, node_id, conf))
            except Exception:
                pass

    def find_drug(name):
        nl = name.lower().strip()
        fw = nl.split()[0] if nl else ""
        for ll, pid, psys, nid, conf in prefix_index.get(fw, []):
            if ll == nl:
                return pid, psys, nid
        best = None
        for ll, pid, psys, nid, conf in prefix_index.get(fw, []):
            if ll.startswith(nl + " ") or ll.startswith(nl + ","):
                if best is None or conf > best[4]:
                    best = (ll, pid, psys, nid, conf)
        return (best[1], best[2], best[3]) if best else None

    print(f"    {len(prefix_index):,} prefix keys indexed")

    # ── Resolve disease keywords ──────────────────────────────
    kw_map = {}
    for kw, icd in DISEASE_KW.items():
        row = conn.execute(
            "SELECT primary_id, primary_system FROM nodes "
            "WHERE primary_id=? LIMIT 1", (icd,)
        ).fetchone()
        if row:
            kw_map[kw] = (row[0], row[1])
    print(f"    {len(kw_map)}/{len(DISEASE_KW)} disease keywords resolved")

    # ── Edge schema ───────────────────────────────────────────
    edge_cols = [r[1] for r in
                 conn.execute("PRAGMA table_info(edges)").fetchall()
                 if r[1] != "id"]
    col = {c: i for i, c in enumerate(edge_cols)}
    ph  = ",".join(["?"] * len(edge_cols))
    cs  = ",".join(edge_cols)

    def make_edge(s, ss, t, ts, rel, conf, src):
        e = [None] * len(edge_cols)
        e[col["source_id"]]               = s
        e[col["source_system"]]           = ss
        e[col["target_id"]]               = t
        e[col["target_system"]]           = ts
        e[col["relationship_type"]]       = rel
        e[col["source_relationship_type"]]= f"drugbank_{rel.lower()}"
        e[col["confidence"]]              = conf
        e[col["primary_source"]]          = src
        e[col["imported_via"]]            = "DrugBank_5.1.10_CSV"
        e[col["study_design"]]            = "clinical_review"
        e[col["source_version"]]          = "5.1.10"
        e[col["loaded_at"]]               = now
        e[col["feedback"]]                = 0
        return tuple(e)

    # ── Parse CSV ─────────────────────────────────────────────
    print("[3/4] Parsing DrugBank CSV...")
    moa_updates = []
    treat_edges = []
    seen        = set()
    matched     = 0

    with open(csv_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            db_id      = row.get("drugbank-id", "").strip()
            name       = row.get("name", "").strip()
            if not db_id or not name:
                continue
            hit = find_drug(name)
            if not hit:
                continue
            pid, psys, nid = hit
            matched += 1

            moa        = row.get("mechanism-of-action", "").strip()
            desc       = row.get("description", "").strip()
            targets    = row.get("targets", "").strip()
            groups     = row.get("groups", "").lower()
            indication = row.get("indication", "").strip()

            pr = conn.execute(
                "SELECT properties FROM nodes WHERE id=?", (nid,)
            ).fetchone()
            if pr:
                props = json.loads(pr[0] or "{}")
                props["drugbank_id"] = db_id
                if moa and not props.get("moa"):
                    props["moa"] = moa[:600]
                if desc and not props.get("description"):
                    props["description"] = desc[:400]
                if targets:
                    props["drugbank_targets"] = targets[:300]
                moa_updates.append((json.dumps(props), nid))

            if indication and "approved" in groups:
                ind_l = indication.lower()
                for kw, (dpid, dpsys) in kw_map.items():
                    if kw in ind_l and (pid, dpid) not in seen:
                        seen.add((pid, dpid))
                        treat_edges.append(
                            make_edge(pid, psys, dpid, dpsys,
                                      "TREATS", 0.80, db_id)
                        )

    print(f"    Matched={matched:,}  MOA={len(moa_updates):,}  "
          f"TREATS={len(treat_edges):,}")

    # ── Commit ────────────────────────────────────────────────
    print("[4/4] Committing to database...")
    conn.execute("BEGIN")
    conn.executemany(
        "UPDATE nodes SET properties=? WHERE id=?", moa_updates
    )
    conn.execute("COMMIT")

    inserted = 0
    conn.execute("BEGIN")
    for e in treat_edges:
        try:
            conn.execute(
                f"INSERT OR IGNORE INTO edges ({cs}) VALUES ({ph})", e
            )
            inserted += 1
        except Exception:
            pass
    conn.execute("COMMIT")
    print(f"    {inserted:,} TREATS edges inserted")

    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.close()
    print("✓ DrugBank load complete")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    run()
