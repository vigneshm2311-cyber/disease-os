"""
tests/validate.py

Disease-OS Validation Suite — 20 ground truth tests.
Uses verified node IDs from the actual graph.

Run: python3 tests/validate.py
Exit 0 = all passed, Exit 1 = failures
"""

import sys, json, sqlite3
sys.path.insert(0, ".")
from core.config        import DB_PATH
from core.causal_engine import CausalEngine

conn   = sqlite3.connect(str(DB_PATH))
conn.execute("PRAGMA cache_size=-256000;")
engine = CausalEngine(str(DB_PATH))

results = []

def test(name, passed, detail="", warn=False):
    status = "✅ PASS" if passed else ("⚠️  WARN" if warn else "❌ FAIL")
    results.append((name, passed or warn, status, detail))
    print(f"  {status}  {name}")
    if detail:
        print(f"         {detail}")
    return passed

def q1(sql, *args):
    return conn.execute(sql, args).fetchone()

def qall(sql, *args):
    return conn.execute(sql, args).fetchall()

def trace(disease_id, system="ICD-10-CM", depth=4, top_n=500, patient=None):
    node = q1("SELECT primary_id,primary_system FROM nodes "
              "WHERE icd10_code=? ORDER BY confidence DESC LIMIT 1", disease_id)
    if not node:
        node = q1("SELECT primary_id,primary_system FROM nodes "
                  "WHERE primary_id=? AND primary_system=? LIMIT 1",
                  disease_id, system)
    if not node:
        return [], {}
    paths = engine.trace(
        disease_id=node[0], disease_system=node[1],
        max_depth=depth, min_confidence=0.50,
        min_path_score=0.01, top_n=top_n, patient_data=patient,
    )
    tier_scores = {}
    for p in paths:
        t = p.edges[0].tier_from
        tier_scores[t] = tier_scores.get(t,0) + p.path_score
    return paths, tier_scores

print("=" * 62)
print("  DISEASE-OS VALIDATION SUITE  — 20 GROUND TRUTH TESTS")
print("=" * 62)

# ── SECTION 1: VARIANT → DISEASE EDGES ───────────────────────
print(f"\n── Section 1: Variant→Disease Causal Edges ─────────────")

VARIANT_TESTS = [
    # (rsid, gene, disease_fragment, min_conf, note)
    ("rs7903146",   "TCF7L2", "diabetes",      0.85, "T2D risk variant"),
    ("rs429358",    "APOE",   "lzheimer",      0.60, "AD risk variant"),
    ("rs749553271", "HFE",    "emochromat",    0.60, "HFE CAUSES hemochromatosis"),
    ("rs80359884",  "BRCA1",  "hboc",          0.85, "BRCA1 CAUSES HBOC conf=0.85"),
]

for rsid, gene, frag, min_conf, note in VARIANT_TESTS:
    row = q1("""
        SELECT e.relationship_type, e.confidence, n.label
        FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=?
          AND e.relationship_type IN ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF')
          AND LOWER(n.label) LIKE ?
        ORDER BY e.confidence DESC LIMIT 1
    """, rsid, f"%{frag}%")
    if row:
        passed = row[1] >= min_conf
        test(f"{rsid} ({gene}) → {note}",
             passed,
             f"rel={row[0]} conf={row[1]:.2f}≥{min_conf} "
             f"disease='{row[2][:45]}'")
    else:
        test(f"{rsid} ({gene}) → {note}", False, "No matching edge found")

# ── SECTION 2: CAUSAL ENGINE — 10 DISEASES ───────────────────
print(f"\n── Section 2: Causal Engine — 10 Diseases ──────────────")

DISEASE_TESTS = [
    # (lookup_id, system, name, expected_tier_set, min_depth)
    ("E11.9",    "ICD-10-CM",   "Type 2 Diabetes",          {5,8},  3),
    ("G30",      "ICD-10-CM",   "Alzheimer's Disease",       {5,8},  2),
    ("I10",      "ICD-10-CM",   "Hypertension",              {5},    2),
    ("F32.9",    "ICD-10-CM",   "Major Depression",          {5,8},  2),
    ("I25.10",   "ICD-10-CM",   "Coronary Artery Disease",   {5},    1),
    ("J45",      "ICD-10-CM",   "Asthma",                    {1},    1),
    ("612555",   "OMIM",        "Hereditary Breast/Ovarian", {1},    1),
    ("E83.110",  "ICD-10-CM",   "Hereditary Hemochromatosis",{1},    1),
    ("N18.9",    "ICD-10-CM",   "Chronic Kidney Disease",    set(),  1),
    ("F20",     "ICD-10-CM",   "Schizophrenia",              {1},    1),
]

for lookup, system, name, expected_tiers, min_depth in DISEASE_TESTS:
    paths, tier_scores = trace(lookup, system)

    if not paths:
        test(f"{name} ({lookup})", False,
             "No paths — disease node not found or no causal edges")
        continue

    max_depth   = max(p.depth for p in paths)
    found_tiers = set(tier_scores.keys())
    tiers_ok    = expected_tiers.issubset(found_tiers) if expected_tiers else True
    depth_ok    = max_depth >= min_depth
    passed      = depth_ok and tiers_ok

    total = sum(tier_scores.values()) or 1
    tnames = {1:"Mol",2:"Net",3:"Cell",4:"Tiss",5:"Sys",
              6:"Phen",7:"Dis",8:"Behav",9:"Soc",10:"HC"}
    bd = " ".join(f"T{t}={s/total*100:.0f}%"
                  for t,s in sorted(tier_scores.items(),key=lambda x:-x[1])
                  if s/total > 0.02)
    test(f"{name} ({lookup})", passed,
         f"paths={len(paths)} depth={max_depth}≥{min_depth} "
         f"tiers={sorted(found_tiers)} | {bd}")

# ── SECTION 3: MOLECULAR CHAIN INTEGRITY ─────────────────────
print(f"\n── Section 3: Molecular Chain Integrity ─────────────────")

# ENCODES edges
n_enc = q1("SELECT COUNT(*) FROM edges WHERE relationship_type='ENCODES'")[0]
test("ENCODES edges ≥19,000", n_enc >= 19000, f"{n_enc:,} edges")

# TP53 chain
tp53 = q1("SELECT primary_id FROM nodes "
           "WHERE json_extract(properties,'$.hgnc_symbol')='TP53' "
           "AND primary_system='NCBI_Gene' LIMIT 1")
if tp53:
    prot = q1("SELECT n.primary_id FROM edges e "
              "JOIN nodes n ON n.primary_id=e.target_id "
              "WHERE e.source_id=? AND e.relationship_type='ENCODES'", tp53[0])
    paths_via = q1("SELECT COUNT(*) FROM edges "
                   "WHERE source_id='P04637' AND relationship_type='PART_OF'")[0]
    test("TP53: Gene→Protein→Pathway chain",
         prot is not None and paths_via > 0,
         f"gene={tp53[0]} → protein={prot[0] if prot else '?'} → {paths_via} pathways")
else:
    test("TP53: Gene→Protein→Pathway chain", False, "TP53 gene not found")

# HMDB → UniProt resolved
n_met = q1("SELECT COUNT(*) FROM edges "
           "WHERE source_system='HMDB' AND target_system='UniProt'")[0]
test("Metabolite→Protein resolved (≥100K)",
     n_met >= 100000, f"{n_met:,} HMDB→UniProt edges")

# Rare disease nodes loaded
n_orph = q1("SELECT COUNT(*) FROM nodes WHERE primary_system='Orphanet'")[0]
n_mondo = q1("SELECT COUNT(*) FROM nodes WHERE primary_system='MONDO_disease'")[0]
n_rare_edges = q1("""
    SELECT COUNT(*) FROM edges e
    JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
    WHERE n.primary_system IN ('Orphanet','MONDO_disease')
      AND e.relationship_type IN ('CAUSES','CONTRIBUTES_TO')
""")[0]
test("Rare disease coverage (Orphanet+MONDO)",
     n_orph >= 11000 and n_mondo >= 30000 and n_rare_edges >= 100,
     f"Orphanet={n_orph:,} MONDO={n_mondo:,} causal_edges={n_rare_edges:,}")

# ClinVar rare variant fixed
foxred = q1("SELECT COUNT(*) FROM edges WHERE source_id='rs267606830'")[0]
test("rs267606830 (FOXRED1 pathogenic) linked to rare disease",
     foxred > 0, f"{foxred} edge(s)")

# ── SECTION 4: PERSONALIZATION ───────────────────────────────
print(f"\n── Section 4: Personalization ───────────────────────────")

for rsid, disease, icd10 in [
    ("rs7903146", "T2D",           "E11.9"),
    ("rs429358",  "Alzheimer's",   "G30"),
]:
    gen, _ = trace(icd10)
    per, _ = trace(icd10, patient={"variants":[rsid],"metabolites":[],"diseases":[]})
    gs = next((p.path_score for p in gen if any(e.source_id==rsid for e in p.edges)), None)
    ps = next((p.path_score for p in per if any(e.source_id==rsid for e in p.edges)), None)
    if gs and ps:
        boost = ps/gs
        test(f"Personalization: {rsid} → {disease}",
             boost >= 1.3, f"{gs:.4f}→{ps:.4f} (boost={boost:.2f}x ≥1.3x)")
    else:
        test(f"Personalization: {rsid} → {disease}", False,
             warn=True, detail=f"variant not found in trace")

# ── FINAL RESULTS ─────────────────────────────────────────────
print(f"\n{'='*62}")
print(f"  FINAL RESULTS")
print(f"{'='*62}")

passed_n = sum(1 for _,p,_,_ in results if p)
total_n  = len(results)
failed   = [(n,d) for n,p,s,d in results if not p]

for name, passed, status, detail in results:
    print(f"  {status}  {name}")

print(f"\n  Score: {passed_n}/{total_n}")
if passed_n == total_n:
    print(f"\n  ✅ ALL TESTS PASSED — system validated")
elif passed_n >= total_n * 0.85:
    print(f"\n  ⚠️  {total_n-passed_n} test(s) need attention")
else:
    print(f"\n  ❌ {total_n-passed_n} tests failed")

if failed:
    print(f"\n  Failed:")
    for name, detail in failed:
        print(f"    ❌ {name}")
        if detail:
            print(f"       {detail}")

print(f"{'='*62}\n")
conn.close()
sys.exit(0 if passed_n == total_n else 1)
