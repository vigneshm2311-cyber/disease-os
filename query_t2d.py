"""
query_t2d.py

End-to-end causal trace for Type 2 Diabetes.
Demonstrates what the graph can answer right now after Phase 1.

Run from project root:
    python3 query_t2d.py
"""

import sys
import json
sys.path.insert(0, ".")

from core.graph_store import GraphStore
from core.config import DB_PATH

gs   = GraphStore(DB_PATH)
conn = gs.conn

SEP  = "─" * 55

def section(title):
    print(f"\n{SEP}")
    print(f"  {title}")
    print(SEP)

# ── 1. Disease node ────────────────────────────────────────
section("1. DISEASE NODE — Type 2 Diabetes")

t2d = gs.find_by_icd10("E11.9")
if t2d:
    n = t2d[0]
    print(f"  Label      : {n['label']}")
    print(f"  Primary ID : {n['primary_id']} ({n['primary_system']})")
    print(f"  UMLS CUI   : {n['xrefs'].get('UMLS_CUI','n/a')}")
    print(f"  SNOMED     : {n['xrefs'].get('SNOMED','n/a')}")
    print(f"  HCC code   : {n['hcc_code']}  <- insurance risk adjustment")
    print(f"  ICD-11     : {n['icd11_code']}  <- forward compatibility")
    print(f"  Confidence : {n['confidence']}")

# ── 2. What causes T2D — CAUSES edges INTO the disease ────
section("2. CAUSAL VARIANTS → T2D  (CAUSES edges, conf ≥ 0.70)")

# T2D is stored under ICD-10 primary, but ClinVar edges point
# to SNOMED or UMLS_CUI — query both
t2d_snomed = "44054006"
t2d_cui    = "C0011860"

causes_snomed = conn.execute("""
    SELECT source_id, source_system, relationship_type,
           confidence, study_design, source_relationship_type
    FROM edges
    WHERE target_id = ?
      AND relationship_type IN ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF')
      AND confidence >= 0.70
    ORDER BY confidence DESC
    LIMIT 10
""", (t2d_snomed,)).fetchall()

causes_cui = conn.execute("""
    SELECT source_id, source_system, relationship_type,
           confidence, study_design, source_relationship_type
    FROM edges
    WHERE target_id = ?
      AND relationship_type IN ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF')
      AND confidence >= 0.70
    ORDER BY confidence DESC
    LIMIT 10
""", (t2d_cui,)).fetchall()

all_causes = causes_snomed + causes_cui
print(f"  Found {len(all_causes)} high-confidence causal edges (showing top 10)")
for row in all_causes[:10]:
    src_id, src_sys, rel, conf, study, raw = row
    # Look up variant label
    node = gs.get_node(src_id, src_sys)
    label = node['label'][:50] if node else src_id
    props = node['properties'] if node else {}
    gene  = props.get('gene_symbol', '?') if isinstance(props, dict) else '?'
    print(f"\n  [{rel}] conf={conf} study={study}")
    print(f"    Variant : {src_id} ({src_sys})")
    print(f"    Gene    : {gene}")
    print(f"    Label   : {label}")
    print(f"    ClinSig : {raw}")

# ── 3. TCF7L2 full picture ─────────────────────────────────
section("3. TCF7L2 GENE — full picture in the graph")

tcf7l2_variants = conn.execute("""
    SELECT primary_id, primary_system, confidence,
           json_extract(properties,'$.clinsig') as clinsig,
           json_extract(properties,'$.review_status') as review
    FROM nodes
    WHERE entity_type = 'Variant'
      AND json_extract(properties,'$.gene_symbol') = 'TCF7L2'
    ORDER BY confidence DESC
    LIMIT 8
""").fetchall()

print(f"  TCF7L2 variants in graph: "
      f"{conn.execute('SELECT COUNT(*) FROM nodes WHERE entity_type=? AND json_extract(properties,?) = ?', ('Variant','$.gene_symbol','TCF7L2')).fetchone()[0]:,}")
print(f"  Top variants by confidence:")
for v in tcf7l2_variants:
    print(f"    {v[0]:15s} conf={v[2]:.2f}  clinsig={v[3]}  review={v[4]}")

# ── 4. BRCA1 — breast cancer causal chain ─────────────────
section("4. BRCA1 VARIANTS → Breast Cancer")

brca1_count = conn.execute("""
    SELECT COUNT(*) FROM nodes
    WHERE entity_type='Variant'
    AND json_extract(properties,'$.gene_symbol')='BRCA1'
""").fetchone()[0]

brca1_causes = conn.execute("""
    SELECT e.source_id, e.relationship_type, e.confidence,
           n.label as target_label
    FROM edges e
    JOIN nodes n ON n.primary_id = e.target_id
    WHERE e.source_id IN (
        SELECT primary_id FROM nodes
        WHERE entity_type='Variant'
        AND json_extract(properties,'$.gene_symbol')='BRCA1'
        AND confidence >= 0.85
    )
    AND e.relationship_type = 'CAUSES'
    ORDER BY e.confidence DESC
    LIMIT 5
""").fetchall()

print(f"  BRCA1 variants in graph : {brca1_count:,}")
print(f"  High-confidence CAUSES edges (expert-panel reviewed):")
for row in brca1_causes:
    print(f"    {row[0]} --[{row[1]}]--> {row[3][:45]}  conf={row[2]}")

# ── 5. Drug → T2D treatment edges ─────────────────────────
section("5. DRUGS THAT TREAT T2D")

drugs = conn.execute("""
    SELECT e.source_id, e.source_system, e.confidence,
           n.label, n.rxnorm_cui
    FROM edges e
    JOIN nodes n ON n.primary_id = e.source_id
                 AND n.primary_system = e.source_system
    WHERE e.target_id IN (?, ?)
      AND e.relationship_type = 'TREATS'
    ORDER BY e.confidence DESC
    LIMIT 10
""", (t2d_snomed, "E11.9")).fetchall()

print(f"  Found {len(drugs)} drugs with TREATS → T2D edges")
for d in drugs[:8]:
    print(f"  [{d[1]}] {d[3][:50]}")
    print(f"    RxNorm={d[4]}  confidence={d[2]}")

# ── 6. Pathways connected to T2D ──────────────────────────
section("6. PATHWAYS RELEVANT TO T2D")

# Find proteins associated with T2D-related genes,
# then look up which pathways they belong to
t2d_pathways = conn.execute("""
    SELECT DISTINCT n.primary_id, n.label, n.tier
    FROM nodes n
    JOIN edges e ON e.target_id = n.primary_id
    WHERE n.entity_type = 'Pathway'
      AND e.source_id IN (
          SELECT primary_id FROM nodes
          WHERE entity_type='Variant'
          AND json_extract(properties,'$.gene_symbol')
              IN ('TCF7L2','PPARG','KCNJ11','ABCC8','SLC30A8','HNF1A')
      )
    LIMIT 10
""").fetchall()

# Alternative: search pathway labels directly
insulin_pathways = conn.execute("""
    SELECT primary_id, label
    FROM nodes
    WHERE entity_type = 'Pathway'
      AND (label LIKE '%nsulin%'
        OR label LIKE '%lucose%'
        OR label LIKE '%lycoly%'
        OR label LIKE '%eta cell%'
        OR label LIKE '%ancrea%')
    ORDER BY label
    LIMIT 12
""").fetchall()

print(f"  Insulin/glucose/beta-cell pathways in Reactome:")
for p in insulin_pathways:
    print(f"    {p[0]}  {p[1]}")

# ── 7. Causal chain summary ────────────────────────────────
section("7. CAUSAL CHAIN SUMMARY — what the graph knows about T2D")

total_t2d_edges = conn.execute("""
    SELECT relationship_type, COUNT(*) as n
    FROM edges
    WHERE target_id IN (?, ?, 'E11.9')
    GROUP BY relationship_type
    ORDER BY n DESC
""", (t2d_snomed, t2d_cui)).fetchall()

print("  All edge types pointing TO T2D concepts:")
for rel, count in total_t2d_edges:
    print(f"    {rel:30s} : {count:,}")

# ── 8. Graph health check ──────────────────────────────────
section("8. GRAPH HEALTH CHECK")

stats = gs.stats()
print(f"  Total nodes      : {stats['total_nodes']:,}")
print(f"  Total edges      : {stats['total_edges']:,}")
print(f"  Nodes by tier    :")
tier_names = {
    1:"Molecular", 2:"Networks", 3:"Cellular",
    4:"Tissue/Organ", 5:"Systemic", 6:"Phenotype",
    7:"Disease", 8:"Behavior", 9:"Social",
    10:"Healthcare", 11:"Population"
}
for tier, count in sorted(stats['nodes_by_tier'].items()):
    name = tier_names.get(tier, "Unknown")
    bar  = "█" * min(count // 10000, 40)
    print(f"    Tier {tier:2d} {name:12s}: {count:>9,}  {bar}")

print(f"\n  Causal edge quality breakdown:")
causal_rels = ['CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF',
               'PROTECTS_AGAINST','ASSOCIATED_WITH']
for rel in causal_rels:
    count = stats['edges_by_relationship'].get(rel, 0)
    print(f"    {rel:30s}: {count:,}")

print(f"\n✓  Phase 1 query complete")
