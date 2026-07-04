"""
query_t2d.py — End-to-end T2D causal trace after Phase 1 + resolver.
Run: python3 query_t2d.py
"""
import sys, json
sys.path.insert(0, ".")
from core.graph_store import GraphStore
from core.config import DB_PATH

gs   = GraphStore(DB_PATH)
conn = gs.conn
SEP  = "─" * 55

def section(title):
    print(f"\n{SEP}\n  {title}\n{SEP}")

# ── 1. Disease node ────────────────────────────────────────
section("1. DISEASE NODE — Type 2 Diabetes")
t2d = gs.find_by_icd10("E11.9")[0]
print(f"  Label      : {t2d['label']}")
print(f"  Primary ID : {t2d['primary_id']} ({t2d['primary_system']})")
print(f"  UMLS CUI   : {t2d['xrefs'].get('UMLS_CUI','n/a')}")
print(f"  SNOMED     : {t2d['xrefs'].get('SNOMED','n/a')}")
print(f"  HCC code   : {t2d['hcc_code']}  <- insurance risk adjustment")
print(f"  ICD-11     : {t2d['icd11_code']}  <- forward compatibility")

# ── 2. Causal variants — correct target IDs after resolver ─
section("2. CAUSAL VARIANTS → T2D CONCEPTS (conf ≥ 0.70)")

rows = conn.execute("""
    SELECT e.source_id,
           json_extract(n_src.properties,'$.gene_symbol') as gene,
           e.relationship_type, e.confidence,
           n_tgt.label as disease, n_tgt.primary_system
    FROM edges e
    JOIN nodes n_src ON n_src.primary_id = e.source_id
    JOIN nodes n_tgt ON n_tgt.primary_id   = e.target_id
                     AND n_tgt.primary_system = e.target_system
    WHERE n_src.entity_type = 'Variant'
      AND n_tgt.label LIKE '%iabetes%'
      AND e.relationship_type IN ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF')
      AND e.confidence >= 0.70
    ORDER BY e.confidence DESC, e.relationship_type
    LIMIT 15
""").fetchall()

print(f"  Found {len(rows)} high-confidence causal edges (showing 15)")
for r in rows:
    print(f"\n  [{r[2]}] conf={r[3]} gene={r[1]}")
    print(f"    Variant : {r[0]}")
    print(f"    Disease : {r[4][:55]} ({r[5]})")

# ── 3. All diabetes concepts with causal edges ─────────────
section("3. DIABETES CONCEPT MAP — all subtypes in graph")

concepts = conn.execute("""
    SELECT n.primary_id, n.primary_system, n.label,
           COUNT(DISTINCT e.source_id) as n_variants,
           SUM(CASE WHEN e.relationship_type='CAUSES' THEN 1 ELSE 0 END) as causes,
           SUM(CASE WHEN e.relationship_type='CONTRIBUTES_TO' THEN 1 ELSE 0 END) as contrib
    FROM nodes n
    JOIN edges e ON e.target_id = n.primary_id
                 AND e.target_system = n.primary_system
    WHERE n.label LIKE '%iabetes%'
      AND e.relationship_type IN ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF')
    GROUP BY n.primary_id, n.primary_system
    ORDER BY n_variants DESC
    LIMIT 12
""").fetchall()

print(f"  {'Label':<50} {'System':<12} {'Variants':>8} {'CAUSES':>7} {'CONTRIB':>7}")
print(f"  {'-'*50} {'-'*12} {'-'*8} {'-'*7} {'-'*7}")
for r in concepts:
    print(f"  {r[2][:50]:<50} {r[1]:<12} {r[3]:>8,} {r[4]:>7,} {r[5]:>7,}")

# ── 4. BRCA1 → Breast Cancer ───────────────────────────────
section("4. BRCA1 VARIANTS → Breast Cancer (after resolver)")

brca1 = conn.execute("""
    SELECT e.source_id, e.relationship_type, e.confidence,
           n_tgt.label, n_tgt.primary_system
    FROM edges e
    JOIN nodes n_src ON n_src.primary_id = e.source_id
    JOIN nodes n_tgt ON n_tgt.primary_id   = e.target_id
                     AND n_tgt.primary_system = e.target_system
    WHERE json_extract(n_src.properties,'$.gene_symbol') = 'BRCA1'
      AND e.relationship_type = 'CAUSES'
      AND e.confidence >= 0.85
    ORDER BY e.confidence DESC
    LIMIT 8
""").fetchall()

total_brca1 = conn.execute("""
    SELECT COUNT(*) FROM edges e
    JOIN nodes n ON n.primary_id = e.source_id
    WHERE json_extract(n.properties,'$.gene_symbol') = 'BRCA1'
      AND e.relationship_type = 'CAUSES'
""").fetchone()[0]

print(f"  Total BRCA1 CAUSES edges : {total_brca1:,}")
print(f"  Expert-panel reviewed (conf≥0.85):")
for r in brca1:
    print(f"    {r[0]} --[{r[1]}]--> {r[3][:45]} ({r[4]}) conf={r[2]}")

# ── 5. GLP-1 pathway → Drug connection ────────────────────
section("5. GLP-1 PATHWAY → DRUG TARGET CONNECTION")

glp1_pathway = conn.execute("""
    SELECT primary_id, label FROM nodes
    WHERE primary_id = 'R-HSA-381676'
""").fetchone()

if glp1_pathway:
    print(f"  Pathway: {glp1_pathway[1]}")
    members = conn.execute("""
        SELECT e.source_id, e.source_system, e.confidence
        FROM edges e
        WHERE e.target_id = 'R-HSA-381676'
          AND e.relationship_type = 'PART_OF'
        LIMIT 10
    """).fetchall()
    print(f"  Proteins in this pathway: {len(members)}")
    for m in members[:5]:
        print(f"    {m[0]} ({m[1]}) conf={m[2]}")

# ── 6. Insulin-related pathways ────────────────────────────
section("6. INSULIN/GLUCOSE PATHWAYS IN REACTOME")

pathways = conn.execute("""
    SELECT primary_id, label FROM nodes
    WHERE entity_type = 'Pathway'
      AND (label LIKE '%nsulin%'
        OR label LIKE '%lucose%'
        OR label LIKE '%lycoly%'
        OR label LIKE '%eta cell%'
        OR label LIKE '%GLP%')
    ORDER BY label
    LIMIT 12
""").fetchall()

for p in pathways:
    n_members = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE target_id=? AND relationship_type='PART_OF'",
        (p[0],)
    ).fetchone()[0]
    print(f"  {p[0]}  {p[1][:50]:<50}  ({n_members} members)")

# ── 7. Full causal summary ─────────────────────────────────
section("7. CAUSAL CHAIN — all edge types into diabetes concepts")

all_edges = conn.execute("""
    SELECT e.relationship_type, COUNT(*) as n
    FROM edges e
    JOIN nodes n ON n.primary_id   = e.target_id
                 AND n.primary_system = e.target_system
    WHERE n.label LIKE '%iabetes%'
    GROUP BY e.relationship_type
    ORDER BY n DESC
""").fetchall()

for rel, n in all_edges:
    print(f"  {rel:30s}: {n:,}")

# ── 8. Graph health ────────────────────────────────────────
section("8. GRAPH HEALTH")

stats = gs.stats()
tier_names = {
    1:"Molecular", 2:"Networks", 3:"Cellular",
    4:"Tissue/Organ", 5:"Systemic", 6:"Phenotype",
    7:"Disease", 8:"Behavior", 9:"Social",
    10:"Healthcare", 11:"Population"
}
print(f"  Total nodes: {stats['total_nodes']:,}  "
      f"Total edges: {stats['total_edges']:,}")
print(f"\n  Nodes by tier:")
for tier, count in sorted(stats['nodes_by_tier'].items()):
    bar  = "█" * min(count // 10000, 35)
    name = tier_names.get(tier,"?")
    print(f"    Tier {tier:2d} {name:12s}: {count:>9,}  {bar}")

print(f"\n  Causal edge quality:")
for rel in ['CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF',
            'PROTECTS_AGAINST','ASSOCIATED_WITH']:
    n = stats['edges_by_relationship'].get(rel, 0)
    print(f"    {rel:30s}: {n:,}")

print(f"\n✓  Phase 1 complete — ready for Phase 2")
