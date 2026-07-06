"""
query_causal.py

Demonstrates the CausalEngine on Type 2 Diabetes.
Shows root-cause tracing, intervention impact, and
personalized path scoring for a patient with TCF7L2 variant.

Run from project root:
    python3 query_causal.py
"""

import sys, json
sys.path.insert(0, ".")

from core.causal_engine import CausalEngine
from core.config import DB_PATH

SEP = "─" * 60
engine = CausalEngine(DB_PATH)

print(f"\n{'═'*60}")
print(f"  Disease-OS Causal Reasoning Engine")
print(f"  Disease: Type 2 Diabetes (E11.9, ICD-10-CM)")
print(f"{'═'*60}")

# ── 1. Basic root cause trace ──────────────────────────────────────────────
print(f"\n{SEP}")
print(f"  1. ROOT CAUSE TRACE — top upstream contributors")
print(SEP)

paths = engine.trace(
    disease_id     = "E11.9",
    disease_system = "ICD-10-CM",
    max_depth      = 3,
    min_confidence = 0.55,
    top_n          = 20,
)

print(f"  Found {len(paths)} causal paths (showing top 10)\n")

for i, path in enumerate(paths[:10]):
    root = path.edges[0]
    print(f"  [{i+1}] score={path.path_score:.4f}  depth={path.depth}")
    print(f"      Root: [{root.tier_from}] {root.source_label[:55]}")
    print(f"      Via : {' → '.join(e.rel_type for e in path.edges)}")
    print(f"      Conf: {' × '.join(f'{e.confidence:.2f}' for e in path.edges)}")
    print()

# ── 2. Tier breakdown of root causes ──────────────────────────────────────
print(f"\n{SEP}")
print(f"  2. ROOT CAUSE TIERS — where do causes originate?")
print(SEP)

from collections import Counter
tier_names = {
    1:"Molecular", 2:"Networks", 3:"Cellular",
    4:"Tissue/Organ", 5:"Systemic", 6:"Phenotype",
    7:"Disease", 8:"Behavior", 9:"Social", 10:"Healthcare"
}

tier_scores = {}
for path in paths:
    root_tier = path.edges[0].tier_from
    tier_scores[root_tier] = tier_scores.get(root_tier, 0) + path.path_score

total = sum(tier_scores.values())
for tier in sorted(tier_scores, key=tier_scores.get, reverse=True):
    pct = tier_scores[tier] / total * 100 if total > 0 else 0
    bar = "█" * int(pct / 2)
    name = tier_names.get(tier, f"Tier {tier}")
    print(f"  Tier {tier:2d} {name:12s}: {pct:5.1f}%  {bar}")

# ── 3. Top root cause nodes ────────────────────────────────────────────────
print(f"\n{SEP}")
print(f"  3. TOP ROOT CAUSE NODES — highest combined path scores")
print(SEP)

node_scores = {}
node_labels = {}
for path in paths:
    root = path.edges[0]
    key  = (root.source_id, root.source_system)
    node_scores[key] = node_scores.get(key, 0) + path.path_score
    node_labels[key] = root.source_label

top_nodes = sorted(node_scores.items(), key=lambda x: x[1], reverse=True)[:10]
for (nid, nsys), score in top_nodes:
    label = node_labels[(nid, nsys)]
    print(f"  {score:.4f}  [{nsys}] {nid:<20} {label[:45]}")

# ── 4. Personalized trace — patient with TCF7L2 variant ───────────────────
print(f"\n{SEP}")
print(f"  4. PERSONALIZED TRACE — patient with rs7903146 (TCF7L2)")
print(SEP)

patient = {
    "variants":    ["rs7903146"],   # T2D risk variant
    "metabolites": ["HMDB0000122"], # Elevated glucose
}

pers_paths = engine.trace(
    disease_id     = "E11.9",
    disease_system = "ICD-10-CM",
    max_depth      = 3,
    min_confidence = 0.55,
    top_n          = 20,
    patient_data   = patient,
)

print(f"  Patient data: variant rs7903146, elevated glucose")
print(f"  Found {len(pers_paths)} personalized paths\n")

print(f"  Top 5 paths boosted by patient data:")
for i, path in enumerate(pers_paths[:5]):
    root = path.edges[0]
    print(f"  [{i+1}] score={path.path_score:.4f}  "
          f"root=[{root.tier_from}] {root.source_label[:50]}")
    print(f"      path: {' → '.join(e.rel_type for e in path.edges)}")
    print()

# ── 5. Intervention impact ─────────────────────────────────────────────────
print(f"\n{SEP}")
print(f"  5. INTERVENTION IMPACT — what if we target specific nodes?")
print(SEP)

# Find the most connected upstream node from our top paths
if top_nodes:
    top_node_id = top_nodes[0][0][0]
    top_node_label = node_labels[top_nodes[0][0]]

    impact = engine.intervention_impact(
        disease_id     = "E11.9",
        disease_system = "ICD-10-CM",
        intervene_on   = top_node_id,
        max_depth      = 3,
    )

    print(f"  Intervening on: {top_node_label[:55]}")
    print(f"  Total paths found   : {impact['total_paths']}")
    print(f"  Paths through target: {impact['affected_paths']}")
    print(f"  Fraction of causal  ")
    print(f"  weight addressable  : {impact['fraction_covered']*100:.1f}%")

# ── 6. Alzheimer's trace ───────────────────────────────────────────────────
print(f"\n{SEP}")
print(f"  6. ALZHEIMER'S DISEASE — root cause trace")
print(SEP)

# Find Alzheimer's node
alz = engine.conn.execute("""
    SELECT primary_id, primary_system, label
    FROM nodes
    WHERE label LIKE '%lzheimer%'
      AND entity_type = 'Disease_clinical'
    ORDER BY length(label)
    LIMIT 1
""").fetchone()

if alz:
    print(f"  Anchor: {alz[2]} ({alz[0]}, {alz[1]})")
    alz_paths = engine.trace(
        disease_id     = alz[0],
        disease_system = alz[1],
        max_depth      = 3,
        min_confidence = 0.55,
        top_n          = 10,
    )
    print(f"  Found {len(alz_paths)} causal paths\n")
    for i, path in enumerate(alz_paths[:5]):
        root = path.edges[0]
        print(f"  [{i+1}] score={path.path_score:.4f}")
        print(f"      Root: [{root.tier_from}] {root.source_label[:55]}")
        print(f"      Via : {' → '.join(e.rel_type for e in path.edges)}")
        print()
else:
    print("  Alzheimer's node not found — check disease loading")

# ── 7. Export top paths as JSON ────────────────────────────────────────────
print(f"\n{SEP}")
print(f"  7. EXPORT — saving top T2D paths to JSON")
print(SEP)

output = {
    "disease":       "Type 2 Diabetes",
    "disease_id":    "E11.9",
    "system":        "ICD-10-CM",
    "n_paths":       len(paths),
    "top_paths":     [p.to_dict() for p in paths[:10]],
    "tier_breakdown": {
        tier_names.get(t, f"Tier {t}"): round(s, 4)
        for t, s in sorted(tier_scores.items())
    },
}

out_path = "data/processed/t2d_causal_trace.json"
with open(out_path, "w") as f:
    json.dump(output, f, indent=2)

print(f"  Saved to {out_path}")
print(f"\n{'═'*60}")
print(f"  Causal reasoning complete")
print(f"{'═'*60}\n")
