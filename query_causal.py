"""
query_causal.py

Full causal reasoning demonstration for Disease-OS.
Shows root-cause tracing across all tiers, personalization,
and intervention impact analysis.

Run: python3 query_causal.py
"""

import sys, json
sys.path.insert(0, ".")
from core.causal_engine import CausalEngine
from core.config import DB_PATH

SEP = "─" * 60
engine = CausalEngine(str(DB_PATH))

print(f"\n{'═'*60}")
print(f"  Disease-OS Causal Reasoning Engine")
print(f"{'═'*60}")

def trace_disease(disease_id, disease_system, disease_name, patient_data=None):
    print(f"\n{SEP}")
    print(f"  DISEASE: {disease_name} ({disease_id})")
    print(SEP)

    paths = engine.trace(
        disease_id     = disease_id,
        disease_system = disease_system,
        max_depth      = 4,
        min_confidence = 0.50,
        min_path_score = 0.01,
        top_n          = 500,
        patient_data   = patient_data,
    )

    # Bucket by depth
    by_depth = {}
    for p in paths:
        by_depth.setdefault(p.depth, []).append(p)

    print(f"  Total paths: {len(paths)}")
    for d in sorted(by_depth):
        print(f"    Depth {d}: {len(by_depth[d])} paths")

    # Show tier origins
    tier_names = {
        1:"Molecular", 2:"Networks", 3:"Cellular",
        4:"Tissue", 5:"Systemic", 6:"Phenotype",
        7:"Disease", 8:"Behavior", 9:"Social", 10:"Healthcare"
    }
    tier_scores = {}
    for p in paths:
        t = p.edges[0].tier_from
        tier_scores[t] = tier_scores.get(t, 0) + p.path_score
    total = sum(tier_scores.values())

    print(f"\n  Root cause tier breakdown:")
    for tier in sorted(tier_scores, key=tier_scores.get, reverse=True):
        pct = tier_scores[tier] / total * 100 if total > 0 else 0
        name = tier_names.get(tier, f"Tier {tier}")
        bar  = "█" * int(pct / 3)
        print(f"    Tier {tier:2d} {name:10s}: {pct:5.1f}%  {bar}")

    # Show best path per tier
    best_per_tier = {}
    for p in paths:
        t = p.edges[0].tier_from
        if t not in best_per_tier or p.path_score > best_per_tier[t].path_score:
            best_per_tier[t] = p

    print(f"\n  Best path per root-cause tier:")
    for tier in sorted(best_per_tier, key=lambda t: best_per_tier[t].path_score, reverse=True):
        p    = best_per_tier[tier]
        root = p.edges[0]
        name = tier_names.get(tier, f"Tier {tier}")
        print(f"\n  Tier {tier} ({name})  score={p.path_score:.4f}  depth={p.depth}")
        for e in p.edges:
            print(f"    [{e.tier_from}] {e.source_label[:45]:<45} "
                  f"--[{e.rel_type}]--> [{e.tier_to}] {e.target_label[:35]}")

    # Deepest paths
    deepest = sorted(paths, key=lambda p: (-p.depth, -p.path_score))
    print(f"\n  Deepest causal chains (top 3):")
    for p in deepest[:3]:
        print(f"\n  depth={p.depth}  score={p.path_score:.4f}")
        for e in p.edges:
            print(f"    [{e.tier_from}] {e.source_label[:50]}")
            print(f"         --[{e.rel_type} conf={e.confidence:.2f}]-->")
        last = p.edges[-1]
        print(f"    [{last.tier_to}] {last.target_label[:50]}")

    return paths, tier_scores


# ── Type 2 Diabetes ────────────────────────────────────────────────────────
t2d_paths, t2d_tiers = trace_disease("E11.9", "ICD-10-CM", "Type 2 Diabetes")

# ── Personalized T2D — patient with rs7903146 + poor sleep ────────────────
print(f"\n{SEP}")
print(f"  PERSONALIZED: T2D patient with rs7903146 + poor sleep")
print(SEP)

patient = {
    "variants":    ["rs7903146"],
    "metabolites": ["HMDB0000122"],
}

pers_paths = engine.trace(
    disease_id="E11.9", disease_system="ICD-10-CM",
    max_depth=4, min_confidence=0.50,
    min_path_score=0.01, top_n=500,
    patient_data=patient,
)

# Find rs7903146 path and the sleep path
rs_path   = next((p for p in pers_paths
                  if any(e.source_id=="rs7903146" for e in p.edges)), None)
sleep_path = next((p for p in pers_paths
                   if any(e.source_id=="POOR_SLEEP" for e in p.edges)
                   and p.depth > 1), None)

if rs_path:
    print(f"\n  Genetic path (boosted by patient variant):")
    print(f"  score={rs_path.path_score:.4f}  depth={rs_path.depth}")
    for e in rs_path.edges:
        print(f"    {e.source_label[:45]} --[{e.rel_type}]--> {e.target_label[:35]}")

if sleep_path:
    print(f"\n  Behavioral path (patient has poor sleep):")
    print(f"  score={sleep_path.path_score:.4f}  depth={sleep_path.depth}")
    for e in sleep_path.edges:
        print(f"    {e.source_label[:45]} --[{e.rel_type}]--> {e.target_label[:35]}")

# ── Alzheimer's ────────────────────────────────────────────────────────────
alz = engine.conn.execute("""
    SELECT primary_id, primary_system, label FROM nodes
    WHERE label LIKE '%lzheimer%' AND entity_type='Disease_clinical'
    ORDER BY length(label) LIMIT 1
""").fetchone()

if alz:
    alz_paths, _ = trace_disease(alz[0], alz[1], f"Alzheimer's Disease")

# ── Hypertension ───────────────────────────────────────────────────────────
hyp_paths, _ = trace_disease("I10", "ICD-10-CM", "Hypertension")

# ── Depression ─────────────────────────────────────────────────────────────
dep_paths, _ = trace_disease("F32.9", "ICD-10-CM", "Major Depression")

# ── Export ─────────────────────────────────────────────────────────────────
print(f"\n{SEP}")
print(f"  EXPORT")
print(SEP)

output = {
    "t2d": {
        "paths": len(t2d_paths),
        "top_paths": [p.to_dict() for p in
                      sorted(t2d_paths,
                             key=lambda p: (-p.depth, -p.path_score))[:15]],
    }
}

with open("data/processed/causal_traces.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"  Saved to data/processed/causal_traces.json")
print(f"\n{'═'*60}")
print(f"  Causal reasoning complete")
print(f"{'═'*60}\n")
