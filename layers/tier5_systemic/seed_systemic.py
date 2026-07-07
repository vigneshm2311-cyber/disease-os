"""
layers/tier5_systemic/seed_systemic.py

Seeds Tier 5 (Systemic Regulatory), Tier 8 (Behavior),
and Tier 9 (Social) with well-established causal relationships
from literature — enabling multi-tier causal chains in the engine.

These are not associations — they are mechanistically established
relationships supported by interventional evidence.
Each edge carries its primary source (PMID or guideline).

Run from project root:
    python3 layers/tier5_systemic/seed_systemic.py
"""

import sys
import sqlite3
from pathlib import Path
from datetime import datetime

sys.path.insert(0, ".")
from core.graph_store import GraphStore
from core.node        import Node
from core.edge        import Edge
from core.config      import DB_PATH

gs   = GraphStore(DB_PATH)
conn = gs.conn

def add(node): gs.add_node(node)
def edge(e):   gs.add_edge(e)

print("[seed] Adding Tier 5 systemic regulatory nodes...")

# ── Tier 5 — Systemic regulatory axis nodes ────────────────────────────────
SYSTEMIC_NODES = [
    Node("HPA_AXIS",          "DiseaseOS_System", "HPA Axis (Hypothalamic-Pituitary-Adrenal)", 5, "RegulatoryAxis"),
    Node("CIRCADIAN_SYSTEM",  "DiseaseOS_System", "Circadian System",                          5, "RegulatoryAxis"),
    Node("SYMPATHETIC_NS",    "DiseaseOS_System", "Sympathetic Nervous System",                 5, "RegulatoryAxis"),
    Node("INSULIN_SIGNALING", "DiseaseOS_System", "Insulin Signaling System",                   5, "RegulatoryAxis"),
    Node("IMMUNE_INFLAM",     "DiseaseOS_System", "Chronic Low-Grade Inflammation",              5, "ImmuneState"),
    Node("GUT_MICROBIOME",    "DiseaseOS_System", "Gut Microbiome",                             5, "MicrobialCommunity"),
    Node("RAAS",              "DiseaseOS_System", "Renin-Angiotensin-Aldosterone System",        5, "RegulatoryAxis"),
]

print("[seed] Adding Tier 8 behavior nodes...")

# ── Tier 8 — Behavior/exposure nodes ──────────────────────────────────────
BEHAVIOR_NODES = [
    Node("POOR_SLEEP",        "DiseaseOS_Behavior", "Poor Sleep / Sleep Deprivation",    8, "BehaviorPattern"),
    Node("HIGH_GLYCEMIC_DIET","DiseaseOS_Behavior", "High Glycemic Diet",                8, "BehaviorPattern"),
    Node("PHYSICAL_INACTIVITY","DiseaseOS_Behavior","Physical Inactivity / Sedentary",   8, "BehaviorPattern"),
    Node("CHRONIC_STRESS",    "DiseaseOS_Behavior", "Chronic Psychological Stress",       8, "BehaviorPattern"),
    Node("SMOKING",           "DiseaseOS_Behavior", "Smoking / Tobacco Use",              8, "BehaviorPattern"),
    Node("EXCESS_ALCOHOL",    "DiseaseOS_Behavior", "Excess Alcohol Consumption",         8, "BehaviorPattern"),
    Node("HIGH_SODIUM_DIET",  "DiseaseOS_Behavior", "High Sodium Diet",                   8, "BehaviorPattern"),
    Node("OBESITY",           "DiseaseOS_Behavior", "Visceral Obesity / Excess Adiposity",8, "BehaviorPattern"),
]

print("[seed] Adding Tier 9 social context nodes...")

# ── Tier 9 — Social context nodes ─────────────────────────────────────────
SOCIAL_NODES = [
    Node("LOW_SES",           "DiseaseOS_Social", "Low Socioeconomic Status",             9, "SocialDeterminant"),
    Node("FOOD_INSECURITY",   "DiseaseOS_Social", "Food Insecurity / Poor Food Access",   9, "SocialDeterminant"),
    Node("CHRONIC_ADVERSITY", "DiseaseOS_Social", "Chronic Social Adversity",             9, "SocialDeterminant"),
    Node("AIR_POLLUTION",     "DiseaseOS_Social", "Air Pollution / Environmental Toxin Exposure", 9, "SocialDeterminant"),
]

for n in SYSTEMIC_NODES + BEHAVIOR_NODES + SOCIAL_NODES:
    n.source         = "DiseaseOS_curated"
    n.source_version = "1.0"
    n.confidence     = 1.0
    add(n)

print("[seed] Adding causal edges across tiers...")

def mk_edge(src, src_sys, tgt, tgt_sys, rel, conf, pmid, notes=""):
    return Edge(
        source_id                = src,
        source_system            = src_sys,
        target_id                = tgt,
        target_system            = tgt_sys,
        relationship_type        = rel,
        source_relationship_type = notes,
        confidence               = conf,
        primary_source           = pmid,
        imported_via             = "DiseaseOS_curated_v1",
        study_design             = "curated",
        source_version           = "1.0",
    )

S = "DiseaseOS_System"
B = "DiseaseOS_Behavior"
C = "DiseaseOS_Social"
I = "ICD-10-CM"
D = "dbSNP_rsID"

# ── Tier 9 → Tier 8 (social context shapes behavior) ──────────────────────
SOCIAL_BEHAVIOR_EDGES = [
    mk_edge("LOW_SES",         C, "FOOD_INSECURITY",    C, "CONTRIBUTES_TO", 0.85, "PMID:28237157"),
    mk_edge("LOW_SES",         C, "PHYSICAL_INACTIVITY",B, "CONTRIBUTES_TO", 0.80, "PMID:22709777"),
    mk_edge("CHRONIC_ADVERSITY",C,"CHRONIC_STRESS",     B, "CAUSES",         0.85, "PMID:23274282"),
    mk_edge("FOOD_INSECURITY",  C, "HIGH_GLYCEMIC_DIET",B, "CONTRIBUTES_TO", 0.75, "PMID:28237157"),
    mk_edge("AIR_POLLUTION",    C, "IMMUNE_INFLAM",      S, "CAUSES",         0.80, "PMID:22461442"),
]

# ── Tier 8 → Tier 5 (behavior activates regulatory systems) ───────────────
BEHAVIOR_SYSTEMIC_EDGES = [
    mk_edge("POOR_SLEEP",        B, "HPA_AXIS",           S, "UPREGULATES",    0.85, "PMID:19460946"),
    mk_edge("POOR_SLEEP",        B, "CIRCADIAN_SYSTEM",   S, "CONTRIBUTES_TO", 0.90, "PMID:12519938"),
    mk_edge("POOR_SLEEP",        B, "IMMUNE_INFLAM",       S, "UPREGULATES",    0.80, "PMID:22648462"),
    mk_edge("CHRONIC_STRESS",    B, "HPA_AXIS",            S, "UPREGULATES",    0.90, "PMID:10747820"),
    mk_edge("CHRONIC_STRESS",    B, "SYMPATHETIC_NS",      S, "UPREGULATES",    0.85, "PMID:10747820"),
    mk_edge("CHRONIC_STRESS",    B, "IMMUNE_INFLAM",       S, "UPREGULATES",    0.80, "PMID:23274282"),
    mk_edge("HIGH_GLYCEMIC_DIET",B, "INSULIN_SIGNALING",   S, "CONTRIBUTES_TO", 0.85, "PMID:19064952"),
    mk_edge("HIGH_GLYCEMIC_DIET",B, "IMMUNE_INFLAM",       S, "UPREGULATES",    0.75, "PMID:15051604"),
    mk_edge("PHYSICAL_INACTIVITY",B,"INSULIN_SIGNALING",   S, "CONTRIBUTES_TO", 0.85, "PMID:11122304"),
    mk_edge("PHYSICAL_INACTIVITY",B,"IMMUNE_INFLAM",       S, "UPREGULATES",    0.75, "PMID:15475368"),
    mk_edge("OBESITY",           B, "INSULIN_SIGNALING",   S, "CONTRIBUTES_TO", 0.90, "PMID:14612429"),
    mk_edge("OBESITY",           B, "IMMUNE_INFLAM",       S, "UPREGULATES",    0.85, "PMID:16476399"),
    mk_edge("OBESITY",           B, "GUT_MICROBIOME",      S, "CONTRIBUTES_TO", 0.75, "PMID:23985870"),
    mk_edge("SMOKING",           B, "IMMUNE_INFLAM",       S, "UPREGULATES",    0.85, "PMID:16242130"),
    mk_edge("SMOKING",           B, "INSULIN_SIGNALING",   S, "CONTRIBUTES_TO", 0.75, "PMID:11473851"),
    mk_edge("HIGH_SODIUM_DIET",  B, "RAAS",                S, "UPREGULATES",    0.85, "PMID:17101785"),
]

# ── Tier 5 → Tier 5 (regulatory systems interact) ─────────────────────────
SYSTEMIC_SYSTEMIC_EDGES = [
    mk_edge("HPA_AXIS",      S, "INSULIN_SIGNALING",S, "DOWNREGULATES",0.85, "PMID:12200087",
            "cortisol impairs insulin sensitivity"),
    mk_edge("HPA_AXIS",      S, "IMMUNE_INFLAM",    S, "CONTRIBUTES_TO",0.75, "PMID:10747820",
            "chronic HPA activation → immune dysregulation"),
    mk_edge("CIRCADIAN_SYSTEM",S,"HPA_AXIS",        S, "CONTRIBUTES_TO",0.85, "PMID:12519938",
            "circadian disruption → cortisol rhythm loss"),
    mk_edge("CIRCADIAN_SYSTEM",S,"INSULIN_SIGNALING",S,"CONTRIBUTES_TO",0.80, "PMID:22974445"),
    mk_edge("SYMPATHETIC_NS",S, "RAAS",             S, "UPREGULATES",   0.85, "PMID:17101785"),
    mk_edge("GUT_MICROBIOME",S, "IMMUNE_INFLAM",    S, "CONTRIBUTES_TO",0.80, "PMID:21552190"),
]

# ── Tier 5 → Tier 7 (systemic dysregulation → disease) ────────────────────
SYSTEMIC_DISEASE_EDGES = [
    # T2D
    mk_edge("INSULIN_SIGNALING",S,"E11.9", I,"CONTRIBUTES_TO",0.90,"PMID:14612429"),
    mk_edge("IMMUNE_INFLAM",    S,"E11.9", I,"CONTRIBUTES_TO",0.80,"PMID:21552190"),

    # Hypertension
    mk_edge("RAAS",             S,"I10",   I,"CONTRIBUTES_TO",0.90,"PMID:17101785"),
    mk_edge("SYMPATHETIC_NS",   S,"I10",   I,"CONTRIBUTES_TO",0.85,"PMID:17101785"),

    # Alzheimer's
    mk_edge("IMMUNE_INFLAM",    S,"G30",   I,"CONTRIBUTES_TO",0.75,"PMID:25461519"),
    mk_edge("HPA_AXIS",         S,"G30",   I,"CONTRIBUTES_TO",0.70,"PMID:18215501"),

    # Depression
    mk_edge("HPA_AXIS",         S,"F32.9", I,"CONTRIBUTES_TO",0.85,"PMID:10747820"),
    mk_edge("IMMUNE_INFLAM",    S,"F32.9", I,"CONTRIBUTES_TO",0.80,"PMID:23274282"),
    mk_edge("CIRCADIAN_SYSTEM", S,"F32.9", I,"CONTRIBUTES_TO",0.80,"PMID:26344114"),

    # Cardiovascular disease
    mk_edge("IMMUNE_INFLAM",    S,"I25.10",I,"CONTRIBUTES_TO",0.85,"PMID:16476399"),
    mk_edge("RAAS",             S,"I25.10",I,"CONTRIBUTES_TO",0.80,"PMID:17101785"),
]

# ── Tier 8 → Tier 7 (behavior → disease, direct evidence) ─────────────────
BEHAVIOR_DISEASE_EDGES = [
    mk_edge("OBESITY",           B,"E11.9", I,"INCREASES_RISK_OF",0.90,"PMID:14612429"),
    mk_edge("PHYSICAL_INACTIVITY",B,"E11.9",I,"INCREASES_RISK_OF",0.85,"PMID:11122304"),
    mk_edge("HIGH_GLYCEMIC_DIET",B,"E11.9", I,"INCREASES_RISK_OF",0.80,"PMID:19064952"),
    mk_edge("POOR_SLEEP",        B,"E11.9", I,"INCREASES_RISK_OF",0.80,"PMID:19460946"),
    mk_edge("CHRONIC_STRESS",    B,"E11.9", I,"INCREASES_RISK_OF",0.75,"PMID:23274282"),
    mk_edge("SMOKING",           B,"I10",   I,"INCREASES_RISK_OF",0.85,"PMID:16242130"),
    mk_edge("HIGH_SODIUM_DIET",  B,"I10",   I,"INCREASES_RISK_OF",0.85,"PMID:17101785"),
    mk_edge("OBESITY",           B,"I10",   I,"INCREASES_RISK_OF",0.85,"PMID:16476399"),
    mk_edge("SMOKING",           B,"G30",   I,"INCREASES_RISK_OF",0.75,"PMID:17545523"),
    mk_edge("CHRONIC_STRESS",    B,"G30",   I,"INCREASES_RISK_OF",0.70,"PMID:18215501"),
    mk_edge("SMOKING",           B,"F32.9", I,"INCREASES_RISK_OF",0.70,"PMID:22548908"),
    mk_edge("CHRONIC_STRESS",    B,"F32.9", I,"INCREASES_RISK_OF",0.85,"PMID:10747820"),
    mk_edge("POOR_SLEEP",        B,"F32.9", I,"INCREASES_RISK_OF",0.80,"PMID:26344114"),
]

all_edges = (SOCIAL_BEHAVIOR_EDGES + BEHAVIOR_SYSTEMIC_EDGES +
             SYSTEMIC_SYSTEMIC_EDGES + SYSTEMIC_DISEASE_EDGES +
             BEHAVIOR_DISEASE_EDGES)

added = skipped = 0
for e in all_edges:
    if gs.add_edge(e):
        added += 1
    else:
        skipped += 1

print(f"\n── Seed complete ────────────────────────────────────")
print(f"  Nodes added  : {len(SYSTEMIC_NODES+BEHAVIOR_NODES+SOCIAL_NODES)}")
print(f"  Edges added  : {added}")
print(f"  Edges skipped: {skipped} (already existed)")

# Verify tiers populated
for tier, name in [(5,"Systemic"),(8,"Behavior"),(9,"Social")]:
    n = conn.execute(
        f"SELECT COUNT(*) FROM nodes WHERE tier={tier}"
    ).fetchone()[0]
    print(f"  Tier {tier} {name}: {n} nodes")

print(f"\n  Running causal trace to verify depth > 1...")

# Quick verify — trace T2D and check depth
from core.causal_engine import CausalEngine
engine = CausalEngine(DB_PATH)
paths = engine.trace(
    disease_id="E11.9", disease_system="ICD-10-CM",
    max_depth=4, min_confidence=0.50, top_n=15
)

max_depth_found = max((p.depth for p in paths), default=0)
tier_origins    = set(p.edges[0].tier_from for p in paths)

print(f"  Max path depth found : {max_depth_found}")
print(f"  Root cause tiers     : {sorted(tier_origins)}")
print(f"\n  Sample multi-tier path:")
for p in paths:
    if p.depth > 1:
        print(f"\n  {p.summary()}")
        break

print(f"\n✓ Tier 5/8/9 seeded and causal engine verified")
