"""
core/gwas_trait_resolver.py

Resolves GWAS Catalog free-text trait labels to Disease-OS disease nodes.

Problem:
  GWAS edges have target_system='GWAS_trait_label' and
  target_id='Type 2 diabetes' — a free-text string.
  Disease nodes have labels like 'Type 2 diabetes mellitus
  without complications' stored under ICD-10-CM primary IDs.
  Direct JOIN fails.

Solution (three-pass matching):
  Pass 1 — Exact match: trait label == node label (case-insensitive)
  Pass 2 — Contains match: node label contains the trait label
  Pass 3 — Keyword match: all significant words in trait appear in node label

For each matched trait label, UPDATE all GWAS edges to point to
the matched node's primary_id + primary_system.

Unmatched traits remain with target_system='GWAS_trait_label'
and are flagged in the output for manual review.

Run from project root:
    python3 core/gwas_trait_resolver.py
"""

import sys
import re
import sqlite3
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))
from core.config import DB_PATH

# Words to ignore when doing keyword matching
STOPWORDS = {
    "the", "a", "an", "of", "in", "with", "and", "or", "for",
    "by", "to", "from", "at", "on", "is", "are", "was", "were",
    "type", "disease", "disorder", "syndrome", "condition",
    "snp", "interaction", "1df", "2df", "gwas",
}

def _keywords(text: str) -> set:
    """Extract significant keywords from a label."""
    words = re.sub(r"[^a-z0-9\s]", " ", text.lower()).split()
    return {w for w in words if w not in STOPWORDS and len(w) > 2}


def build_disease_index(conn) -> list:
    """
    Load all disease/phenotype/finding nodes into memory
    as a searchable index.
    Returns list of (primary_id, primary_system, label, label_lower, keywords)
    """
    print("[trait-resolver] Building disease node index...")
    rows = conn.execute("""
        SELECT primary_id, primary_system, label
        FROM nodes
        WHERE entity_type IN (
            'Disease_clinical', 'ClinicalFinding', 'LabTest', 'Phenotype'
        )
          AND label IS NOT NULL
          AND length(label) > 3
    """).fetchall()

    index = []
    for pid, psys, label in rows:
        index.append((
            pid, psys, label,
            label.lower().strip(),
            _keywords(label),
        ))

    print(f"[trait-resolver] {len(index):,} disease nodes indexed")
    return index


def resolve_trait(trait: str, index: list) -> tuple | None:
    """
    Try three passes to match a GWAS trait string to a disease node.
    Returns (primary_id, primary_system, label, match_type) or None.

    Cleans the trait first — removes SNP interaction qualifiers etc.
    """
    # Clean the trait label
    clean = re.sub(
        r"\(SNP x SNP interaction.*?\)|"
        r"\(adjusted.*?\)|"
        r"\(.*?df\)|"
        r"\(age.*?\)|"
        r"\(sex.*?\)",
        "", trait, flags=re.IGNORECASE
    ).strip().rstrip(",").strip()
    clean_lower = clean.lower().strip()

    if not clean_lower or len(clean_lower) < 4:
        return None

    # Pass 1 — exact match (case-insensitive)
    for pid, psys, label, label_lower, _ in index:
        if label_lower == clean_lower:
            return (pid, psys, label, "exact")

    # Pass 2 — node label contains the cleaned trait
    matches = []
    for pid, psys, label, label_lower, _ in index:
        if clean_lower in label_lower:
            # Prefer shorter labels (more specific match)
            matches.append((len(label), pid, psys, label))
    if matches:
        matches.sort()
        _, pid, psys, label = matches[0]
        return (pid, psys, label, "contains")

    # Pass 3 — keyword overlap (>= 80% of trait keywords in node label)
    trait_kws = _keywords(clean)
    if len(trait_kws) < 2:
        return None

    best_score = 0.0
    best_match = None
    for pid, psys, label, label_lower, node_kws in index:
        if not node_kws:
            continue
        overlap = len(trait_kws & node_kws)
        score   = overlap / len(trait_kws)
        if score >= 0.80 and score > best_score:
            best_score = score
            best_match = (pid, psys, label, "keyword")

    return best_match


def resolve_all(conn) -> dict:
    """
    Main resolution loop:
    1. Get distinct GWAS trait labels from edges
    2. Match each to a disease node
    3. UPDATE edges with matched primary_id + primary_system
    4. Report resolution rate
    """
    # Get all distinct unresolved GWAS trait labels
    traits = conn.execute("""
        SELECT target_id, COUNT(*) as n
        FROM edges
        WHERE target_system = 'GWAS_trait_label'
        GROUP BY target_id
        ORDER BY n DESC
    """).fetchall()

    total_traits  = len(traits)
    total_edges   = sum(r[1] for r in traits)
    print(f"\n[trait-resolver] {total_traits:,} distinct trait labels "
          f"covering {total_edges:,} edges")

    if total_traits == 0:
        print("[trait-resolver] Nothing to resolve.")
        return {}

    # Build disease index
    index = build_disease_index(conn)

    # Resolve each trait
    resolved   = {}   # trait_label -> (pid, psys, label, match_type)
    unresolved = []

    print(f"[trait-resolver] Resolving traits...")
    for trait, n_edges in traits:
        match = resolve_trait(trait, index)
        if match:
            resolved[trait] = match
        else:
            unresolved.append((trait, n_edges))

    print(f"  Resolved   : {len(resolved):,} / {total_traits:,} traits")
    print(f"  Unresolved : {len(unresolved):,} traits")

    # UPDATE edges for resolved traits
    print(f"\n[trait-resolver] Updating edges...")
    edges_updated = 0

    for trait, (pid, psys, label, match_type) in resolved.items():
        conn.execute("BEGIN")
        conn.execute("""
            UPDATE edges
            SET target_id     = ?,
                target_system = ?
            WHERE target_id     = ?
              AND target_system = 'GWAS_trait_label'
        """, (pid, psys, trait))
        n = conn.execute("SELECT changes()").fetchone()[0]
        conn.execute("COMMIT")
        edges_updated += n

    print(f"  {edges_updated:,} edges updated to primary node IDs")

    # Report unresolved traits (for manual review)
    remaining = conn.execute("""
        SELECT COUNT(*) FROM edges WHERE target_system = 'GWAS_trait_label'
    """).fetchone()[0]

    return {
        "total_traits":    total_traits,
        "resolved_traits": len(resolved),
        "unresolved_traits": len(unresolved),
        "edges_updated":   edges_updated,
        "edges_remaining": remaining,
        "unresolved_list": unresolved[:20],  # top 20 for review
    }


def validate(conn):
    print("\n── Validation ───────────────────────────────────────")

    # T2D risk variants now resolvable
    t2d = conn.execute("""
        SELECT e.source_id, e.effect_size, e.confidence,
               e.population_context, n.label
        FROM edges e
        JOIN nodes n ON n.primary_id   = e.target_id
                     AND n.primary_system = e.target_system
        WHERE e.relationship_type = 'INCREASES_RISK_OF'
          AND n.label LIKE '%iabetes%'
          AND e.effect_size IS NOT NULL
        ORDER BY e.confidence DESC, e.effect_size DESC
        LIMIT 8
    """).fetchall()

    print(f"  T2D INCREASES_RISK_OF edges resolved: {len(t2d)}")
    for r in t2d:
        ef = f"OR={r[1]:.2f}" if r[1] else "?"
        print(f"    {r[0]:<15} {ef:<10} conf={r[2]} pop={r[3]}")
        print(f"      -> {r[4][:55]}")

    # Alzheimer's
    alz = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN nodes n ON n.primary_id   = e.target_id
                     AND n.primary_system = e.target_system
        WHERE e.relationship_type = 'INCREASES_RISK_OF'
          AND n.label LIKE '%lzheimer%'
    """).fetchone()[0]
    print(f"\n  Alzheimer's risk edges resolved : {alz:,}")

    # Remaining unresolved
    remaining = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE target_system = 'GWAS_trait_label'
    """).fetchone()[0]
    print(f"  GWAS_trait_label edges remaining: {remaining:,}")


if __name__ == "__main__":
    print(f"[Disease-OS] GWAS Trait Resolver")
    print(f"  DB: {DB_PATH}")

    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")

    t0 = datetime.now()

    summary = resolve_all(conn)

    conn.execute("PRAGMA synchronous=NORMAL;")
    t1 = datetime.now()

    print(f"\n── Summary ──────────────────────────────────────────")
    print(f"  Traits resolved  : {summary.get('resolved_traits',0):,} / "
          f"{summary.get('total_traits',0):,}")
    print(f"  Edges updated    : {summary.get('edges_updated',0):,}")
    print(f"  Edges remaining  : {summary.get('edges_remaining',0):,}")
    print(f"  Time             : {int((t1-t0).total_seconds())}s")

    print(f"\n  Top unresolved traits (manual review):")
    for trait, n in summary.get("unresolved_list", [])[:15]:
        print(f"    [{n:4d} edges]  {trait}")

    validate(conn)

    print(f"\n✓  GWAS trait resolution complete")
    conn.close()
