"""
core/resolver.py

Cross-reference resolver for Disease-OS.

Problem: ClinVar/UMLS edges point to UMLS_CUI targets/sources.
         Nodes are stored under domain-appropriate primary IDs.
         JOINs fail silently because IDs don't match.

Fix strategy (delete-then-reinsert):
  1. Build cui_map: CUI -> (primary_id, primary_system)
  2. For edges with UMLS_CUI source/target:
     a. Compute the resolved version in a temp table
     b. DELETE the original UMLS_CUI edges
     c. INSERT OR IGNORE the resolved versions
     — duplicates are silently dropped, no constraint errors

Safe to re-run — idempotent.

Run from project root:
    python3 core/resolver.py
"""

import sys
import sqlite3
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))
from core.config import DB_PATH

CUI_MAP_SCHEMA = """
CREATE TABLE IF NOT EXISTS cui_map (
    cui            TEXT PRIMARY KEY,
    primary_id     TEXT NOT NULL,
    primary_system TEXT NOT NULL,
    label          TEXT,
    tier           INTEGER
);
CREATE INDEX IF NOT EXISTS idx_cui_map_cui
    ON cui_map(cui);
"""


def build_cui_map(conn) -> int:
    """
    Scan all nodes, extract UMLS_CUI from xrefs, write to cui_map.
    Returns number of CUI mappings built.
    """
    print("[resolver] Building CUI map from node xrefs...")

    conn.execute("DELETE FROM cui_map")   # clear previous run
    conn.execute("""
        INSERT OR REPLACE INTO cui_map (cui, primary_id, primary_system, label, tier)
        SELECT
            json_extract(xrefs, '$.UMLS_CUI'),
            primary_id,
            primary_system,
            label,
            tier
        FROM nodes
        WHERE json_extract(xrefs, '$.UMLS_CUI') IS NOT NULL
    """)
    conn.commit()

    n = conn.execute("SELECT COUNT(*) FROM cui_map").fetchone()[0]
    print(f"[resolver] {n:,} CUI mappings built")
    return n


def resolve_targets(conn) -> dict:
    """
    Resolve edges where target_system = 'UMLS_CUI'.
    Strategy: build resolved set in temp table,
              delete originals, reinsert resolved with OR IGNORE.
    """
    print("\n[resolver] Resolving edge TARGETS...")

    # How many need resolving?
    total = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE target_system='UMLS_CUI'"
    ).fetchone()[0]
    print(f"  {total:,} edges with UMLS_CUI targets")

    if total == 0:
        return {"resolved": 0, "unresolved": 0, "duplicates_dropped": 0}

    resolvable = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN cui_map c ON c.cui = e.target_id
        WHERE e.target_system = 'UMLS_CUI'
    """).fetchone()[0]
    print(f"  {resolvable:,} resolvable  |  "
          f"{total - resolvable:,} CUIs not in graph (will be deleted)")

    # Step 1 — create temp table of resolved edges
    conn.execute("DROP TABLE IF EXISTS _resolved_edges")
    conn.execute("""
        CREATE TEMP TABLE _resolved_edges AS
        SELECT
            e.id,
            e.source_id,
            e.source_system,
            c.primary_id     AS target_id,
            c.primary_system AS target_system,
            e.relationship_type,
            e.source_relationship_type,
            e.effect_size,
            e.effect_unit,
            e.direction,
            e.confidence,
            e.feedback,
            e.feedback_notes,
            e.primary_source,
            e.imported_via,
            e.study_design,
            e.population_context,
            e.tissue_context,
            e.species,
            e.typical_latency,
            e.source_version,
            e.loaded_at
        FROM edges e
        JOIN cui_map c ON c.cui = e.target_id
        WHERE e.target_system = 'UMLS_CUI'
    """)

    resolved_count = conn.execute(
        "SELECT COUNT(*) FROM _resolved_edges"
    ).fetchone()[0]

    # Step 2 — delete ALL UMLS_CUI target edges (resolved + unresolved)
    conn.execute("BEGIN")
    conn.execute("DELETE FROM edges WHERE target_system = 'UMLS_CUI'")

    # Step 3 — reinsert resolved edges, silently drop duplicates
    conn.execute("""
        INSERT OR IGNORE INTO edges (
            source_id, source_system,
            target_id, target_system,
            relationship_type, source_relationship_type,
            effect_size, effect_unit, direction,
            confidence, feedback, feedback_notes,
            primary_source, imported_via, study_design,
            population_context, tissue_context, species,
            typical_latency, source_version, loaded_at
        )
        SELECT
            source_id, source_system,
            target_id, target_system,
            relationship_type, source_relationship_type,
            effect_size, effect_unit, direction,
            confidence, feedback, feedback_notes,
            primary_source, imported_via, study_design,
            population_context, tissue_context, species,
            typical_latency, source_version, loaded_at
        FROM _resolved_edges
    """)
    conn.execute("COMMIT")

    inserted = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE target_system != 'UMLS_CUI'
          AND loaded_at IN (SELECT loaded_at FROM _resolved_edges LIMIT 1)
    """).fetchone()[0]

    duplicates_dropped = resolved_count - conn.execute("""
        SELECT changes()
    """).fetchone()[0]

    # Remaining unresolved
    still_umls = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE target_system='UMLS_CUI'"
    ).fetchone()[0]

    print(f"  {resolved_count:,} resolved and reinserted")
    print(f"  {total - resolvable:,} unresolvable edges dropped")
    print(f"  {still_umls:,} UMLS_CUI target edges remaining")

    return {
        "resolved":          resolved_count,
        "unresolved":        total - resolvable,
        "still_umls":        still_umls,
    }


def resolve_sources(conn) -> dict:
    """
    Resolve edges where source_system = 'UMLS_CUI'.
    Same delete-then-reinsert pattern.
    """
    print("\n[resolver] Resolving edge SOURCES...")

    total = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE source_system='UMLS_CUI'"
    ).fetchone()[0]
    print(f"  {total:,} edges with UMLS_CUI sources")

    if total == 0:
        return {"resolved": 0, "unresolved": 0}

    resolvable = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN cui_map c ON c.cui = e.source_id
        WHERE e.source_system = 'UMLS_CUI'
    """).fetchone()[0]
    print(f"  {resolvable:,} resolvable  |  "
          f"{total - resolvable:,} not in graph (will be dropped)")

    conn.execute("DROP TABLE IF EXISTS _resolved_src")
    conn.execute("""
        CREATE TEMP TABLE _resolved_src AS
        SELECT
            e.id,
            c.primary_id     AS source_id,
            c.primary_system AS source_system,
            e.target_id,
            e.target_system,
            e.relationship_type,
            e.source_relationship_type,
            e.effect_size,
            e.effect_unit,
            e.direction,
            e.confidence,
            e.feedback,
            e.feedback_notes,
            e.primary_source,
            e.imported_via,
            e.study_design,
            e.population_context,
            e.tissue_context,
            e.species,
            e.typical_latency,
            e.source_version,
            e.loaded_at
        FROM edges e
        JOIN cui_map c ON c.cui = e.source_id
        WHERE e.source_system = 'UMLS_CUI'
    """)

    conn.execute("BEGIN")
    conn.execute("DELETE FROM edges WHERE source_system = 'UMLS_CUI'")
    conn.execute("""
        INSERT OR IGNORE INTO edges (
            source_id, source_system,
            target_id, target_system,
            relationship_type, source_relationship_type,
            effect_size, effect_unit, direction,
            confidence, feedback, feedback_notes,
            primary_source, imported_via, study_design,
            population_context, tissue_context, species,
            typical_latency, source_version, loaded_at
        )
        SELECT
            source_id, source_system,
            target_id, target_system,
            relationship_type, source_relationship_type,
            effect_size, effect_unit, direction,
            confidence, feedback, feedback_notes,
            primary_source, imported_via, study_design,
            population_context, tissue_context, species,
            typical_latency, source_version, loaded_at
        FROM _resolved_src
    """)
    conn.execute("COMMIT")

    still_umls = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE source_system='UMLS_CUI'"
    ).fetchone()[0]
    print(f"  {resolvable:,} resolved and reinserted")
    print(f"  {total - resolvable:,} unresolvable edges dropped")
    print(f"  {still_umls:,} UMLS_CUI source edges remaining")

    return {"resolved": resolvable, "unresolved": total - resolvable}


def validate(conn):
    print("\n── Validation ───────────────────────────────────────")

    # BRCA1 → Breast Cancer chain
    brca1 = conn.execute("""
        SELECT e.source_id, e.relationship_type,
               n.primary_id, n.primary_system, n.label
        FROM edges e
        JOIN nodes n ON n.primary_id   = e.target_id
                     AND n.primary_system = e.target_system
        WHERE e.source_id IN (
            SELECT primary_id FROM nodes
            WHERE entity_type = 'Variant'
              AND json_extract(properties,'$.gene_symbol') = 'BRCA1'
        )
          AND e.relationship_type = 'CAUSES'
        LIMIT 5
    """).fetchall()

    print(f"  BRCA1 CAUSES chains now resolvable: {len(brca1)}")
    for row in brca1[:3]:
        print(f"    {row[0]} --[{row[1]}]--> "
              f"{row[2]} ({row[3]})")
        print(f"    Target: {row[4][:55]}")

    # T2D edges
    t2d = conn.execute("""
        SELECT e.relationship_type, COUNT(*) as n
        FROM edges e
        JOIN nodes n ON n.primary_id   = e.target_id
                     AND n.primary_system = e.target_system
        WHERE n.label LIKE '%iabetes%'
          AND e.relationship_type IN
              ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF','TREATS')
        GROUP BY e.relationship_type
        ORDER BY n DESC
    """).fetchall()

    print(f"\n  Edges into diabetes concepts after resolution:")
    for rel, n in t2d:
        print(f"    {rel:25s} : {n:,}")

    # Full edge breakdown
    print(f"\n── Edge breakdown ───────────────────────────────────")
    for rel, n in conn.execute("""
        SELECT relationship_type, COUNT(*) as n
        FROM edges GROUP BY relationship_type ORDER BY n DESC
    """).fetchall():
        print(f"  {rel:30s}: {n:,}")

    # Total graph state
    total_nodes = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    total_edges = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    remaining   = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE source_system='UMLS_CUI' OR target_system='UMLS_CUI'
    """).fetchone()[0]

    print(f"\n── Graph state ──────────────────────────────────────")
    print(f"  Total nodes              : {total_nodes:,}")
    print(f"  Total edges              : {total_edges:,}")
    print(f"  UMLS_CUI edges remaining : {remaining:,}")


if __name__ == "__main__":
    print(f"[Disease-OS] Cross-reference resolver")
    print(f"  DB: {DB_PATH}")

    conn = sqlite3.connect(str(DB_PATH))

    for stmt in CUI_MAP_SCHEMA.strip().split(";"):
        if stmt.strip():
            conn.execute(stmt)
    conn.commit()

    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-512000;")   # 512MB — this is a heavy pass
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    n_mapped       = build_cui_map(conn)
    target_summary = resolve_targets(conn)
    source_summary = resolve_sources(conn)

    conn.execute("PRAGMA synchronous=NORMAL;")

    t1 = datetime.now()

    print(f"\n── Summary ──────────────────────────────────────────")
    print(f"  CUI mappings built    : {n_mapped:,}")
    print(f"  Target edges resolved : {target_summary['resolved']:,}")
    print(f"  Source edges resolved : {source_summary['resolved']:,}")
    print(f"  Time                  : {int((t1-t0).total_seconds())}s")

    validate(conn)
    print(f"\n✓  Resolution complete")
    conn.close()
