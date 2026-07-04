"""
layers/tier2_networks/load_string.py

Loads STRING v12.0 protein-protein interactions into Disease-OS Tier 2.

Run from project root:
    python3 layers/tier2_networks/load_string.py
"""

import sys
import sqlite3
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.graph_store import GraphStore
from core.edge        import Edge
from core.config      import DB_PATH
from core.sources.string_db import StringSource

BATCH_SIZE = 10_000


def bulk_insert_edges(conn, edges: list):
    if not edges:
        return 0
    sample = Edge(
        source_id="x", source_system="x",
        target_id="y", target_system="y",
        relationship_type="ASSOCIATED_WITH",
        primary_source="x", imported_via="x",
        source_version="x",
    )
    cols = list(sample.to_dict().keys())
    ph   = ", ".join(["?"] * len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO edges ({', '.join(cols)}) VALUES ({ph})",
            [list(e.to_dict().values()) for e in edges]
        )
        conn.execute("COMMIT")
        return len(edges)
    except Exception as e:
        conn.execute("ROLLBACK")
        raise e


def validate(conn):
    print("\n── Validation ───────────────────────────────────────")

    # Count new INTERACTS_WITH edges
    total = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='INTERACTS_WITH'"
    ).fetchone()[0]
    print(f"  INTERACTS_WITH edges : {total:,}")

    # Sample high-confidence interactions
    samples = conn.execute("""
        SELECT source_id, target_id, confidence, effect_size
        FROM edges
        WHERE relationship_type = 'INTERACTS_WITH'
          AND confidence >= 0.90
        LIMIT 5
    """).fetchall()
    print(f"  High-confidence (≥0.90) samples:")
    for r in samples:
        print(f"    {r[0]} <-> {r[1]} "
              f"conf={r[2]} score={int(r[3]*1000)}")

    # Gene symbols that appear most (hub proteins)
    hubs = conn.execute("""
        SELECT source_id, COUNT(*) as n
        FROM edges
        WHERE relationship_type = 'INTERACTS_WITH'
          AND source_system = 'HGNC_Symbol'
        GROUP BY source_id
        ORDER BY n DESC
        LIMIT 8
    """).fetchall()
    print(f"\n  Most connected proteins (hubs):")
    for gene, n in hubs:
        print(f"    {gene:<15} {n:,} interactions")

    # Full stats
    print(f"\n── Graph stats ──────────────────────────────────────")
    total_nodes = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    total_edges = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    print(f"  Total nodes : {total_nodes:,}")
    print(f"  Total edges : {total_edges:,}")

    rel_counts = conn.execute("""
        SELECT relationship_type, COUNT(*) as n
        FROM edges
        GROUP BY relationship_type
        ORDER BY n DESC
    """).fetchall()
    for rel, n in rel_counts:
        print(f"  {rel:30s}: {n:,}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading STRING v12.0 into Tier 2")
    print(f"  DB: {DB_PATH}")

    gs   = GraphStore(DB_PATH)
    conn = gs.conn
    source = StringSource()

    # Performance settings
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    print(f"\n[1/1] Loading edges...")
    batch       = []
    total       = 0

    for edge in source.edges():
        if edge is None:
            continue
        batch.append(edge)
        if len(batch) >= BATCH_SIZE:
            bulk_insert_edges(conn, batch)
            total += len(batch)
            batch  = []
            if total % 500_000 == 0:
                print(f"  {total:,} edges inserted...")

    if batch:
        bulk_insert_edges(conn, batch)
        total += len(batch)

    conn.execute("PRAGMA synchronous=NORMAL;")
    t1 = datetime.now()

    print(f"\n  {total:,} edges loaded in {int((t1-t0).total_seconds())}s")

    validate(conn)

    # Commit
    gs.conn.commit()

    print(f"\n✓  STRING loaded into {DB_PATH}")
