"""
layers/tier7_disease/load_open_targets.py

Loads Open Targets 26.06 into Disease-OS.

Run from project root:
    python3 layers/tier7_disease/load_open_targets.py
"""

import sys
from pathlib import Path
from datetime import datetime

sys.path.insert(0, ".")
from core.graph_store import GraphStore
from core.node        import Node
from core.edge        import Edge
from core.config      import DB_PATH
from core.sources.open_targets import OpenTargetsSource

BATCH_SIZE = 5_000


def bulk_insert_nodes(conn, nodes):
    if not nodes:
        return
    sample = Node("x", "x", "x", 0, "Gene")
    cols   = list(sample.to_dict().keys())
    ph     = ", ".join(["?"] * len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO nodes ({', '.join(cols)}) VALUES ({ph})",
            [list(n.to_dict().values()) for n in nodes]
        )
        conn.execute("COMMIT")
    except Exception as e:
        conn.execute("ROLLBACK")
        raise e


def bulk_insert_edges(conn, edges):
    if not edges:
        return
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
    except Exception as e:
        conn.execute("ROLLBACK")
        raise e


def validate(gs):
    conn = gs.conn
    print("\n── Validation ───────────────────────────────────────")

    # OT edges loaded
    ot_edges = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE imported_via LIKE 'OpenTargets%'
    """).fetchone()[0]
    print(f"  Open Targets edges    : {ot_edges:,}")

    # High-confidence hits
    high_conf = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE imported_via LIKE 'OpenTargets%'
          AND confidence >= 0.50
    """).fetchone()[0]
    print(f"  High-confidence (≥0.5): {high_conf:,}")

    # T2D targets
    t2d_targets = conn.execute("""
        SELECT e.source_id, e.confidence
        FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id
                     AND n.primary_system=e.target_system
        WHERE n.label LIKE '%type 2 diabetes%'
          AND e.imported_via LIKE 'OpenTargets%'
          AND e.confidence >= 0.30
        ORDER BY e.confidence DESC
        LIMIT 8
    """).fetchall()
    print(f"\n  Top T2D targets from Open Targets:")
    for ensembl_id, conf in t2d_targets:
        print(f"    {ensembl_id:<25} score={conf:.3f}")

    # Alzheimer's targets
    alz_count = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id
                     AND n.primary_system=e.target_system
        WHERE n.label LIKE '%lzheimer%'
          AND e.imported_via LIKE 'OpenTargets%'
    """).fetchone()[0]
    print(f"\n  Alzheimer's gene targets  : {alz_count:,}")

    # Full stats
    s = gs.stats()
    print(f"\n── Graph stats ──────────────────────────────────────")
    print(f"  Total nodes : {s['total_nodes']:,}")
    print(f"  Total edges : {s['total_edges']:,}")
    print(f"  ASSOCIATED_WITH edges: {s['edges_by_relationship'].get('ASSOCIATED_WITH',0):,}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading Open Targets 26.06")
    print(f"  DB: {DB_PATH}")

    gs     = GraphStore(DB_PATH)
    conn   = gs.conn
    source = OpenTargetsSource()

    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    # Step 1 — disease nodes
    print(f"\n[1/2] Loading OT disease nodes...")
    node_batch = []
    node_total = 0

    for node in source.nodes():
        node_batch.append(node)
        if len(node_batch) >= BATCH_SIZE:
            bulk_insert_nodes(conn, node_batch)
            node_total += len(node_batch)
            node_batch  = []
            if node_total % 10_000 == 0:
                print(f"  {node_total:,} nodes...")

    if node_batch:
        bulk_insert_nodes(conn, node_batch)
        node_total += len(node_batch)

    t1 = datetime.now()
    print(f"  {node_total:,} disease nodes in {int((t1-t0).total_seconds())}s")

    # Step 2 — association edges
    print(f"\n[2/2] Loading OT association edges...")
    edge_batch = []
    edge_total = 0

    for edge in source.edges():
        if edge is None:
            continue
        edge_batch.append(edge)
        if len(edge_batch) >= BATCH_SIZE:
            bulk_insert_edges(conn, edge_batch)
            edge_total += len(edge_batch)
            edge_batch  = []
            if edge_total % 50_000 == 0:
                print(f"  {edge_total:,} edges...")

    if edge_batch:
        bulk_insert_edges(conn, edge_batch)
        edge_total += len(edge_batch)

    t2 = datetime.now()
    print(f"  {edge_total:,} edges in {int((t2-t1).total_seconds())}s")

    conn.execute("PRAGMA synchronous=NORMAL;")
    print(f"\n[Disease-OS] Total: {int((t2-t0).total_seconds())}s")

    validate(gs)
    print(f"\n✓  Open Targets loaded into {DB_PATH}")
