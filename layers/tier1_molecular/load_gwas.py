"""
layers/tier1_molecular/load_gwas.py

Loads GWAS Catalog genome-wide significant associations.

Run from project root:
    python3 layers/tier1_molecular/load_gwas.py
"""

import sys
import sqlite3
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.graph_store import GraphStore
from core.node        import Node
from core.edge        import Edge
from core.config      import DB_PATH
from core.sources.gwas_catalog import GWASCatalogSource

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

    # Total GWAS edges
    gwas_edges = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE imported_via LIKE 'GWAS_Catalog%'
    """).fetchone()[0]
    print(f"  GWAS edges loaded     : {gwas_edges:,}")

    # T2D associations
    t2d = conn.execute("""
        SELECT source_id, effect_size, confidence,
               population_context, target_id
        FROM edges
        WHERE relationship_type = 'INCREASES_RISK_OF'
          AND target_id LIKE '%iabetes%'
        ORDER BY confidence DESC, effect_size DESC
        LIMIT 8
    """).fetchall()
    print(f"\n  T2D risk variants (top 8 by confidence):")
    for r in t2d:
        ef = f"OR={r[1]:.2f}" if r[1] else "OR=?"
        print(f"    {r[0]:<15} {ef:<12} conf={r[2]} pop={r[3]}")
        print(f"      trait: {r[4][:50]}")

    # Alzheimer's associations
    alz = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE relationship_type = 'INCREASES_RISK_OF'
          AND (target_id LIKE '%lzheimer%'
            OR target_id LIKE '%dementia%')
    """).fetchone()[0]
    print(f"\n  Alzheimer's/dementia risk variants : {alz:,}")

    # Population breakdown
    pops = conn.execute("""
        SELECT population_context, COUNT(*) as n
        FROM edges
        WHERE imported_via LIKE 'GWAS_Catalog%'
        GROUP BY population_context
        ORDER BY n DESC
        LIMIT 8
    """).fetchall()
    print(f"\n  Edges by ancestry population:")
    for pop, n in pops:
        print(f"    {pop:<30}: {n:,}")

    # Full stats
    print(f"\n── Graph stats ──────────────────────────────────────")
    s = gs.stats()
    print(f"  Total nodes : {s['total_nodes']:,}")
    print(f"  Total edges : {s['total_edges']:,}")
    for rel, n in s['edges_by_relationship'].items():
        print(f"  {rel:30s}: {n:,}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading GWAS Catalog")
    print(f"  DB: {DB_PATH}")

    gs     = GraphStore(DB_PATH)
    conn   = gs.conn
    source = GWASCatalogSource()

    # Step 1 — preprocess
    source.preprocess()

    # Performance
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    # Step 2 — load variant nodes
    print(f"\n[1/2] Loading GWAS variant nodes...")
    node_batch = []
    node_total = 0
    for node in source.nodes():
        node_batch.append(node)
        if len(node_batch) >= BATCH_SIZE:
            bulk_insert_nodes(conn, node_batch)
            node_total += len(node_batch)
            node_batch  = []
    if node_batch:
        bulk_insert_nodes(conn, node_batch)
        node_total += len(node_batch)
    t1 = datetime.now()
    print(f"  {node_total:,} variant nodes in {int((t1-t0).total_seconds())}s")

    # Step 3 — load association edges
    print(f"\n[2/2] Loading GWAS association edges...")
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
            if edge_total % 100_000 == 0:
                print(f"  {edge_total:,} edges loaded...")
    if edge_batch:
        bulk_insert_edges(conn, edge_batch)
        edge_total += len(edge_batch)

    t2 = datetime.now()
    print(f"  {edge_total:,} edges in {int((t2-t1).total_seconds())}s")

    conn.execute("PRAGMA synchronous=NORMAL;")
    print(f"\n[Disease-OS] Total: {int((t2-t0).total_seconds())}s")

    validate(gs)
    print(f"\n✓  GWAS Catalog loaded into {DB_PATH}")
