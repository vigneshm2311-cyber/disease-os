"""
layers/tier2_networks/load_reactome.py

Loads Reactome human pathways into Disease-OS Tier 2.

Run from project root:
    python3 layers/tier2_networks/load_reactome.py
"""

import sys
import sqlite3
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.graph_store import GraphStore
from core.node        import Node
from core.edge        import Edge
from core.config      import DB_PATH
from core.sources.reactome import ReactomeSource

BATCH_SIZE = 2_000


def bulk_insert_nodes(conn, nodes: list):
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


def bulk_insert_edges(conn, edges: list):
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


def load(gs: GraphStore):
    conn   = gs.conn
    source = ReactomeSource()

    # Performance settings
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-128000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    # ── Nodes ──────────────────────────────────────────────────────────
    print("\n[1/2] Loading Reactome pathway nodes...")
    t0          = datetime.now()
    node_batch  = []
    node_total  = 0

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
    print(f"  {node_total:,} pathway nodes loaded in {int((t1-t0).total_seconds())}s")

    # ── Edges ──────────────────────────────────────────────────────────
    print("\n[2/2] Loading Reactome edges...")
    edge_batch   = []
    edge_total   = 0
    edge_skipped = 0

    for edge in source.edges():
        if edge is None:
            continue
        edge_batch.append(edge)
        if len(edge_batch) >= BATCH_SIZE:
            bulk_insert_edges(conn, edge_batch)
            edge_total  += len(edge_batch)
            edge_batch   = []
            if edge_total % 50_000 == 0:
                print(f"  {edge_total:,} edges loaded...")

    if edge_batch:
        bulk_insert_edges(conn, edge_batch)
        edge_total += len(edge_batch)

    t2 = datetime.now()
    print(f"  {edge_total:,} edges loaded in {int((t2-t1).total_seconds())}s")

    # Restore safe mode
    conn.execute("PRAGMA synchronous=NORMAL;")

    return node_total, edge_total


def validate(gs: GraphStore):
    print("\n── Validation ───────────────────────────────────────")

    # Pathway node lookup
    node = gs.get_node("R-HSA-70171", "Reactome")
    if node:
        print(f"  Pathway R-HSA-70171 -> {node['label']}")
    else:
        results = gs.search_label("Glycolysis", limit=2)
        for r in results:
            if r["primary_system"] == "Reactome":
                print(f"  Glycolysis pathway  -> {r['label']} ({r['primary_id']})")
                break

    # Insulin signaling pathway
    results = gs.search_label("Insulin", limit=5)
    reactome_hits = [r for r in results if r["primary_system"] == "Reactome"]
    if reactome_hits:
        print(f"  Insulin pathways    -> {len(reactome_hits)} found")
        print(f"    e.g. {reactome_hits[0]['label']}")

    # Edges from a pathway
    edges = gs.edges_to("R-HSA-109582", relationship_type="PART_OF")
    print(f"  Sub-pathways of R-HSA-109582 (Hemostasis) -> {len(edges)} child pathways")

    # Full stats
    print(f"\n── Graph stats ──────────────────────────────────────")
    s = gs.stats()
    print(f"  Total nodes      : {s['total_nodes']:,}")
    print(f"  Total edges      : {s['total_edges']:,}")
    print(f"  Nodes by tier    : {s['nodes_by_tier']}")
    print(f"  Edges by rel     : {s['edges_by_relationship']}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading Reactome into Tier 2")
    print(f"  DB: {DB_PATH}")

    gs = GraphStore(DB_PATH)
    node_total, edge_total = load(gs)

    print(f"\n[Disease-OS] Reactome load complete")
    print(f"  Nodes: {node_total:,}  Edges: {edge_total:,}")

    validate(gs)

    print(f"\n✓  Reactome loaded into {DB_PATH}")
