"""
layers/tier1_molecular/load_clinvar.py

Loads ClinVar variant nodes and variant->disease edges into Disease-OS.

Run from project root:
    python3 layers/tier1_molecular/load_clinvar.py
"""

import sys
import csv
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node        import Node
from core.edge        import Edge
from core.graph_store import GraphStore
from core.config      import DB_PATH, PROCESSED_DIR
from core.sources.clinvar import ClinVarSource

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
    print("\n── Validation ───────────────────────────────────────")

    # Known pathogenic variant
    node = gs.get_node("rs7903146", "dbSNP_rsID")
    if node:
        print(f"  rs7903146         -> {node['label'][:60]}")
        print(f"  ClinSig           -> {node['properties'].get('clinsig','?')}")
        print(f"  Gene              -> {node['properties'].get('gene_symbol','?')}")
        print(f"  Confidence        -> {node['confidence']}")

    # BRCA1 variants
    conn = gs.conn
    brca1 = conn.execute(
        "SELECT COUNT(*) FROM nodes "
        "WHERE entity_type='Variant' "
        "AND json_extract(properties,'$.gene_symbol')='BRCA1'"
    ).fetchone()[0]
    print(f"\n  BRCA1 variants    -> {brca1:,} loaded")

    # CAUSES edges — the most valuable
    causes = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='CAUSES'"
    ).fetchone()[0]
    contrib = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='CONTRIBUTES_TO'"
    ).fetchone()[0]
    risk = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='INCREASES_RISK_OF'"
    ).fetchone()[0]
    print(f"\n  CAUSES edges      -> {causes:,}")
    print(f"  CONTRIBUTES_TO    -> {contrib:,}")
    print(f"  INCREASES_RISK_OF -> {risk:,}")

    # Sample causal chain
    edges = gs.edges_from("rs7903146", relationship_type="INCREASES_RISK_OF")
    if edges:
        e = edges[0]
        print(f"\n  Sample edge:")
        print(f"    {e['source_id']} --[{e['relationship_type']}]"
              f"--> {e['target_id']}")
        print(f"    confidence={e['confidence']} "
              f"study={e['study_design']}")

    # Full stats
    print(f"\n── Graph stats ──────────────────────────────────────")
    s = gs.stats()
    print(f"  Total nodes      : {s['total_nodes']:,}")
    print(f"  Total edges      : {s['total_edges']:,}")
    print(f"  Nodes by tier    : {s['nodes_by_tier']}")
    print(f"  Edges by rel     : {s['edges_by_relationship']}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading ClinVar into Tier 1")
    print(f"  DB: {DB_PATH}")

    gs     = GraphStore(DB_PATH)
    source = ClinVarSource()
    conn   = gs.conn

    # Step 1 — preprocess raw file (skips if already done)
    source.preprocess()

    # Step 2 — bulk load performance settings
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    # Step 3 — load variant nodes
    print(f"\n[1/2] Loading variant nodes...")
    node_batch = []
    node_total = 0

    for node in source.nodes():
        node_batch.append(node)
        if len(node_batch) >= BATCH_SIZE:
            bulk_insert_nodes(conn, node_batch)
            node_total += len(node_batch)
            node_batch  = []
            if node_total % 100_000 == 0:
                print(f"  {node_total:,} variants loaded...")

    if node_batch:
        bulk_insert_nodes(conn, node_batch)
        node_total += len(node_batch)

    t1 = datetime.now()
    print(f"  {node_total:,} variant nodes in {int((t1-t0).total_seconds())}s")

    # Step 4 — load variant->disease edges
    print(f"\n[2/2] Loading variant->disease edges...")
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
    print(f"\n✓  ClinVar loaded into {DB_PATH}")
