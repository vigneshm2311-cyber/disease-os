"""
layers/tier1_molecular/load_hmdb.py

Loads HMDB v5.0 metabolites into Disease-OS Tier 1.

Run from project root:
    python3 layers/tier1_molecular/load_hmdb.py
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
from core.sources.hmdb import HMDBSource

BATCH_SIZE = 2_000


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

    # Total HMDB nodes
    total = conn.execute(
        "SELECT COUNT(*) FROM nodes WHERE source='HMDB'"
    ).fetchone()[0]
    print(f"  HMDB metabolite nodes : {total:,}")

    # Sample metabolites
    samples = conn.execute("""
        SELECT primary_id, label, confidence,
               json_extract(properties,'$.formula') as formula,
               json_extract(properties,'$.biofluids') as biofluids
        FROM nodes
        WHERE source = 'HMDB'
        ORDER BY confidence DESC
        LIMIT 5
    """).fetchall()
    print(f"\n  Sample metabolites:")
    for r in samples:
        print(f"    {r[0]}  {r[1][:40]:<40}  "
              f"formula={r[3]}  conf={r[2]}")
        if r[4]:
            import json
            bf = json.loads(r[4])
            if bf:
                print(f"      biofluids: {', '.join(bf[:3])}")

    # Biofluid edges
    biofluid_edges = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='DETECTED_IN'"
    ).fetchone()[0]
    print(f"\n  DETECTED_IN edges : {biofluid_edges:,}")

    # Disease association edges
    disease_edges = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE source_system='HMDB'
          AND relationship_type='ASSOCIATED_WITH'
    """).fetchone()[0]
    print(f"  Disease ASSOC edges : {disease_edges:,}")

    # T2D-associated metabolites
    t2d_mets = conn.execute("""
        SELECT COUNT(*) FROM edges
        WHERE source_system = 'HMDB'
          AND relationship_type = 'ASSOCIATED_WITH'
          AND (LOWER(target_id) LIKE '%diabetes%'
            OR LOWER(target_id) LIKE '%glucose%'
            OR LOWER(target_id) LIKE '%insulin%')
    """).fetchone()[0]
    print(f"  T2D-related metabolite edges: {t2d_mets:,}")

    # Common clinical metabolites
    for hmdb_id, name in [
        ("HMDB0000122", "Glucose"),
        ("HMDB0000243", "Pyruvate"),
        ("HMDB0000190", "Lactate"),
        ("HMDB0000695", "Cortisol"),
        ("HMDB0000517", "Creatinine"),
    ]:
        node = gs.get_node(hmdb_id, "HMDB")
        status = f"✅ {node['label']}" if node else "❌ not found"
        print(f"  {hmdb_id}: {status}")

    # Full stats
    print(f"\n── Graph stats ──────────────────────────────────────")
    s = gs.stats()
    print(f"  Total nodes : {s['total_nodes']:,}")
    print(f"  Total edges : {s['total_edges']:,}")
    print(f"  Tier 1 nodes: {s['nodes_by_tier'].get(1,0):,}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading HMDB v5.0 into Tier 1")
    print(f"  DB: {DB_PATH}")

    gs     = GraphStore(DB_PATH)
    conn   = gs.conn
    source = HMDBSource()

    # Step 1 — preprocess XML (skips if already done)
    source.preprocess()

    # Performance settings
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    # Step 2 — load metabolite nodes
    print(f"\n[1/2] Loading metabolite nodes...")
    node_batch = []
    node_total = 0

    for node in source.nodes():
        node_batch.append(node)
        if len(node_batch) >= BATCH_SIZE:
            bulk_insert_nodes(conn, node_batch)
            node_total += len(node_batch)
            node_batch  = []
            if node_total % 50_000 == 0:
                print(f"  {node_total:,} metabolites loaded...")

    if node_batch:
        bulk_insert_nodes(conn, node_batch)
        node_total += len(node_batch)

    t1 = datetime.now()
    print(f"  {node_total:,} metabolite nodes in {int((t1-t0).total_seconds())}s")

    # Step 3 — load edges
    print(f"\n[2/2] Loading metabolite edges...")
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
    print(f"\n✓  HMDB loaded into {DB_PATH}")
