"""
layers/tier1_molecular/load_uniprot.py

Loads UniProt Swiss-Prot human proteins into Disease-OS.

Run from project root:
    python3 layers/tier1_molecular/load_uniprot.py
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
from core.sources.uniprot import UniProtSource

BATCH_SIZE = 2_000


def bulk_nodes(conn, nodes):
    if not nodes:
        return
    sample = Node("x","x","x",0,"Gene")
    cols   = list(sample.to_dict().keys())
    ph     = ",".join(["?"]*len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO nodes ({','.join(cols)}) VALUES ({ph})",
            [list(n.to_dict().values()) for n in nodes]
        )
        conn.execute("COMMIT")
    except Exception as e:
        conn.execute("ROLLBACK")
        raise e


def bulk_edges(conn, edges):
    if not edges:
        return
    sample = Edge(
        source_id="x",source_system="x",
        target_id="y",target_system="y",
        relationship_type="ENCODES",
        primary_source="x",imported_via="x",source_version="x"
    )
    cols = list(sample.to_dict().keys())
    ph   = ",".join(["?"]*len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO edges ({','.join(cols)}) VALUES ({ph})",
            [list(e.to_dict().values()) for e in edges]
        )
        conn.execute("COMMIT")
    except Exception as e:
        conn.execute("ROLLBACK")
        raise e


def validate(gs, conn):
    print("\n── Validation ───────────────────────────────────────")

    # TP53 protein
    p53 = gs.get_node("P04637", "UniProt")
    if p53:
        print(f"  TP53 protein      : {p53['label'][:50]}")
        print(f"  Function preview  : {(p53.get('definition') or '')[:80]}")
        print(f"  NCBI Gene xref    : {p53['xrefs'].get('NCBI_Gene','n/a')}")

    # GLP1R protein (top T2D drug target)
    glp1r = conn.execute("""
        SELECT n.primary_id, n.label
        FROM nodes n
        WHERE n.primary_system = 'UniProt'
          AND n.label LIKE '%glucagon%peptide%receptor%'
        LIMIT 1
    """).fetchone()
    if glp1r:
        print(f"\n  GLP1R protein     : {glp1r[1][:55]}")

    # ENCODES edges
    encodes = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='ENCODES'"
    ).fetchone()[0]
    print(f"\n  ENCODES edges     : {encodes:,}")

    # HMDB pending resolution
    pending = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE target_system='UniProt_pending'"
    ).fetchone()[0]
    resolved_hmdb = conn.execute(
        "SELECT COUNT(*) FROM edges "
        "WHERE imported_via LIKE 'HMDB%' AND target_system='UniProt'"
    ).fetchone()[0]
    print(f"  HMDB→UniProt resolved: {resolved_hmdb:,}")
    print(f"  Still pending       : {pending:,}")

    # Full stats
    s = gs.stats()
    print(f"\n── Graph stats ──────────────────────────────────────")
    print(f"  Total nodes : {s['total_nodes']:,}")
    print(f"  Total edges : {s['total_edges']:,}")
    print(f"  Tier 1 nodes: {s['nodes_by_tier'].get(1,0):,}")
    print(f"\n  Edge types:")
    for rel, n in s['edges_by_relationship'].items():
        print(f"    {rel:<25}: {n:,}")


if __name__ == "__main__":
    print(f"[Disease-OS] Loading UniProt Swiss-Prot")
    print(f"  DB: {DB_PATH}")

    gs     = GraphStore(DB_PATH)
    conn   = gs.conn
    source = UniProtSource()

    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    # Step 1 — protein nodes
    print(f"\n[1/3] Loading protein nodes...")
    node_batch = []
    node_total = 0
    for node in source.nodes():
        node_batch.append(node)
        if len(node_batch) >= BATCH_SIZE:
            bulk_nodes(conn, node_batch)
            node_total += len(node_batch)
            node_batch  = []
    if node_batch:
        bulk_nodes(conn, node_batch)
        node_total += len(node_batch)
    t1 = datetime.now()
    print(f"  {node_total:,} protein nodes in {int((t1-t0).total_seconds())}s")

    # Step 2 — ENCODES + PART_OF edges
    print(f"\n[2/3] Loading edges...")
    edge_batch = []
    edge_total = 0
    for edge in source.edges():
        if edge is None:
            continue
        edge_batch.append(edge)
        if len(edge_batch) >= BATCH_SIZE:
            bulk_edges(conn, edge_batch)
            edge_total += len(edge_batch)
            edge_batch  = []
    if edge_batch:
        bulk_edges(conn, edge_batch)
        edge_total += len(edge_batch)
    t2 = datetime.now()
    print(f"  {edge_total:,} edges in {int((t2-t1).total_seconds())}s")

    # Step 3 — resolve staged HMDB edges
    print(f"\n[3/3] Resolving staged HMDB→UniProt edges...")
    hmdb_result = source.resolve_hmdb_pending(conn)

    conn.execute("PRAGMA synchronous=NORMAL;")
    t3 = datetime.now()
    print(f"\n[Disease-OS] Total: {int((t3-t0).total_seconds())}s")

    validate(gs, conn)
    print(f"\n✓  UniProt Swiss-Prot loaded into {DB_PATH}")
