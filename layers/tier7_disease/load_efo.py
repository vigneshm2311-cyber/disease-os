"""
layers/tier7_disease/load_efo.py

Loads EFO v3.91 nodes and resolves GWAS trait label edges.

Run from project root:
    python3 layers/tier7_disease/load_efo.py
"""

import sys, sqlite3
from pathlib import Path
from datetime import datetime

sys.path.insert(0, ".")
from core.graph_store import GraphStore
from core.node        import Node
from core.edge        import Edge
from core.config      import DB_PATH
from core.sources.efo import EFOSource

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


if __name__ == "__main__":
    print(f"[Disease-OS] Loading EFO v3.91")
    print(f"  DB: {DB_PATH}")

    gs     = GraphStore(DB_PATH)
    conn   = gs.conn
    source = EFOSource()

    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")

    t0 = datetime.now()

    # Step 1 — preprocess OBO
    source.preprocess()

    # Step 2 — load EFO nodes
    print(f"\n[1/2] Loading EFO/MONDO/HP nodes...")
    batch = []
    total = 0

    for node in source.nodes():
        batch.append(node)
        if len(batch) >= BATCH_SIZE:
            bulk_insert_nodes(conn, batch)
            total += len(batch)
            batch  = []
            if total % 20_000 == 0:
                print(f"  {total:,} nodes loaded...")

    if batch:
        bulk_insert_nodes(conn, batch)
        total += len(batch)

    t1 = datetime.now()
    print(f"  {total:,} EFO nodes in {int((t1-t0).total_seconds())}s")

    # Step 3 — resolve GWAS trait labels
    print(f"\n[2/2] Resolving GWAS trait label edges using EFO index...")
    result = source.resolve_gwas_edges(conn)

    conn.execute("PRAGMA synchronous=NORMAL;")
    t2 = datetime.now()

    print(f"\n── Results ──────────────────────────────────────────")
    print(f"  EFO nodes loaded     : {total:,}")
    print(f"  GWAS traits resolved : {result['traits_resolved']:,} / {result['traits_total']:,}")
    print(f"  Edges updated        : {result['edges_updated']:,}")
    print(f"  Still unresolved     : {result['remaining']:,}")
    print(f"  Total time           : {int((t2-t0).total_seconds())}s")

    # Validation
    print(f"\n── Validation ───────────────────────────────────────")

    t2d = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.relationship_type='INCREASES_RISK_OF'
          AND n.label LIKE '%iabetes%'
    """).fetchone()[0]
    print(f"  T2D risk edges resolved     : {t2d:,}")

    alz = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.relationship_type='INCREASES_RISK_OF'
          AND n.label LIKE '%lzheimer%'
    """).fetchone()[0]
    print(f"  Alzheimer risk edges        : {alz:,}")

    bmi = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.relationship_type='INCREASES_RISK_OF'
          AND n.label LIKE '%body mass%'
    """).fetchone()[0]
    print(f"  BMI risk edges resolved     : {bmi:,}")

    height = conn.execute("""
        SELECT COUNT(*) FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.relationship_type='INCREASES_RISK_OF'
          AND n.label LIKE '%height%'
    """).fetchone()[0]
    print(f"  Height risk edges resolved  : {height:,}")

    remaining = conn.execute(
        "SELECT COUNT(*) FROM edges WHERE target_system='GWAS_trait_label'"
    ).fetchone()[0]
    print(f"  Still unresolved (GWAS)     : {remaining:,}")

    s = gs.stats()
    print(f"\n── Graph stats ──────────────────────────────────────")
    print(f"  Total nodes : {s['total_nodes']:,}")
    print(f"  Total edges : {s['total_edges']:,}")
    print(f"\n✓  EFO loaded into {DB_PATH}")
