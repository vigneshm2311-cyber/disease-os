"""
layers/tier1_molecular/load_umls.py

Bulk loads UMLS preprocessed TSVs into Disease-OS graph database.
Uses executemany() with large transaction batches for speed.

Run from project root:
    python3 layers/tier1_molecular/load_umls.py
"""

import sys, json, sqlite3
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node        import Node
from core.edge        import Edge
from core.graph_store import GraphStore
from core.config      import DB_PATH, PROCESSED_DIR
from core.sources.umls import UMLSSource

CONCEPTS_TSV  = PROCESSED_DIR / "umls_concepts.tsv"
RELATIONS_TSV = PROCESSED_DIR / "umls_relations.tsv"
BATCH_SIZE    = 5_000


def get_node_cols():
    sample = Node("x", "x", "x", 0, "Gene")
    return list(sample.to_dict().keys())


def get_edge_cols():
    sample = Edge(
        source_id="x", source_system="x",
        target_id="y", target_system="y",
        relationship_type="ASSOCIATED_WITH",
        primary_source="x", imported_via="x",
        source_version="x",
    )
    return list(sample.to_dict().keys())


def bulk_insert(conn, table, cols, batch):
    """Insert a batch inside one transaction. Skips duplicates silently."""
    if not batch:
        return
    col_str      = ", ".join(cols)
    placeholders = ", ".join(["?"] * len(cols))
    conn.execute("BEGIN")
    try:
        conn.executemany(
            f"INSERT OR IGNORE INTO {table} ({col_str}) VALUES ({placeholders})",
            batch
        )
        conn.execute("COMMIT")
    except Exception as e:
        conn.execute("ROLLBACK")
        raise e


def load_nodes(conn) -> dict:
    source  = UMLSSource()
    cols    = get_node_cols()
    batch   = []
    total   = 0

    print(f"[nodes] Streaming concepts TSV -> nodes table...")
    for node in source.nodes():
        batch.append(list(node.to_dict().values()))
        if len(batch) >= BATCH_SIZE:
            bulk_insert(conn, "nodes", cols, batch)
            total += len(batch)
            batch  = []
            if total % 100_000 == 0:
                print(f"  {total:,} nodes loaded...")

    if batch:
        bulk_insert(conn, "nodes", cols, batch)
        total += len(batch)

    print(f"  {total:,} nodes total")
    return {"nodes_loaded": total}


def load_edges(conn) -> dict:
    source  = UMLSSource()
    cols    = get_edge_cols()
    batch   = []
    total   = 0
    skipped = 0

    print(f"[edges] Streaming relations TSV -> edges table...")
    for edge in source.edges():
        if edge is None:
            skipped += 1
            continue
        batch.append(list(edge.to_dict().values()))
        if len(batch) >= BATCH_SIZE:
            bulk_insert(conn, "edges", cols, batch)
            total  += len(batch)
            batch   = []
            if total % 500_000 == 0:
                print(f"  {total:,} edges loaded, {skipped:,} skipped...")

    if batch:
        bulk_insert(conn, "edges", cols, batch)
        total += len(batch)

    print(f"  {total:,} edges total, {skipped:,} skipped")
    return {"edges_loaded": total, "edges_skipped": skipped}


def validate(gs):
    print("\n── Validation ───────────────────────────────────────")

    # Disease by ICD-10
    hits = gs.find_by_icd10("E11.9")
    if hits:
        n = hits[0]
        print(f"  ICD-10 E11.9    -> {n['label']}")
        print(f"  UMLS CUI        -> {n['xrefs'].get('UMLS_CUI','n/a')}")
    else:
        print("  ICD-10 E11.9    -> not found")

    # Drug by RxNorm
    hits = gs.find_by_rxnorm("860975")
    if hits:
        print(f"  RxNorm 860975   -> {hits[0]['label']}")
    else:
        print("  RxNorm 860975   -> not found")

    # Lab by LOINC
    hits = gs.find_by_loinc("4548-4")
    if hits:
        print(f"  LOINC 4548-4    -> {hits[0]['label']}")
    else:
        print("  LOINC 4548-4    -> not found")

    # Label search
    results = gs.search_label("diabetes", limit=5)
    print(f"\n  Label search 'diabetes' -> {len(results)} hits (top 3):")
    for r in results[:3]:
        print(f"    [{r['primary_system']}] {r['primary_id']} — {r['label']}")

    # Stats
    print(f"\n── Graph stats ──────────────────────────────────────")
    s = gs.stats()
    print(f"  Total nodes      : {s['total_nodes']:,}")
    print(f"  Total edges      : {s['total_edges']:,}")
    print(f"  Nodes by tier    : {s['nodes_by_tier']}")
    print(f"  Nodes by type    : {s['nodes_by_entity_type']}")
    print(f"  Edges by rel     : {s['edges_by_relationship']}")


if __name__ == "__main__":
    if not CONCEPTS_TSV.exists() or not RELATIONS_TSV.exists():
        print("ERROR: Run preprocess() first.")
        print("  python3 -c \"from core.sources.umls import UMLSSource; UMLSSource().preprocess()\"")
        sys.exit(1)

    print(f"[Disease-OS] Loading UMLS 2025AB")
    print(f"  DB           : {DB_PATH}")
    print(f"  Concepts TSV : {CONCEPTS_TSV.stat().st_size // 1_000_000}MB")
    print(f"  Relations TSV: {RELATIONS_TSV.stat().st_size // 1_000_000}MB")

    gs   = GraphStore(DB_PATH)
    conn = gs.conn

    # Performance pragmas for bulk load
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")
    conn.execute("PRAGMA temp_store=MEMORY;")

    t0 = datetime.now()

    print(f"\n[1/2] Loading nodes...")
    node_stats = load_nodes(conn)
    t1 = datetime.now()
    print(f"  Done in {int((t1-t0).total_seconds())}s")

    print(f"\n[2/2] Loading edges...")
    edge_stats = load_edges(conn)
    t2 = datetime.now()
    print(f"  Done in {int((t2-t1).total_seconds())}s")

    # Restore safe mode
    conn.execute("PRAGMA synchronous=NORMAL;")

    print(f"\n[Disease-OS] Total time: {int((t2-t0).total_seconds())}s")
    validate(gs)
    print(f"\n✓  UMLS 2025AB loaded into {DB_PATH}")
