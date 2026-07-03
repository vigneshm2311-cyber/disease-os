"""
graph_store.py — SQLite-backed node and edge storage.
Deduplicates on insert. Merges xrefs across sources.
Indexes optimised for insurance code lookups.
"""
import json, sqlite3
from pathlib import Path
from core.node import Node
from core.edge import Edge

NODE_SCHEMA = """
CREATE TABLE IF NOT EXISTS nodes (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    primary_id      TEXT NOT NULL,
    primary_system  TEXT NOT NULL,
    label           TEXT NOT NULL,
    tier            INTEGER NOT NULL,
    entity_type     TEXT NOT NULL,
    xrefs           TEXT DEFAULT '{}',
    icd10_code      TEXT,
    icd11_code      TEXT,
    snomed_code     TEXT,
    loinc_code      TEXT,
    rxnorm_cui      TEXT,
    cpt_code        TEXT,
    hcc_code        TEXT,
    ndc_code        TEXT,
    synonyms        TEXT DEFAULT '[]',
    definition      TEXT,
    properties      TEXT DEFAULT '{}',
    source          TEXT NOT NULL,
    source_version  TEXT NOT NULL,
    confidence      REAL NOT NULL DEFAULT 1.0,
    loaded_at       TEXT NOT NULL,
    UNIQUE(primary_id, primary_system)
);"""

EDGE_SCHEMA = """
CREATE TABLE IF NOT EXISTS edges (
    id                       INTEGER PRIMARY KEY AUTOINCREMENT,
    source_id                TEXT NOT NULL,
    source_system            TEXT NOT NULL,
    target_id                TEXT NOT NULL,
    target_system            TEXT NOT NULL,
    relationship_type        TEXT NOT NULL,
    source_relationship_type TEXT DEFAULT '',
    effect_size              REAL,
    effect_unit              TEXT,
    direction                TEXT,
    confidence               REAL NOT NULL DEFAULT 0.5,
    feedback                 INTEGER NOT NULL DEFAULT 0,
    feedback_notes           TEXT DEFAULT '',
    primary_source           TEXT NOT NULL,
    imported_via             TEXT NOT NULL,
    study_design             TEXT DEFAULT 'unknown',
    population_context       TEXT DEFAULT 'general',
    tissue_context           TEXT DEFAULT 'systemic',
    species                  TEXT DEFAULT 'human',
    typical_latency          TEXT,
    source_version           TEXT NOT NULL,
    loaded_at                TEXT NOT NULL,
    UNIQUE(source_id, target_id, relationship_type, primary_source)
);"""

INDEXES = [
    "CREATE INDEX IF NOT EXISTS idx_nodes_entity_type ON nodes(entity_type);",
    "CREATE INDEX IF NOT EXISTS idx_nodes_tier        ON nodes(tier);",
    "CREATE INDEX IF NOT EXISTS idx_nodes_icd10       ON nodes(icd10_code);",
    "CREATE INDEX IF NOT EXISTS idx_nodes_rxnorm      ON nodes(rxnorm_cui);",
    "CREATE INDEX IF NOT EXISTS idx_nodes_loinc       ON nodes(loinc_code);",
    "CREATE INDEX IF NOT EXISTS idx_edges_source      ON edges(source_id, source_system);",
    "CREATE INDEX IF NOT EXISTS idx_edges_target      ON edges(target_id, target_system);",
    "CREATE INDEX IF NOT EXISTS idx_edges_rel_type    ON edges(relationship_type);",
    "CREATE INDEX IF NOT EXISTS idx_edges_confidence  ON edges(confidence);",
]


class GraphStore:

    def __init__(self, db_path):
        self.db_path = str(db_path)
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.execute("PRAGMA journal_mode=WAL;")
        self.conn.execute("PRAGMA synchronous=NORMAL;")
        self.conn.execute("PRAGMA cache_size=-64000;")   # 64MB cache
        self.conn.execute(NODE_SCHEMA)
        self.conn.execute(EDGE_SCHEMA)
        for idx in INDEXES:
            self.conn.execute(idx)
        self.conn.commit()

    # ── Internal ───────────────────────────────────────────────────────

    def _parse_node_row(self, row, cols) -> dict:
        d = dict(zip(cols, row))
        d["xrefs"]      = json.loads(d["xrefs"])
        d["synonyms"]   = json.loads(d["synonyms"])
        d["properties"] = json.loads(d["properties"])
        return d

    # ── Nodes ──────────────────────────────────────────────────────────

    def add_node(self, node: Node) -> str:
        d = node.to_dict()
        existing = self.conn.execute(
            "SELECT id, xrefs, synonyms FROM nodes "
            "WHERE primary_id=? AND primary_system=?",
            (node.primary_id, node.primary_system)
        ).fetchone()
        if existing:
            merged_xrefs    = {**json.loads(existing[1]), **node.xrefs}
            merged_synonyms = list(set(json.loads(existing[2]) + node.synonyms))
            self.conn.execute(
                "UPDATE nodes SET xrefs=?, synonyms=? WHERE id=?",
                (json.dumps(merged_xrefs), json.dumps(merged_synonyms), existing[0])
            )
        else:
            cols = ", ".join(d.keys())
            placeholders = ", ".join(["?"] * len(d))
            self.conn.execute(
                f"INSERT INTO nodes ({cols}) VALUES ({placeholders})",
                list(d.values())
            )
        self.conn.commit()
        return node.primary_id

    def add_nodes_batch(self, nodes: list, batch_size: int = 500):
        """Bulk insert for large loads — commits every batch_size nodes."""
        for i in range(0, len(nodes), batch_size):
            for node in nodes[i:i+batch_size]:
                self.add_node(node)

    def get_node(self, primary_id: str, primary_system: str) -> dict | None:
        cur = self.conn.execute(
            "SELECT * FROM nodes WHERE primary_id=? AND primary_system=?",
            (primary_id, primary_system)
        )
        row = cur.fetchone()
        if not row:
            return None
        return self._parse_node_row(row, [c[0] for c in cur.description])

    def find_by_icd10(self, icd10_code: str) -> list:
        cur = self.conn.execute(
            "SELECT * FROM nodes WHERE icd10_code=?", (icd10_code,)
        )
        cols = [c[0] for c in cur.description]
        return [self._parse_node_row(r, cols) for r in cur.fetchall()]

    def find_by_rxnorm(self, rxnorm_cui: str) -> list:
        cur = self.conn.execute(
            "SELECT * FROM nodes WHERE rxnorm_cui=?", (rxnorm_cui,)
        )
        cols = [c[0] for c in cur.description]
        return [self._parse_node_row(r, cols) for r in cur.fetchall()]

    def find_by_loinc(self, loinc_code: str) -> list:
        cur = self.conn.execute(
            "SELECT * FROM nodes WHERE loinc_code=?", (loinc_code,)
        )
        cols = [c[0] for c in cur.description]
        return [self._parse_node_row(r, cols) for r in cur.fetchall()]

    def find_by_umls_cui(self, cui: str) -> list:
        cur = self.conn.execute(
            "SELECT * FROM nodes WHERE json_extract(xrefs,'$.UMLS_CUI')=?",
            (cui,)
        )
        cols = [c[0] for c in cur.description]
        return [self._parse_node_row(r, cols) for r in cur.fetchall()]

    def search_label(self, term: str, limit: int = 20) -> list:
        cur = self.conn.execute(
            "SELECT * FROM nodes WHERE label LIKE ? LIMIT ?",
            (f"%{term}%", limit)
        )
        cols = [c[0] for c in cur.description]
        return [self._parse_node_row(r, cols) for r in cur.fetchall()]

    # ── Edges ──────────────────────────────────────────────────────────

    def add_edge(self, edge: Edge) -> bool:
        d = edge.to_dict()
        try:
            cols = ", ".join(d.keys())
            placeholders = ", ".join(["?"] * len(d))
            self.conn.execute(
                f"INSERT INTO edges ({cols}) VALUES ({placeholders})",
                list(d.values())
            )
            self.conn.commit()
            return True
        except sqlite3.IntegrityError:
            return False

    def edges_from(self, source_id: str,
                   relationship_type: str = None,
                   min_confidence: float = 0.0) -> list:
        q = "SELECT * FROM edges WHERE source_id=? AND confidence>=?"
        p = [source_id, min_confidence]
        if relationship_type:
            q += " AND relationship_type=?"
            p.append(relationship_type)
        cur = self.conn.execute(q, p)
        cols = [c[0] for c in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    def edges_to(self, target_id: str,
                 relationship_type: str = None,
                 min_confidence: float = 0.0) -> list:
        q = "SELECT * FROM edges WHERE target_id=? AND confidence>=?"
        p = [target_id, min_confidence]
        if relationship_type:
            q += " AND relationship_type=?"
            p.append(relationship_type)
        cur = self.conn.execute(q, p)
        cols = [c[0] for c in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    # ── Stats ──────────────────────────────────────────────────────────

    def stats(self) -> dict:
        return {
            "total_nodes":           self.conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0],
            "total_edges":           self.conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0],
            "nodes_by_tier":         dict(self.conn.execute("SELECT tier, COUNT(*) FROM nodes GROUP BY tier ORDER BY tier").fetchall()),
            "nodes_by_entity_type":  dict(self.conn.execute("SELECT entity_type, COUNT(*) FROM nodes GROUP BY entity_type ORDER BY 2 DESC").fetchall()),
            "edges_by_relationship": dict(self.conn.execute("SELECT relationship_type, COUNT(*) FROM edges GROUP BY relationship_type ORDER BY 2 DESC").fetchall()),
            "edges_by_study_design": dict(self.conn.execute("SELECT study_design, COUNT(*) FROM edges GROUP BY study_design ORDER BY 2 DESC").fetchall()),
        }
