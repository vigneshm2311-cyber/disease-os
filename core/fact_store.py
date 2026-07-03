"""
fact_store.py — Raw fact ingestion layer.
One Fact = one measurement about one subject at one point in time.
"""
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

@dataclass
class Fact:
    subject_id:     str
    concept_system: str
    concept_code:   str
    concept_label:  str
    value_text:     Optional[str]   = None
    value_numeric:  Optional[float] = None
    unit:           Optional[str]   = None
    observed_at:    Optional[str]   = None
    source:         str             = "unknown"
    confidence:     float           = 1.0

SCHEMA = """
CREATE TABLE IF NOT EXISTS facts (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    subject_id      TEXT NOT NULL,
    concept_system  TEXT NOT NULL,
    concept_code    TEXT NOT NULL,
    concept_label   TEXT NOT NULL,
    value_text      TEXT,
    value_numeric   REAL,
    unit            TEXT,
    observed_at     TEXT NOT NULL,
    recorded_at     TEXT NOT NULL,
    source          TEXT NOT NULL,
    confidence      REAL NOT NULL DEFAULT 1.0
);"""

class FactStore:
    def __init__(self, db_path: str):
        self.db_path = db_path
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(db_path)
        self.conn.execute(SCHEMA)
        self.conn.commit()

    def add_fact(self, fact: Fact) -> int:
        observed_at = fact.observed_at or datetime.now(timezone.utc).isoformat()
        recorded_at = datetime.now(timezone.utc).isoformat()
        cur = self.conn.execute(
            "INSERT INTO facts (subject_id,concept_system,concept_code,"
            "concept_label,value_text,value_numeric,unit,observed_at,"
            "recorded_at,source,confidence) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            (fact.subject_id, fact.concept_system, fact.concept_code,
             fact.concept_label, fact.value_text, fact.value_numeric,
             fact.unit, observed_at, recorded_at, fact.source, fact.confidence)
        )
        self.conn.commit()
        return cur.lastrowid

    def facts_for_subject(self, subject_id: str) -> list:
        cur = self.conn.execute(
            "SELECT * FROM facts WHERE subject_id=? ORDER BY observed_at",
            (subject_id,)
        )
        cols = [c[0] for c in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]
