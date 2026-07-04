"""
core/id_registry.py

Human-readable ID registry for Disease-OS.

Answers three questions instantly without SQL:
  1. "What did this CUI resolve to?"
  2. "What is the canonical primary ID for this concept?"
  3. "What other IDs does this node have?"

The registry is built from the cui_map table and exported
to a TSV file at data/processed/id_registry.tsv for inspection.

Also provides a Python lookup class used by all source adapters
so they can resolve IDs at load time rather than post-hoc.

Run directly to rebuild the registry:
    python3 core/id_registry.py
"""

import sys
import csv
import sqlite3
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from core.config import DB_PATH, PROCESSED_DIR

REGISTRY_TSV = PROCESSED_DIR / "id_registry.tsv"


def export_registry(conn) -> int:
    """
    Export cui_map + node xrefs into a human-readable TSV.
    Columns:
      old_id        — the CUI or alternate ID being translated
      old_system    — which vocabulary the old_id belongs to
      primary_id    — canonical primary ID now used in the graph
      primary_system— vocabulary of the primary ID
      label         — human-readable name
      tier          — which Disease-OS tier this node belongs to
      all_xrefs     — pipe-separated list of all other known IDs
    """
    print("[registry] Building ID registry from database...")

    rows = conn.execute("""
        SELECT
            c.cui            as old_id,
            'UMLS_CUI'       as old_system,
            c.primary_id,
            c.primary_system,
            c.label,
            c.tier,
            n.xrefs
        FROM cui_map c
        LEFT JOIN nodes n ON n.primary_id     = c.primary_id
                          AND n.primary_system = c.primary_system
        ORDER BY c.tier, c.primary_system, c.primary_id
    """).fetchall()

    with open(REGISTRY_TSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "old_id", "old_system",
            "primary_id", "primary_system",
            "label", "tier", "all_xrefs"
        ])
        for row in rows:
            old_id, old_system, pid, psys, label, tier, xrefs_json = row
            # Flatten xrefs into pipe-separated key=value pairs
            import json
            xrefs = json.loads(xrefs_json) if xrefs_json else {}
            xrefs_str = " | ".join(
                f"{k}:{v}" for k, v in xrefs.items()
                if k != "UMLS_CUI"   # CUI already in old_id column
            )
            writer.writerow([
                old_id, old_system,
                pid, psys,
                (label or "")[:100],
                tier,
                xrefs_str,
            ])

    print(f"[registry] {len(rows):,} entries written to {REGISTRY_TSV}")
    return len(rows)


class IDRegistry:
    """
    In-memory lookup table for ID resolution.
    Used by source adapters to resolve IDs at load time.

    Usage:
        registry = IDRegistry()
        primary = registry.resolve("C0011860", "UMLS_CUI")
        # returns ("E11.9", "ICD-10-CM") or None
    """

    def __init__(self):
        self._cui_to_primary: dict[str, tuple] = {}
        self._primary_to_xrefs: dict[tuple, dict] = {}
        self._loaded = False

    def load(self, conn):
        """Load from database into memory."""
        import json
        rows = conn.execute("""
            SELECT c.cui, c.primary_id, c.primary_system,
                   n.xrefs
            FROM cui_map c
            LEFT JOIN nodes n ON n.primary_id     = c.primary_id
                              AND n.primary_system = c.primary_system
        """).fetchall()

        for cui, pid, psys, xrefs_json in rows:
            self._cui_to_primary[cui] = (pid, psys)
            xrefs = json.loads(xrefs_json) if xrefs_json else {}
            self._primary_to_xrefs[(pid, psys)] = xrefs

        self._loaded = True
        return len(rows)

    def resolve_cui(self, cui: str) -> tuple | None:
        """
        Resolve a UMLS CUI to (primary_id, primary_system).
        Returns None if CUI not in graph.
        """
        return self._cui_to_primary.get(cui)

    def get_xrefs(self, primary_id: str, primary_system: str) -> dict:
        """Return all cross-references for a node."""
        return self._primary_to_xrefs.get((primary_id, primary_system), {})

    def resolve_any(self, id_value: str, system: str) -> tuple | None:
        """
        Resolve any ID to primary.
        Currently handles UMLS_CUI; extensible to other systems.
        """
        if system == "UMLS_CUI":
            return self.resolve_cui(id_value)
        # Future: add SNOMED -> primary, OMIM -> primary etc.
        return None

    @property
    def size(self) -> int:
        return len(self._cui_to_primary)


def sample_lookups(conn):
    """Show sample lookups to confirm registry works."""
    print("\n── Sample ID resolutions ────────────────────────────")

    samples = [
        ("C0011860", "Type 2 Diabetes"),
        ("C0027672", "BRCA1-related Breast Cancer"),
        ("C0006142", "Breast Cancer"),
        ("C0002395", "Alzheimer's Disease"),
        ("C0020538", "Hypertension"),
        ("C0021640", "Insulin"),
        ("C0025598", "Metformin"),
        ("C0004057", "Aspirin"),
    ]

    rows = []
    for cui, description in samples:
        row = conn.execute("""
            SELECT c.cui, c.primary_id, c.primary_system, c.label, c.tier
            FROM cui_map c
            WHERE c.cui = ?
        """, (cui,)).fetchone()

        if row:
            rows.append(row + (description,))

    print(f"  {'Description':<30} {'Old CUI':<12} {'New Primary ID':<15} "
          f"{'System':<15} {'Tier'}")
    print(f"  {'-'*30} {'-'*12} {'-'*15} {'-'*15} {'-'*4}")
    for row in rows:
        cui, pid, psys, label, tier, desc = row
        print(f"  {desc:<30} {cui:<12} {pid:<15} {psys:<15} {tier}")

    # Tier breakdown of resolved IDs
    print(f"\n── Resolutions by tier ──────────────────────────────")
    tier_counts = conn.execute("""
        SELECT tier, primary_system, COUNT(*) as n
        FROM cui_map
        GROUP BY tier, primary_system
        ORDER BY tier, n DESC
    """).fetchall()

    current_tier = None
    for tier, system, n in tier_counts:
        if tier != current_tier:
            tier_names = {
                1:"Molecular", 2:"Networks", 3:"Cellular",
                4:"Tissue/Organ", 5:"Systemic", 6:"Phenotype",
                7:"Disease", 8:"Behavior", 9:"Social",
                10:"Healthcare", 11:"Population"
            }
            print(f"\n  Tier {tier} — {tier_names.get(tier,'?')}:")
            current_tier = tier
        print(f"    {system:<20}: {n:,}")


if __name__ == "__main__":
    print(f"[Disease-OS] ID Registry builder")
    print(f"  DB      : {DB_PATH}")
    print(f"  Output  : {REGISTRY_TSV}")

    conn = sqlite3.connect(str(DB_PATH))

    n = export_registry(conn)
    sample_lookups(conn)

    print(f"\n── How to use the registry ──────────────────────────")
    print(f"  TSV file (open in Excel/Numbers/any editor):")
    print(f"    {REGISTRY_TSV}")
    print(f"\n  Columns:")
    print(f"    old_id        — original UMLS CUI")
    print(f"    old_system    — always UMLS_CUI")
    print(f"    primary_id    — what the graph uses now")
    print(f"    primary_system— which vocabulary")
    print(f"    label         — human-readable name")
    print(f"    tier          — Disease-OS tier (1-11)")
    print(f"    all_xrefs     — all other IDs for same concept")
    print(f"\n  In Python:")
    print(f"    from core.id_registry import IDRegistry")
    print(f"    r = IDRegistry()")
    print(f"    r.load(conn)")
    print(f"    r.resolve_cui('C0011860')  # -> ('E11.9', 'ICD-10-CM')")
    print(f"\n✓  Registry built: {n:,} entries")
    conn.close()
