"""
core/sources/efo.py

EFO (Experimental Factor Ontology) v3.91 adapter for Disease-OS.

Input : data/raw/efo/efo.obo (84MB OBO format)

What EFO adds:
  1. Resolves 598K unresolved GWAS trait label edges
     by providing a name/synonym index for GWAS-relevant terms
  2. Adds EFO disease/phenotype nodes with cross-references
     to MONDO, HP, Orphanet, ICD10, OMIM

Strategy:
  - Parse OBO file, extract all [Term] blocks
  - Keep only terms with efo:gwas_trait="true" OR
    from disease/phenotype namespaces (EFO, MONDO, HP, Orphanet)
  - Build name+synonym index for GWAS trait label matching
  - Write nodes to graph
  - Run GWAS edge resolution using the index
"""

import re
import sys
import json
import sqlite3
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS, PROCESSED_DIR

RAW_DIR  = Path.home() / "disease-os" / "data" / "raw" / "efo"
EFO_OBO  = RAW_DIR / "efo.obo"
EFO_TSV  = PROCESSED_DIR / "efo_terms.tsv"

# Prefixes we load as nodes
TARGET_PREFIXES = {"efo", "MONDO", "HP", "Orphanet"}

# Entity type mapping by prefix
PREFIX_ENTITY = {
    "efo":      ("Disease_scientific", 7),
    "MONDO":    ("Disease_scientific", 7),
    "HP":       ("Phenotype",          6),
    "Orphanet": ("Disease_scientific", 7),
}


def _parse_obo(path: Path) -> list[dict]:
    """
    Parse OBO file into list of term dicts.
    Returns only terms from TARGET_PREFIXES.
    """
    terms = []
    current = {}
    in_term = False

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            if line == "[Term]":
                if current and _keep_term(current):
                    terms.append(current)
                current  = {"synonyms": [], "xrefs": [], "is_gwas": False}
                in_term  = True
                continue

            if line.startswith("[") and line != "[Term]":
                in_term = False
                continue

            if not in_term or not line:
                continue

            if line.startswith("id: "):
                current["id"] = line[4:].strip()

            elif line.startswith("name: "):
                current["name"] = line[6:].strip()

            elif line.startswith("def: "):
                # Extract text between first pair of quotes
                m = re.match(r'def: "([^"]*)"', line)
                if m:
                    current["def"] = m.group(1)[:500]

            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "([^"]*)"', line)
                if m:
                    current["synonyms"].append(m.group(1).strip())

            elif line.startswith("xref: "):
                current["xrefs"].append(line[6:].strip())

            elif 'gwas_trait "true"' in line:
                current["is_gwas"] = True

            elif line.startswith("is_obsolete: true"):
                current["obsolete"] = True

    # Don't forget last term
    if current and _keep_term(current):
        terms.append(current)

    return terms


def _keep_term(term: dict) -> bool:
    """Keep term if it's from a target prefix and not obsolete."""
    if term.get("obsolete"):
        return False
    tid = term.get("id", "")
    prefix = tid.split(":")[0]
    return prefix in TARGET_PREFIXES


def _extract_xref_codes(xrefs: list) -> dict:
    """
    Extract cross-reference codes from xref list.
    e.g. "ICD10WHO:G30" -> {"ICD10": "G30"}
    """
    result = {}
    for xref in xrefs:
        xref = xref.strip()
        if ":" not in xref:
            continue
        system, _, code = xref.partition(":")
        system = system.strip()
        code   = code.strip().split(" ")[0]  # remove trailing comments
        if system in ("ICD10WHO", "ICD10CM", "ICD10"):
            result["ICD-10"] = code
        elif system == "OMIM":
            result["OMIM"] = code
        elif system == "NCIT":
            result["NCI"] = code
        elif system == "DOID":
            result["DOID"] = code
        elif system in ("MESH", "MeSH"):
            result["MeSH"] = code
        elif system == "SNOMEDCT":
            result["SNOMED"] = code
        elif system.startswith("MONDO"):
            result["MONDO"] = code
    return result


class EFOSource(BaseSource):
    """
    Loads EFO v3.91 into Disease-OS and resolves GWAS trait labels.

    Usage:
        source = EFOSource()
        source.preprocess()
        source.load_into(graph_store)
        source.resolve_gwas_edges(conn)
    """

    source_name    = "EFO"
    source_version = "3.91.0"

    def __init__(self):
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        self._terms: list[dict] = []

    def preprocess(self, force: bool = False):
        """Parse OBO and write compact TSV."""
        if not force and EFO_TSV.exists():
            print("[EFO] Processed file exists — skipping.")
            return

        print(f"[EFO] Parsing {EFO_OBO.name} ({EFO_OBO.stat().st_size//1_000_000}MB)...")
        terms = _parse_obo(EFO_OBO)

        gwas_terms = sum(1 for t in terms if t.get("is_gwas"))
        print(f"[EFO] {len(terms):,} terms parsed | {gwas_terms:,} GWAS-tagged")

        import csv
        with open(EFO_TSV, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "id", "name", "definition", "is_gwas",
                "synonyms_json", "xref_codes_json",
            ])
            for t in terms:
                writer.writerow([
                    t.get("id", ""),
                    t.get("name", ""),
                    t.get("def", ""),
                    "1" if t.get("is_gwas") else "0",
                    json.dumps(t.get("synonyms", [])),
                    json.dumps(_extract_xref_codes(t.get("xrefs", []))),
                ])

        print(f"[EFO] Written to {EFO_TSV.name}")

    def _load_terms(self):
        """Load processed TSV into memory."""
        if self._terms:
            return
        import csv
        with open(EFO_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                self._terms.append({
                    "id":       row["id"],
                    "name":     row["name"],
                    "def":      row["definition"],
                    "is_gwas":  row["is_gwas"] == "1",
                    "synonyms": json.loads(row["synonyms_json"] or "[]"),
                    "xrefs":    json.loads(row["xref_codes_json"] or "{}"),
                })

    def nodes(self) -> Generator[Node, None, None]:
        """Yield EFO/MONDO/HP disease nodes."""
        if not EFO_TSV.exists():
            raise FileNotFoundError("Run preprocess() first.")
        self._load_terms()

        for t in self._terms:
            tid    = t["id"]
            prefix = tid.split(":")[0]
            entity_type, tier = PREFIX_ENTITY.get(prefix, ("Disease_scientific", 7))

            # Build primary ID — use EFO ID as primary
            primary_id     = tid.replace("efo:", "EFO:")
            primary_system = "EFO"

            xrefs = dict(t["xrefs"])
            # Add ICD-10 as explicit field if present
            icd10 = xrefs.get("ICD-10")

            yield Node(
                primary_id     = primary_id,
                primary_system = primary_system,
                label          = t["name"],
                tier           = tier,
                entity_type    = entity_type,
                xrefs          = xrefs,
                synonyms       = t["synonyms"][:20],
                definition     = t["def"] or None,
                icd10_code     = icd10,
                properties     = {"is_gwas_trait": t["is_gwas"]},
                source         = self.source_name,
                source_version = self.source_version,
                confidence     = 0.90,
            )

    def edges(self) -> Generator[Edge, None, None]:
        """No structural edges from EFO itself — edges come from resolution."""
        return
        yield

    def normalize_confidence(self, raw_value=None) -> float:
        return 0.90

    def build_gwas_index(self) -> dict:
        """
        Build {label_lower: (efo_id, efo_system, canonical_name)} index
        for matching GWAS trait labels.
        Includes both names and synonyms.
        """
        self._load_terms()
        index = {}

        for t in self._terms:
            if not t["name"]:
                continue

            efo_id = t["id"].replace("efo:", "EFO:")
            val    = (efo_id, "EFO", t["name"])

            # Index by name
            index[t["name"].lower().strip()] = val

            # Index by each synonym
            for syn in t["synonyms"]:
                if syn and len(syn) > 3:
                    index[syn.lower().strip()] = val

        print(f"[EFO] GWAS index built: {len(index):,} name/synonym entries")
        return index

    def resolve_gwas_edges(self, conn: sqlite3.Connection) -> dict:
        """
        Resolve remaining GWAS_trait_label edges using EFO name/synonym index.
        Uses delete-then-reinsert to avoid UNIQUE constraint violations.
        """
        index = self.build_gwas_index()

        # Get distinct unresolved trait labels
        traits = conn.execute("""
            SELECT target_id, COUNT(*) as n
            FROM edges
            WHERE target_system = 'GWAS_trait_label'
            GROUP BY target_id
            ORDER BY n DESC
        """).fetchall()

        print(f"[EFO] Resolving {len(traits):,} unresolved GWAS trait labels...")

        resolved_map = {}   # trait_label -> (efo_id, system, name)
        unresolved   = []

        for trait_label, n_edges in traits:
            clean = re.sub(
                r'\s*\(SNP x SNP.*?\)|\s*\(.*?df\)|\s*\(age.*?\)|'
                r'\s*\(sex.*?\)|\s*\(joint.*?\)|\s*\(adjusted.*?\)',
                '', trait_label, flags=re.IGNORECASE
            ).strip().rstrip(",").strip().lower()

            if not clean or len(clean) < 3:
                unresolved.append((trait_label, n_edges))
                continue

            # Exact match
            match = index.get(clean)

            # Partial match — try removing parenthetical qualifiers
            if not match:
                base = re.sub(r'\s*\([^)]*\)', '', clean).strip()
                if base:
                    match = index.get(base)

            if match:
                resolved_map[trait_label] = match
            else:
                unresolved.append((trait_label, n_edges))

        print(f"[EFO] Matched {len(resolved_map):,} / {len(traits):,} traits")
        print(f"[EFO] Unresolved: {len(unresolved):,}")

        # Get edge column names
        cols = [r[1] for r in conn.execute(
            "PRAGMA table_info(edges)"
        ).fetchall() if r[1] != "id"]
        col_str = ", ".join(cols)
        ph      = ", ".join(["?"] * len(cols))

        edges_updated = 0

        for trait_label, (efo_id, sys, name) in resolved_map.items():
            # Fetch originals
            rows = conn.execute(
                f"SELECT {col_str} FROM edges "
                f"WHERE target_id=? AND target_system='GWAS_trait_label'",
                (trait_label,)
            ).fetchall()

            if not rows:
                continue

            col_list = cols
            tid_idx  = col_list.index("target_id")
            tsys_idx = col_list.index("target_system")

            resolved = []
            for r in rows:
                r = list(r)
                r[tid_idx]  = efo_id
                r[tsys_idx] = sys
                resolved.append(tuple(r))

            conn.execute("BEGIN")
            conn.execute(
                "DELETE FROM edges WHERE target_id=? AND target_system='GWAS_trait_label'",
                (trait_label,)
            )
            conn.executemany(
                f"INSERT OR IGNORE INTO edges ({col_str}) VALUES ({ph})",
                resolved
            )
            conn.execute("COMMIT")
            edges_updated += len(rows)

        remaining = conn.execute(
            "SELECT COUNT(*) FROM edges WHERE target_system='GWAS_trait_label'"
        ).fetchone()[0]

        return {
            "traits_resolved": len(resolved_map),
            "traits_total":    len(traits),
            "edges_updated":   edges_updated,
            "remaining":       remaining,
        }
