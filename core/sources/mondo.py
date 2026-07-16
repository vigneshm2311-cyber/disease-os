"""
core/sources/mondo.py

MONDO (Monarch Disease Ontology) loader for Disease-OS.
Input: data/raw/mondo.obo (51MB OBO format)

Adds 20K+ disease nodes with MONDO IDs as primary,
cross-referenced to OMIM, DOID, ICD-10, Orphanet, MeSH.

Primary ID: MONDO:0000004
Primary system: MONDO_disease
Tier: 7 (Disease label)

Filters out obsolete terms.
"""

import sys, re, json
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource

RAW = Path.home() / "disease-os/data/raw/mondo.obo"


def _parse_mondo(path: Path) -> list[dict]:
    print(f"[MONDO] Parsing {path.name}...")
    terms   = []
    current = {}
    in_term = False

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()

            if line == "[Term]":
                if current and not current.get("obsolete") \
                        and current.get("id","").startswith("MONDO:"):
                    terms.append(current)
                current = {"synonyms": [], "xrefs": {}, "is_a": []}
                in_term = True
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
                m = re.match(r'def: "([^"]*)"', line)
                if m:
                    current["def"] = m.group(1)[:400]

            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "([^"]*)"', line)
                if m:
                    current["synonyms"].append(m.group(1).strip())

            elif line.startswith("xref: "):
                xref = line[6:].split("{")[0].strip()
                if ":" in xref:
                    sys_raw, _, code = xref.partition(":")
                    sys_raw = sys_raw.strip()
                    code    = code.strip()
                    if sys_raw == "OMIM" or sys_raw == "OMIMPS":
                        current["xrefs"]["OMIM"] = code
                    elif sys_raw in ("ICD10", "ICD-10", "ICD10CM"):
                        current["xrefs"]["ICD-10"] = code
                    elif sys_raw == "Orphanet":
                        current["xrefs"]["Orphanet"] = code
                    elif sys_raw in ("MeSH", "MESH"):
                        current["xrefs"]["MeSH"] = code
                    elif sys_raw == "DOID":
                        current["xrefs"]["DOID"] = code
                    elif sys_raw == "UMLS":
                        current["xrefs"]["UMLS_CUI"] = code
                    elif sys_raw == "MedGen":
                        current["xrefs"]["MedGen"] = code
                    elif sys_raw == "EFO":
                        current["xrefs"]["EFO"] = code

            elif line.startswith("is_a: "):
                parent = line[6:].split("!")[0].strip()
                if parent.startswith("MONDO:"):
                    current["is_a"].append(parent)

            elif line.startswith("is_obsolete: true"):
                current["obsolete"] = True

    # Last term
    if current and not current.get("obsolete") \
            and current.get("id","").startswith("MONDO:"):
        terms.append(current)

    print(f"[MONDO] {len(terms):,} non-obsolete terms parsed")
    return terms


class MondoSource(BaseSource):
    source_name    = "MONDO"
    source_version = "2026-06"

    def __init__(self):
        self._terms = None

    def _load(self):
        if self._terms is None:
            self._terms = _parse_mondo(RAW)

    def nodes(self) -> Generator[Node, None, None]:
        self._load()
        for t in self._terms:
            if not t.get("name"):
                continue
            xrefs  = dict(t["xrefs"])
            icd10  = xrefs.get("ICD-10")
            orpha  = xrefs.get("Orphanet")
            # Add Orphanet as full ID if present
            if orpha:
                xrefs["Orphanet_full"] = f"Orphanet:{orpha}"

            yield Node(
                primary_id     = t["id"],
                primary_system = "MONDO_disease",
                label          = t["name"],
                tier           = 7,
                entity_type    = "Disease_scientific",
                xrefs          = xrefs,
                synonyms       = t["synonyms"][:15],
                definition     = t.get("def"),
                icd10_code     = icd10,
                source         = self.source_name,
                source_version = self.source_version,
                confidence     = 0.90,
            )

    def edges(self) -> Generator[Edge, None, None]:
        """ISA edges within MONDO hierarchy."""
        self._load()
        for t in self._terms:
            if not t.get("name"):
                continue
            for parent in t.get("is_a", []):
                try:
                    yield Edge(
                        source_id                = t["id"],
                        source_system            = "MONDO_disease",
                        target_id                = parent,
                        target_system            = "MONDO_disease",
                        relationship_type        = "ISA",
                        source_relationship_type = "mondo_is_a",
                        confidence               = 1.0,
                        primary_source           = f"MONDO_{self.source_version}",
                        imported_via             = f"MONDO_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                except ValueError:
                    continue

    def normalize_confidence(self, raw_value=None) -> float:
        return 0.90
