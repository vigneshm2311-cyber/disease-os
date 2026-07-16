"""
core/sources/orphanet.py

Orphanet rare disease loader for Disease-OS.
Input: data/raw/orphanet.xml (52MB, CC-BY-4.0)

Adds 11,645 rare disease nodes with Orphanet IDs as primary,
cross-referenced to OMIM, ICD-10, MONDO.

Primary ID: Orphanet:166024 (OrphaCode with prefix)
Primary system: Orphanet
Tier: 7 (Disease label)
"""

import sys, json
from pathlib import Path
from typing import Generator
from xml.etree import ElementTree as ET

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource

RAW = Path.home() / "disease-os/data/raw/orphanet.xml"


def _parse_orphanet(path: Path) -> list[dict]:
    print(f"[Orphanet] Parsing {path.name}...")
    tree  = ET.parse(str(path))
    root  = tree.getroot()
    terms = []

    for disorder in root.iter("Disorder"):
        orpha_code = disorder.findtext("OrphaCode")
        if not orpha_code:
            continue

        name_el = disorder.find("Name[@lang='en']")
        name    = name_el.text.strip() if name_el is not None and name_el.text else ""
        if not name:
            continue

        # Disorder type
        dtype_el = disorder.find("DisorderType/Name[@lang='en']")
        dtype    = dtype_el.text.strip() if dtype_el is not None else "Disease"

        # Synonyms
        synonyms = []
        for syn in disorder.findall("SynonymList/Synonym[@lang='en']"):
            if syn.text:
                synonyms.append(syn.text.strip())

        # External cross-references
        xrefs    = {}
        icd10    = None
        for xref in disorder.findall(
                "ExternalReferenceList/ExternalReference"):
            source_el = xref.find("Source")
            ref_el    = xref.find("Reference")
            if source_el is None or ref_el is None:
                continue
            source = source_el.text.strip() if source_el.text else ""
            ref    = ref_el.text.strip()    if ref_el.text    else ""
            if not source or not ref:
                continue
            if source == "OMIM":
                xrefs["OMIM"] = ref
            elif source in ("ICD-10","ICD10"):
                xrefs["ICD-10"] = ref
                icd10 = ref
            elif source == "MeSH":
                xrefs["MeSH"] = ref
            elif source == "UMLS":
                xrefs["UMLS_CUI"] = ref
            elif source == "MONDO":
                xrefs["MONDO"] = ref
            elif source == "MedGen":
                xrefs["MedGen"] = ref

        terms.append({
            "orpha_code": orpha_code,
            "name":       name,
            "dtype":      dtype,
            "synonyms":   synonyms[:15],
            "xrefs":      xrefs,
            "icd10":      icd10,
        })

    print(f"[Orphanet] {len(terms):,} disorders parsed")
    return terms


class OrphanetSource(BaseSource):
    source_name    = "Orphanet"
    source_version = "2026-06"

    def __init__(self):
        self._terms = None

    def _load(self):
        if self._terms is None:
            self._terms = _parse_orphanet(RAW)

    def nodes(self) -> Generator[Node, None, None]:
        self._load()
        for t in self._terms:
            yield Node(
                primary_id     = f"Orphanet:{t['orpha_code']}",
                primary_system = "Orphanet",
                label          = t["name"],
                tier           = 7,
                entity_type    = "Disease_scientific",
                xrefs          = t["xrefs"],
                synonyms       = t["synonyms"],
                icd10_code     = t["icd10"],
                source         = self.source_name,
                source_version = self.source_version,
                confidence     = 0.90,
            )

    def edges(self) -> Generator[Edge, None, None]:
        # ISA edges handled separately via disorder hierarchy
        return
        yield

    def normalize_confidence(self, raw_value=None) -> float:
        return 0.90
