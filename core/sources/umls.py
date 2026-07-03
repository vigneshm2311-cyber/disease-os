"""
core/sources/umls.py

UMLS 2025AB adapter for Disease-OS.

File sizes:
  MRCONSO.RRF  2.1 GB  — concepts and their codes across all SABs
  MRREL.RRF    5.7 GB  — relationships between concepts
  MRSTY.RRF    201 MB  — semantic types per concept
  MRDEF.RRF    131 MB  — definitions

Strategy: pre-process raw RRFs once into compact TSV files in
data/processed/. Every subsequent run reads from those small files.

Pipeline:
  1. _build_semtype_index()  MRSTY   -> {CUI: {tui, sty}}  in memory
  2. _extract_concepts()     MRCONSO -> processed/umls_concepts.tsv
  3. _extract_relations()    MRREL   -> processed/umls_relations.tsv
  4. nodes()                 reads concepts TSV -> yields Node objects
  5. edges()                 reads relations TSV -> yields Edge objects
"""

import csv
import sys
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import (
    UMLS_MRCONSO, UMLS_MRREL, UMLS_MRSTY,
    PROCESSED_DIR, SOURCE_VERSIONS,
    UMLS_TARGET_SABS, UMLS_TARGET_SEMTYPES,
    CONFIDENCE_RULES,
)

CONCEPTS_TSV  = PROCESSED_DIR / "umls_concepts.tsv"
RELATIONS_TSV = PROCESSED_DIR / "umls_relations.tsv"

# MRCONSO column indices
C_CUI      = 0
C_LAT      = 1
C_ISPREF   = 6
C_SAB      = 11
C_CODE     = 13
C_STR      = 14
C_SUPPRESS = 16

# MRREL column indices
R_CUI1     = 0
R_REL      = 3
R_CUI2     = 4
R_RELA     = 7
R_SAB      = 10
R_SUPPRESS = 14

# MRSTY column indices
T_CUI = 0
T_TUI = 1
T_STY = 3

# Semantic type -> (entity_type, tier, preferred_primary_system)
SEMTYPE_TO_ENTITY = {
    "T047": ("Disease_clinical",  7, "ICD-10-CM"),
    "T048": ("Disease_clinical",  7, "ICD-10-CM"),
    "T191": ("Disease_clinical",  7, "ICD-10-CM"),
    "T046": ("Disease_clinical",  7, "ICD-10-CM"),
    "T184": ("ClinicalFinding",   6, "SNOMED"),
    "T033": ("ClinicalFinding",   6, "SNOMED"),
    "T201": ("ClinicalFinding",   6, "SNOMED"),
    "T059": ("LabTest",           6, "LOINC"),
    "T034": ("ClinicalFinding",   6, "LOINC"),
    "T121": ("Drug_clinical",    10, "RxNorm"),
    "T200": ("Drug_clinical",    10, "RxNorm"),
    "T116": ("Protein",           1, "UniProt"),
    "T028": ("Gene",              1, "NCBI_Gene"),
    "T086": ("Gene",              1, "NCBI_Gene"),
    "T023": ("Anatomy",           4, "UBERON"),
    "T025": ("CellType",          3, "CellOntology"),
    "T043": ("CellType",          3, "CellOntology"),
    "T044": ("Protein",           1, "UniProt"),
    "T045": ("Gene",              1, "NCBI_Gene"),
    "T038": ("ClinicalFinding",   6, "SNOMED"),
    "T031": ("Metabolite",        1, "HMDB"),
    "T109": ("Metabolite",        1, "HMDB"),
    "T123": ("Protein",           1, "UniProt"),
    "T058": ("Procedure",        10, "CPT"),
    "T074": ("Procedure",        10, "CPT"),
}

SAB_TO_SYSTEM = {
    "ICD10CM":     "ICD-10-CM",
    "SNOMEDCT_US": "SNOMED",
    "LNC":         "LOINC",
    "RXNORM":      "RxNorm",
    "HPO":         "HPO",
    "MONDO":       "MONDO",
    "OMIM":        "OMIM",
    "NCI":         "NCI",
    "MSH":         "MeSH",
    "CPT":         "CPT",
    "GO":          "GO",
    "FMA":         "FMA",
}

UMLS_REL_MAP = {
    "CHD":                    "ISA",
    "PAR":                    "ISA",
    "RB":                     "ISA",
    "RN":                     "ISA",
    "isa":                    "ISA",
    "part_of":                "PART_OF",
    "has_part":               "PART_OF",
    "component_of":           "PART_OF",
    "causative_agent_of":     "CAUSES",
    "has_causative_agent":    "CAUSES",
    "may_treat":              "TREATS",
    "may_prevent":            "PROTECTS_AGAINST",
    "may_be_treated_by":      "TREATS",
    "has_finding_site":       "LOCATED_IN",
    "finding_site_of":        "LOCATED_IN",
    "has_manifestation":      "HAS_SYMPTOM",
    "manifestation_of":       "HAS_SYMPTOM",
    "has_sign_or_symptom":    "HAS_SYMPTOM",
    "sign_or_symptom_of":     "HAS_SYMPTOM",
    "interprets":             "INDICATES",
    "is_interpreted_by":      "INDICATES",
    "diagnoses":              "DIAGNOSED_BY",
    "diagnosed_by":           "DIAGNOSED_BY",
    "gene_product_of":        "ENCODES",
    "has_gene_product":       "ENCODES",
    "encodes":                "ENCODES",
    "encoded_by":             "ENCODES",
    "positively_regulates":   "UPREGULATES",
    "negatively_regulates":   "DOWNREGULATES",
    "inhibits":               "INHIBITS",
    "activated_by":           "ACTIVATES",
    "inhibited_by":           "INHIBITS",
    "has_location":           "LOCATED_IN",
    "location_of":            "LOCATED_IN",
    "associated_with":        "ASSOCIATED_WITH",
    "AQ":                     "ASSOCIATED_WITH",
    "QB":                     "ASSOCIATED_WITH",
    "RO":                     "ASSOCIATED_WITH",
    "SY":                     "ASSOCIATED_WITH",
}


class UMLSSource(BaseSource):
    """
    Streams UMLS 2025AB into Disease-OS Node and Edge objects.

    Usage:
        source = UMLSSource()
        source.preprocess()           # run once — processes raw RRFs
        source.load_into(graph_store) # load nodes + edges into graph
    """

    source_name    = "UMLS"
    source_version = SOURCE_VERSIONS["UMLS"]

    def __init__(self):
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # ── Public API ─────────────────────────────────────────────────────

    def preprocess(self, force: bool = False):
        """
        One-time preprocessing pipeline.
        Skips if processed files already exist unless force=True.
        """
        if not force and CONCEPTS_TSV.exists() and RELATIONS_TSV.exists():
            print("[UMLS] Processed files already exist — skipping.")
            print(f"       Delete {PROCESSED_DIR} and rerun to force rebuild.")
            return

        print("[UMLS] Step 1/3 — Building semantic type index from MRSTY...")
        semtype_index = self._build_semtype_index()
        print(f"[UMLS]   -> {len(semtype_index):,} CUIs with target semantic types")

        print("[UMLS] Step 2/3 — Extracting concepts from MRCONSO...")
        n_concepts = self._extract_concepts(semtype_index)
        print(f"[UMLS]   -> {n_concepts:,} concept rows written to {CONCEPTS_TSV.name}")

        print("[UMLS] Step 3/3 — Building known CUI set for relation filtering...")
        known_cuis = set()
        with open(CONCEPTS_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                known_cuis.add(row["cui"])
        print(f"[UMLS]   -> {len(known_cuis):,} unique CUIs")

        print("[UMLS] Step 3/3 — Extracting relations from MRREL...")
        n_relations = self._extract_relations(known_cuis)
        print(f"[UMLS]   -> {n_relations:,} relations written to {RELATIONS_TSV.name}")
        print("[UMLS] Preprocessing complete.")

    def nodes(self) -> Generator[Node, None, None]:
        """Yield Node objects from the preprocessed concepts TSV."""
        if not CONCEPTS_TSV.exists():
            raise FileNotFoundError(
                f"{CONCEPTS_TSV} not found. Run UMLSSource().preprocess() first."
            )
        current_cui  = None
        current_rows = []

        with open(CONCEPTS_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row["cui"] != current_cui:
                    if current_rows:
                        node = self._rows_to_node(current_rows)
                        if node:
                            yield node
                    current_cui  = row["cui"]
                    current_rows = [row]
                else:
                    current_rows.append(row)
            if current_rows:
                node = self._rows_to_node(current_rows)
                if node:
                    yield node

    def edges(self) -> Generator[Edge, None, None]:
        """Yield Edge objects from the preprocessed relations TSV."""
        if not RELATIONS_TSV.exists():
            raise FileNotFoundError(
                f"{RELATIONS_TSV} not found. Run UMLSSource().preprocess() first."
            )
        with open(RELATIONS_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                edge = self._row_to_edge(row)
                if edge:
                    yield edge

    def normalize_confidence(self, raw_value=None) -> float:
        return CONFIDENCE_RULES["UMLS_MRREL"]

    # ── Preprocessing passes ───────────────────────────────────────────

    def _build_semtype_index(self) -> dict:
        """Pass 1: MRSTY -> {CUI: {tui, sty}} for target semtypes only."""
        index = {}
        with open(UMLS_MRSTY, encoding="utf-8") as f:
            for i, line in enumerate(f):
                parts = line.rstrip("\n").split("|")
                if len(parts) < 4:
                    continue
                cui = parts[T_CUI]
                tui = parts[T_TUI]
                sty = parts[T_STY]
                if tui in UMLS_TARGET_SEMTYPES:
                    if cui not in index:
                        index[cui] = {"tui": tui, "sty": sty}
                if i % 1_000_000 == 0 and i > 0:
                    print(f"    MRSTY: {i:,} lines, {len(index):,} matching CUIs...")
        return index

    def _extract_concepts(self, semtype_index: dict) -> int:
        """
        Pass 2: MRCONSO -> concepts TSV.
        Filters: English + not obsolete + target SAB + target semtype CUI.
        """
        n_written = 0
        with open(UMLS_MRCONSO, encoding="utf-8") as fin, \
             open(CONCEPTS_TSV, "w", encoding="utf-8", newline="") as fout:
            writer = csv.writer(fout, delimiter="\t")
            writer.writerow(["cui", "tui", "sty", "sab", "code", "label", "ispref"])
            for i, line in enumerate(fin):
                parts = line.rstrip("\n").split("|")
                if len(parts) < 17:
                    continue
                cui      = parts[C_CUI]
                lat      = parts[C_LAT]
                ispref   = parts[C_ISPREF]
                sab      = parts[C_SAB]
                code     = parts[C_CODE]
                label    = parts[C_STR]
                suppress = parts[C_SUPPRESS]

                if lat      != "ENG":                continue
                if suppress == "O":                  continue
                if sab      not in UMLS_TARGET_SABS: continue
                if cui      not in semtype_index:    continue

                sem = semtype_index[cui]
                writer.writerow([
                    cui, sem["tui"], sem["sty"],
                    sab, code, label, ispref
                ])
                n_written += 1
                if i % 1_000_000 == 0 and i > 0:
                    print(f"    MRCONSO: {i:,} lines read, {n_written:,} rows written...")
        return n_written

    def _extract_relations(self, known_cuis: set) -> int:
        """
        Pass 3: MRREL -> relations TSV.
        Filters: target SABs + both CUIs known + not obsolete + has RELA.
        """
        n_written = 0
        with open(UMLS_MRREL, encoding="utf-8") as fin, \
             open(RELATIONS_TSV, "w", encoding="utf-8", newline="") as fout:
            writer = csv.writer(fout, delimiter="\t")
            writer.writerow(["cui1", "rel", "cui2", "rela", "sab"])
            for i, line in enumerate(fin):
                parts = line.rstrip("\n").split("|")
                if len(parts) < 15:
                    continue
                cui1     = parts[R_CUI1]
                rel      = parts[R_REL]
                cui2     = parts[R_CUI2]
                rela     = parts[R_RELA]
                sab      = parts[R_SAB]
                suppress = parts[R_SUPPRESS]

                if suppress == "O":                  continue
                if sab not in UMLS_TARGET_SABS:      continue
                if cui1 not in known_cuis:           continue
                if cui2 not in known_cuis:           continue
                if not rela and rel not in ("CHD","PAR","RB","RN"): continue

                writer.writerow([cui1, rel, cui2, rela or rel, sab])
                n_written += 1
                if i % 2_000_000 == 0 and i > 0:
                    print(f"    MRREL: {i:,} lines read, {n_written:,} rows written...")
        return n_written

    # ── Node construction ──────────────────────────────────────────────

    def _rows_to_node(self, rows: list) -> Node | None:
        """
        Convert all MRCONSO rows for one CUI into a single Node.
        Picks best label, collects all SAB codes as xrefs,
        assigns primary_id based on semantic type preference.
        """
        cui = rows[0]["cui"]
        tui = rows[0]["tui"]
        sty = rows[0]["sty"]

        entity_info = SEMTYPE_TO_ENTITY.get(tui)
        if not entity_info:
            return None
        entity_type, tier, preferred_system = entity_info

        codes    = {}
        label    = None
        synonyms = []

        for row in rows:
            sab    = row["sab"]
            code   = row["code"]
            text   = row["label"]
            ispref = row["ispref"]
            system = SAB_TO_SYSTEM.get(sab, sab)

            if system not in codes:
                codes[system] = code

            if ispref == "Y" and label is None:
                label = text
            elif text != label and text not in synonyms:
                synonyms.append(text)

        if not label:
            label = rows[0]["label"]

        # Pick primary_id from preferred system, fall back to first available
        primary_id     = codes.get(preferred_system)
        primary_system = preferred_system
        if not primary_id:
            for sys, code in codes.items():
                primary_id     = code
                primary_system = sys
                break
        if not primary_id:
            return None

        xrefs = {"UMLS_CUI": cui}
        xrefs.update({
            sys: code for sys, code in codes.items()
            if sys != primary_system
        })

        return Node(
            primary_id     = primary_id,
            primary_system = primary_system,
            label          = label[:500],
            tier           = tier,
            entity_type    = entity_type,
            xrefs          = xrefs,
            icd10_code     = codes.get("ICD-10-CM"),
            snomed_code    = codes.get("SNOMED"),
            loinc_code     = codes.get("LOINC"),
            rxnorm_cui     = codes.get("RxNorm"),
            cpt_code       = codes.get("CPT"),
            synonyms       = synonyms[:20],
            properties     = {"umls_sty": sty, "umls_tui": tui},
            source         = self.source_name,
            source_version = self.source_version,
            confidence     = 0.9,
        )

    # ── Edge construction ──────────────────────────────────────────────

    def _row_to_edge(self, row: dict) -> Edge | None:
        """Convert one MRREL row into an Edge."""
        rela     = row["rela"]
        sab      = row["sab"]
        rel_type = UMLS_REL_MAP.get(rela, "ASSOCIATED_WITH")

        try:
            return Edge(
                source_id                = row["cui1"],
                source_system            = "UMLS_CUI",
                target_id                = row["cui2"],
                target_system            = "UMLS_CUI",
                relationship_type        = rel_type,
                source_relationship_type = rela,
                confidence               = self.normalize_confidence(),
                primary_source           = f"UMLS_{self.source_version}_{sab}",
                imported_via             = f"UMLS_MRREL_{self.source_version}",
                study_design             = "curated",
                source_version           = self.source_version,
            )
        except ValueError:
            return None
