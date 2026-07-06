"""
core/sources/hmdb.py
HMDB v5.0 adapter for Disease-OS.
Namespace: http://www.hmdb.ca — stripped at parse time.
"""

import csv, sys, json
from pathlib import Path
from typing import Generator
from xml.etree import ElementTree as ET

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS, PROCESSED_DIR

RAW_DIR  = Path.home() / "disease-os" / "data" / "raw" / "hmdb"
HMDB_XML = RAW_DIR / "hmdb_metabolites.xml"
HMDB_TSV = PROCESSED_DIR / "hmdb_metabolites.tsv"

TARGET_STATUS = {"quantified", "detected", "expected"}

CLINICAL_BIOFLUIDS = {
    "blood", "plasma", "serum", "urine", "cerebrospinal fluid",
    "csf", "saliva", "feces", "bile", "sweat", "breast milk",
    "amniotic fluid", "aqueous humour",
}

SKIP_DISEASE_TERMS = {
    "physiological effect", "health effect", "health condition",
    "disposition", "process", "role", "biological function",
    "industrial application", "naturally occurring process",
    "biological process", "biochemical pathway", "source",
    "biological location", "tissue and substructures",
    "biofluid and excreta", "subcellular",
}

NS = "{http://www.hmdb.ca}"


def _t(tag: str) -> str:
    """Add namespace prefix."""
    return f"{NS}{tag}"


def _text(elem, tag: str) -> str:
    node = elem.find(_t(tag))
    return node.text.strip() if node is not None and node.text else ""


def _parse_metabolite(elem) -> dict | None:
    """
    Parse one <metabolite> element.
    Tags still have namespace prefix — use _t() helper.
    """
    status = _text(elem, "status").lower()
    if status not in TARGET_STATUS:
        return None

    accession = _text(elem, "accession")
    if not accession:
        return None

    # Synonyms
    synonyms = [
        s.text.strip()
        for s in elem.findall(f"{_t('synonyms')}/{_t('synonym')}")
        if s.text and s.text.strip()
    ]

    # Cross-references
    xrefs = {}
    for xref_tag, xref_key in [
        ("chebi_id",          "ChEBI"),
        ("kegg_id",           "KEGG"),
        ("pubchem_compound_id","PubChem"),
        ("drugbank_id",        "DrugBank"),
        ("cas_registry_number","CAS"),
        ("inchikey",           "InChIKey"),
    ]:
        val = _text(elem, xref_tag)
        if val:
            xrefs[xref_key] = val

    # Diseases — from dedicated <diseases><disease><name> section
    diseases = []
    diseases_elem = elem.find(_t("diseases"))
    if diseases_elem is not None:
        for d in diseases_elem.findall(_t("disease")):
            name_elem = d.find(_t("name"))
            if name_elem is not None and name_elem.text:
                name = name_elem.text.strip()
                if name and name.lower() not in SKIP_DISEASE_TERMS:
                    diseases.append(name)

    # Biofluids — from <normal_concentrations><concentration><biofluid>
    biofluids = set()
    normal_elem = elem.find(_t("normal_concentrations"))
    if normal_elem is not None:
        for c in normal_elem.findall(_t("concentration")):
            bf_elem = c.find(_t("biofluid"))
            if bf_elem is not None and bf_elem.text:
                bf = bf_elem.text.strip().lower()
                if bf in CLINICAL_BIOFLUIDS:
                    biofluids.add(bf)

    # Also check abnormal concentrations for additional biofluids
    abnormal_elem = elem.find(_t("abnormal_concentrations"))
    if abnormal_elem is not None:
        for c in abnormal_elem.findall(_t("concentration")):
            bf_elem = c.find(_t("biofluid"))
            if bf_elem is not None and bf_elem.text:
                bf = bf_elem.text.strip().lower()
                if bf in CLINICAL_BIOFLUIDS:
                    biofluids.add(bf)

    # Pathways — from <biological_properties><pathways><pathway><name>
    pathways = []
    bio_elem = elem.find(_t("biological_properties"))
    if bio_elem is not None:
        pathways_elem = bio_elem.find(_t("pathways"))
        if pathways_elem is not None:
            for p in pathways_elem.findall(_t("pathway")):
                pname = p.find(_t("name"))
                if pname is not None and pname.text:
                    pathways.append(pname.text.strip())

    # Protein associations
    proteins = []
    prot_elem = elem.find(_t("protein_associations"))
    if prot_elem is not None:
        for p in prot_elem.findall(_t("protein")):
            uniprot = p.find(_t("uniprot_id"))
            if uniprot is not None and uniprot.text:
                proteins.append(uniprot.text.strip())

    return {
        "accession":   accession,
        "name":        _text(elem, "name"),
        "description": _text(elem, "description")[:500],
        "formula":     _text(elem, "chemical_formula"),
        "mass":        _text(elem, "average_molecular_weight"),
        "smiles":      _text(elem, "smiles"),
        "status":      status,
        "xrefs":       xrefs,
        "synonyms":    synonyms[:15],
        "diseases":    list(set(diseases))[:20],
        "biofluids":   list(biofluids)[:10],
        "pathways":    list(set(pathways))[:10],
        "proteins":    proteins[:20],
    }


class HMDBSource(BaseSource):

    source_name    = "HMDB"
    source_version = SOURCE_VERSIONS["HMDB"]

    def __init__(self):
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    def preprocess(self, force: bool = False):
        if not force and HMDB_TSV.exists():
            print("[HMDB] Processed file exists — skipping.")
            return

        size_gb = HMDB_XML.stat().st_size / 1_000_000_000
        print(f"[HMDB] Preprocessing {HMDB_XML.name} ({size_gb:.1f}GB)...")

        n_read = n_written = n_skipped = 0

        with open(HMDB_TSV, "w", encoding="utf-8", newline="") as fout:
            writer = csv.writer(fout, delimiter="\t")
            writer.writerow([
                "accession", "name", "description", "formula",
                "mass", "smiles", "status",
                "xrefs_json", "synonyms_json", "diseases_json",
                "biofluids_json", "pathways_json", "proteins_json",
            ])

            # Use iterparse on "end" event for metabolite
            # Tags keep their namespace — _t() helper handles this
            context = ET.iterparse(str(HMDB_XML), events=("end",))

            for event, elem in context:
                if elem.tag != f"{NS}metabolite":
                    continue

                n_read += 1
                m = _parse_metabolite(elem)
                elem.clear()

                if m is None:
                    n_skipped += 1
                    continue

                writer.writerow([
                    m["accession"], m["name"], m["description"],
                    m["formula"], m["mass"], m["smiles"], m["status"],
                    json.dumps(m["xrefs"]),
                    json.dumps(m["synonyms"]),
                    json.dumps(m["diseases"]),
                    json.dumps(m["biofluids"]),
                    json.dumps(m["pathways"]),
                    json.dumps(m["proteins"]),
                ])
                n_written += 1

                if n_written % 10_000 == 0:
                    print(f"  {n_read:,} read | {n_written:,} written | "
                          f"{n_skipped:,} skipped...")

        print(f"\n[HMDB] Done: {n_read:,} read | "
              f"{n_written:,} written | {n_skipped:,} skipped")

    def nodes(self) -> Generator[Node, None, None]:
        if not HMDB_TSV.exists():
            raise FileNotFoundError("Run HMDBSource().preprocess() first.")

        with open(HMDB_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if not row["accession"] or not row["name"]:
                    continue
                try:
                    mass = float(row["mass"]) if row["mass"] else None
                except ValueError:
                    mass = None

                xrefs = json.loads(row["xrefs_json"] or "{}")

                yield Node(
                    primary_id     = row["accession"],
                    primary_system = "HMDB",
                    label          = row["name"],
                    tier           = 1,
                    entity_type    = "Metabolite",
                    xrefs          = xrefs,
                    synonyms       = json.loads(row["synonyms_json"] or "[]"),
                    definition     = row["description"] or None,
                    properties     = {k: v for k, v in {
                        "formula":   row["formula"],
                        "mass_da":   mass,
                        "smiles":    row["smiles"],
                        "status":    row["status"],
                        "biofluids": json.loads(row["biofluids_json"] or "[]"),
                    }.items() if v},
                    source         = self.source_name,
                    source_version = self.source_version,
                    confidence     = 0.9 if row["status"] == "quantified" else 0.75,
                )

    def edges(self) -> Generator[Edge, None, None]:
        if not HMDB_TSV.exists():
            raise FileNotFoundError("Run HMDBSource().preprocess() first.")

        with open(HMDB_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                accession = row["accession"]
                if not accession:
                    continue

                # DETECTED_IN edges for biofluids
                for bf in json.loads(row["biofluids_json"] or "[]"):
                    try:
                        yield Edge(
                            source_id                = accession,
                            source_system            = "HMDB",
                            target_id                = bf,
                            target_system            = "HMDB_biofluid",
                            relationship_type        = "DETECTED_IN",
                            source_relationship_type = "hmdb_biofluid",
                            confidence               = 0.85,
                            primary_source           = f"HMDB_{self.source_version}",
                            imported_via             = f"HMDB_{self.source_version}",
                            study_design             = "curated",
                            source_version           = self.source_version,
                        )
                    except ValueError:
                        continue

                # ASSOCIATED_WITH edges for diseases
                for disease in json.loads(row["diseases_json"] or "[]"):
                    try:
                        yield Edge(
                            source_id                = accession,
                            source_system            = "HMDB",
                            target_id                = disease,
                            target_system            = "HMDB_disease_label",
                            relationship_type        = "ASSOCIATED_WITH",
                            source_relationship_type = "hmdb_disease",
                            confidence               = 0.65,
                            primary_source           = f"HMDB_{self.source_version}",
                            imported_via             = f"HMDB_{self.source_version}",
                            study_design             = "curated",
                            source_version           = self.source_version,
                        )
                    except ValueError:
                        continue

                # PART_OF edges for pathways (links to Reactome where names match)
                for pathway in json.loads(row["pathways_json"] or "[]"):
                    try:
                        yield Edge(
                            source_id                = accession,
                            source_system            = "HMDB",
                            target_id                = pathway,
                            target_system            = "HMDB_pathway_label",
                            relationship_type        = "PART_OF",
                            source_relationship_type = "hmdb_pathway",
                            confidence               = 0.80,
                            primary_source           = f"HMDB_{self.source_version}",
                            imported_via             = f"HMDB_{self.source_version}",
                            study_design             = "curated",
                            source_version           = self.source_version,
                        )
                    except ValueError:
                        continue

                # INTERACTS_WITH edges for protein associations
                for uniprot_id in json.loads(row["proteins_json"] or "[]"):
                    try:
                        yield Edge(
                            source_id                = accession,
                            source_system            = "HMDB",
                            target_id                = uniprot_id,
                            target_system            = "UniProt",
                            relationship_type        = "INTERACTS_WITH",
                            source_relationship_type = "hmdb_protein_association",
                            confidence               = 0.75,
                            primary_source           = f"HMDB_{self.source_version}",
                            imported_via             = f"HMDB_{self.source_version}",
                            study_design             = "curated",
                            source_version           = self.source_version,
                        )
                    except ValueError:
                        continue

    def normalize_confidence(self, raw_value=None) -> float:
        return 0.85
