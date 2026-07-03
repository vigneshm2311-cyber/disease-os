"""
config.py
Single source of truth for all paths and constants.
"""
import os
from pathlib import Path

ROOT     = Path(os.environ.get("DISEASE_OS_ROOT", Path.home() / "disease-os"))
DATA_DIR = ROOT / "data"
DB_PATH  = DATA_DIR / "diseaseos.db"

# processed intermediate files (filtered subsets extracted from raw RRFs)
PROCESSED_DIR = DATA_DIR / "processed"

# ── UMLS ──────────────────────────────────────────────────────────────
# Points at your local 2025AB release.
# Change UMLS_DIR if you move the files.
UMLS_DIR     = Path.home() / "Downloads" / "2025AB" / "META"
UMLS_MRCONSO = UMLS_DIR / "MRCONSO.RRF"
UMLS_MRREL   = UMLS_DIR / "MRREL.RRF"
UMLS_MRSTY   = UMLS_DIR / "MRSTY.RRF"
UMLS_MRDEF   = UMLS_DIR / "MRDEF.RRF"
UMLS_MRSAT   = UMLS_DIR / "MRSAT.RRF"

# ── Source versions ────────────────────────────────────────────────────
SOURCE_VERSIONS = {
    "UMLS":         "2025AB",
    "REACTOME":     "88",
    "UNIPROT":      "2024_03",
    "CLINVAR":      "2024-06",
    "STRING":       "12.0",
    "GWAS_CATALOG": "2024-06",
    "HMDB":         "5.0",
    "KEGG":         "2024-06",
    "CHEBI":        "231",
    "DRUGBANK":     "5.1.12",
    "HPO":          "2024-04",
    "MONDO":        "2024-05",
}

# ── UMLS streaming config ──────────────────────────────────────────────
# Source vocabularies we actually want from MRCONSO/MRREL.
# Filtering to these SABs cuts the 2.1GB MRCONSO down ~90%.
# MRREL (5.7GB) filtered similarly.
# MRSAT (8.9GB) — we only extract specific ATN attributes, never full load.
UMLS_TARGET_SABS = {
    "SNOMEDCT_US",  # clinical findings, anatomy, procedures
    "ICD10CM",      # diagnosis billing
    "ICD11",        # forward compatibility
    "LNC",          # LOINC — lab tests
    "RXNORM",       # drugs
    "HPO",          # phenotypes
    "MONDO",        # disease concepts
    "OMIM",         # gene-disease
    "MSH",          # MeSH — broad biomedical
    "NCI",          # NCI thesaurus
    "MEDLINEPLUS",  # patient-facing definitions
    "CPT",          # procedures
    "GO",           # gene ontology (partial in UMLS)
    "FMA",          # anatomy
}

# Semantic types (from MRSTY) we want to keep.
# Everything else filtered out — cuts noise significantly.
UMLS_TARGET_SEMTYPES = {
    "T047",  # Disease or Syndrome
    "T048",  # Mental or Behavioral Dysfunction
    "T191",  # Neoplastic Process
    "T046",  # Pathologic Function
    "T121",  # Pharmacologic Substance
    "T116",  # Amino Acid, Peptide, or Protein
    "T028",  # Gene or Genome
    "T086",  # Nucleotide Sequence
    "T059",  # Laboratory Procedure
    "T034",  # Laboratory or Test Result
    "T023",  # Body Part, Organ, or Organ Component
    "T025",  # Cell
    "T043",  # Cell Function
    "T044",  # Molecular Function
    "T045",  # Genetic Function
    "T038",  # Biologic Function
    "T031",  # Body Substance
    "T109",  # Organic Chemical
    "T123",  # Biologically Active Substance
    "T033",  # Finding
    "T184",  # Sign or Symptom
    "T201",  # Clinical Attribute
    "T200",  # Clinical Drug
    "T074",  # Medical Device
    "T058",  # Health Care Activity
}

# MRSAT attributes we extract (everything else in 8.9GB MRSAT ignored)
UMLS_TARGET_ATN = {
    "ICD10CM",       # ICD-10-CM code on a concept
    "ICD10",         # ICD-10 code
    "RXCUI",         # RxNorm CUI
    "LNC",           # LOINC code
    "CPT",           # CPT code
    "HCPCS",         # HCPCS code
    "NUI",           # NCI code
    "HCC",           # HCC risk adjustment code (if present)
    "OMIM_NUMBER",   # OMIM ID
    "MONDO_ID",      # MONDO ID
    "HPO",           # HPO ID
}

# ── Confidence normalization ───────────────────────────────────────────
CONFIDENCE_RULES = {
    "CLINVAR_STARS": {
        0: 0.30,
        1: 0.55,
        2: 0.70,
        3: 0.85,
        4: 0.95,
    },
    "STRING_SCORE":      lambda s: round(min(s / 1000, 1.0), 3),
    "GWAS_PVALUE":       lambda p: 0.65 if p < 5e-8 else 0.50 if p < 1e-5 else 0.35,
    "REACTOME_CURATED":  0.85,
    "UNIPROT_SWISSPROT": 0.90,
    "UMLS_MRREL":        0.60,
    "LITERATURE":        0.55,
}

# ── Primary ID systems per entity type ────────────────────────────────
PRIMARY_ID_SYSTEMS = {
    "Disease_clinical":   "ICD-10-CM",
    "Disease_scientific": "MONDO",
    "ClinicalFinding":    "SNOMED",
    "LabTest":            "LOINC",
    "Drug_clinical":      "RxNorm",
    "Drug_chemical":      "ChEBI",
    "Gene":               "NCBI_Gene",
    "Protein":            "UniProt",
    "Variant":            "dbSNP_rsID",
    "Metabolite":         "HMDB",
    "Pathway":            "Reactome",
    "CellType":           "CellOntology",
    "Anatomy":            "UBERON",
    "Phenotype":          "HPO",
    "Procedure":          "CPT",
    "SocialDeterminant":  "ICD-10-Z",
}
