"""
core/sources/clinvar.py

ClinVar adapter for Disease-OS.

Input : data/raw/clinvar/variant_summary.txt (3.7GB)
Output: two compact TSVs in data/processed/
  clinvar_variants.tsv      — one row per variant  -> Variant nodes
  clinvar_var_disease.tsv   — variant-disease pairs -> edges

Filters:
  Assembly = GRCh38 only
  Origin   = germline / inherited
  GeneID   != -1
  ClinSig  in target set
  Phenotype not "not provided" / "not specified"

ReviewStatus -> confidence:
  practice guideline                                   -> 0.95
  reviewed by expert panel                             -> 0.85
  criteria provided, multiple submitters, no conflicts -> 0.70
  criteria provided, single submitter                  -> 0.55
  criteria provided, conflicting interpretations       -> 0.40
  no assertion criteria provided                       -> 0.35
  no assertion provided                                -> 0.30

ClinicalSignificance -> canonical edge type:
  Pathogenic / Pathogenic/Likely pathogenic -> CAUSES
  Likely pathogenic                         -> CONTRIBUTES_TO
  risk factor / likely risk allele          -> INCREASES_RISK_OF
  protective / benign / likely benign       -> PROTECTS_AGAINST
  uncertain significance / conflicting      -> ASSOCIATED_WITH
"""

import csv
import json
import sys
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS, PROCESSED_DIR

RAW_DIR         = Path.home() / "disease-os" / "data" / "raw" / "clinvar"
VARIANT_SUMMARY = RAW_DIR / "variant_summary.txt"
VARIANTS_TSV    = PROCESSED_DIR / "clinvar_variants.tsv"
VAR_DISEASE_TSV = PROCESSED_DIR / "clinvar_var_disease.tsv"

# Column indices (0-based)
C_ALLELE_ID  = 0
C_TYPE       = 1
C_NAME       = 2
C_GENE_ID    = 3
C_GENE_SYM   = 4
C_CLIN_SIG   = 6
C_RS         = 9
C_PHENO_IDS  = 12
C_PHENO_LIST = 13
C_ORIGIN     = 14
C_ASSEMBLY   = 16
C_CHROM      = 18
C_START      = 19
C_REF        = 21
C_ALT        = 22
C_REVIEW     = 24
C_N_SUBMIT   = 25
C_VAR_ID     = 30

CLINSIG_TO_REL = {
    "pathogenic":                                          "CAUSES",
    "likely pathogenic":                                   "CONTRIBUTES_TO",
    "pathogenic/likely pathogenic":                        "CAUSES",
    "pathogenic, risk factor":                             "CAUSES",
    "risk factor":                                         "INCREASES_RISK_OF",
    "likely risk allele":                                  "INCREASES_RISK_OF",
    "protective":                                          "PROTECTS_AGAINST",
    "benign":                                              "PROTECTS_AGAINST",
    "likely benign":                                       "PROTECTS_AGAINST",
    "uncertain significance":                              "ASSOCIATED_WITH",
    "conflicting interpretations of pathogenicity":        "ASSOCIATED_WITH",
}

REVIEW_CONFIDENCE = {
    "practice guideline":                                    0.95,
    "reviewed by expert panel":                              0.85,
    "criteria provided, multiple submitters, no conflicts":  0.70,
    "criteria provided, single submitter":                   0.55,
    "criteria provided, conflicting interpretations":        0.40,
    "no assertion criteria provided":                        0.35,
    "no assertion provided":                                 0.30,
    "no classifications from unflagged records":             0.25,
}

TARGET_CLINSIG = set(CLINSIG_TO_REL.keys())

SKIP_PHENOTYPES = {
    "not provided", "not specified", "see cases",
    "allelic variant", "variant of unknown significance",
}


def _parse_phenotype_ids(pheno_ids_str: str) -> dict:
    result = {}
    if not pheno_ids_str or pheno_ids_str == "-":
        return result
    for part in pheno_ids_str.split(","):
        part = part.strip()
        if ":" not in part:
            continue
        system, _, code = part.partition(":")
        system = system.strip()
        code   = code.strip()
        if system == "MedGen":
            result["UMLS_CUI"] = code
        elif system == "OMIM":
            result["OMIM"] = code
        elif system == "SNOMED CT":
            result["SNOMED"] = code
        elif system == "Orphanet":
            result["Orphanet"] = code
        elif system == "Human Phenotype Ontology":
            result["HPO"] = code
    return result


class ClinVarSource(BaseSource):

    source_name    = "CLINVAR"
    source_version = SOURCE_VERSIONS["CLINVAR"]

    def __init__(self):
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    def preprocess(self, force: bool = False):
        if not force and VARIANTS_TSV.exists() and VAR_DISEASE_TSV.exists():
            print("[CLINVAR] Processed files exist — skipping.")
            print(f"          Delete {PROCESSED_DIR}/clinvar_* to force rebuild.")
            return

        if not VARIANT_SUMMARY.exists():
            raise FileNotFoundError(f"Not found: {VARIANT_SUMMARY}")

        size_mb = VARIANT_SUMMARY.stat().st_size // 1_000_000
        print(f"[CLINVAR] Preprocessing {VARIANT_SUMMARY.name} ({size_mb}MB)...")
        print("[CLINVAR] Filters: GRCh38 + germline + pathogenic/risk/benign")

        n_read = n_variants = n_pairs = n_skipped = 0

        with open(VARIANT_SUMMARY, encoding="utf-8") as fin, \
             open(VARIANTS_TSV,    "w", encoding="utf-8", newline="") as fv, \
             open(VAR_DISEASE_TSV, "w", encoding="utf-8", newline="") as fd:

            var_writer  = csv.writer(fv, delimiter="\t")
            pair_writer = csv.writer(fd, delimiter="\t")

            var_writer.writerow([
                "primary_id", "primary_system", "allele_id", "variation_id",
                "variant_type", "hgvs_name", "gene_id", "gene_symbol",
                "chromosome", "start", "ref", "alt",
                "clinsig", "review_status", "n_submitters",
            ])
            pair_writer.writerow([
                "variant_primary_id", "variant_system",
                "pheno_ids_json", "pheno_name",
                "clinsig", "review_status", "rel_type", "confidence",
                "gene_id", "gene_symbol",
            ])

            reader = csv.reader(fin, delimiter="\t")
            next(reader)   # skip header

            for row in reader:
                n_read += 1
                if n_read % 500_000 == 0:
                    print(f"  {n_read:,} rows | "
                          f"{n_variants:,} variants | "
                          f"{n_pairs:,} pairs")

                if len(row) < 31:
                    n_skipped += 1
                    continue

                # ── Filters ────────────────────────────────────────────
                if row[C_ASSEMBLY].strip() != "GRCh38":
                    n_skipped += 1
                    continue

                origin = row[C_ORIGIN].strip().lower()
                if "germline" not in origin and "inherited" not in origin:
                    n_skipped += 1
                    continue

                gene_id = row[C_GENE_ID].strip()
                if gene_id in ("-1", "", "-"):
                    n_skipped += 1
                    continue

                raw_clinsig = row[C_CLIN_SIG].strip().lower()
                clinsig     = None
                if raw_clinsig in TARGET_CLINSIG:
                    clinsig = raw_clinsig
                else:
                    for target in TARGET_CLINSIG:
                        if target in raw_clinsig:
                            clinsig = target
                            break
                if clinsig is None:
                    n_skipped += 1
                    continue

                pheno_list = row[C_PHENO_LIST].strip().lower()
                if any(skip in pheno_list for skip in SKIP_PHENOTYPES):
                    n_skipped += 1
                    continue

                # ── Primary ID ─────────────────────────────────────────
                rs_raw = row[C_RS].strip()
                if rs_raw and rs_raw not in ("-1", "-", ""):
                    primary_id     = f"rs{rs_raw}"
                    primary_system = "dbSNP_rsID"
                else:
                    primary_id     = f"CA{row[C_ALLELE_ID].strip()}"
                    primary_system = "ClinVar_AlleleID"

                # ── Variant row ────────────────────────────────────────
                var_writer.writerow([
                    primary_id, primary_system,
                    row[C_ALLELE_ID].strip(),
                    row[C_VAR_ID].strip(),
                    row[C_TYPE].strip(),
                    row[C_NAME].strip()[:300],
                    gene_id,
                    row[C_GENE_SYM].strip(),
                    row[C_CHROM].strip(),
                    row[C_START].strip(),
                    row[C_REF].strip()[:50],
                    row[C_ALT].strip()[:50],
                    clinsig,
                    row[C_REVIEW].strip(),
                    row[C_N_SUBMIT].strip(),
                ])
                n_variants += 1

                # ── Disease pair rows ──────────────────────────────────
                review_status = row[C_REVIEW].strip().lower()
                confidence    = REVIEW_CONFIDENCE.get(review_status, 0.35)
                rel_type      = CLINSIG_TO_REL.get(clinsig, "ASSOCIATED_WITH")
                pheno_ids     = _parse_phenotype_ids(row[C_PHENO_IDS].strip())

                if not pheno_ids:
                    continue

                for pheno_name in row[C_PHENO_LIST].strip().split("|"):
                    pheno_name = pheno_name.strip()
                    if not pheno_name or pheno_name.lower() in SKIP_PHENOTYPES:
                        continue
                    pair_writer.writerow([
                        primary_id, primary_system,
                        json.dumps(pheno_ids),
                        pheno_name[:200],
                        clinsig, review_status,
                        rel_type, confidence,
                        gene_id, row[C_GENE_SYM].strip(),
                    ])
                    n_pairs += 1

        print(f"\n[CLINVAR] Done:")
        print(f"  Rows read    : {n_read:,}")
        print(f"  Rows skipped : {n_skipped:,}")
        print(f"  Variants     : {n_variants:,}  -> {VARIANTS_TSV.name}")
        print(f"  Disease pairs: {n_pairs:,}  -> {VAR_DISEASE_TSV.name}")

    def nodes(self) -> Generator[Node, None, None]:
        if not VARIANTS_TSV.exists():
            raise FileNotFoundError("Run ClinVarSource().preprocess() first.")
        seen = set()
        with open(VARIANTS_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                pk = (row["primary_id"], row["primary_system"])
                if pk in seen:
                    continue
                seen.add(pk)
                yield Node(
                    primary_id     = row["primary_id"],
                    primary_system = row["primary_system"],
                    label          = row["hgvs_name"][:200] or row["primary_id"],
                    tier           = 1,
                    entity_type    = "Variant",
                    xrefs          = {
                        "ClinVar_AlleleID": row["allele_id"],
                        "ClinVar_VarID":    row["variation_id"],
                    },
                    properties     = {
                        "variant_type":  row["variant_type"],
                        "gene_id":       row["gene_id"],
                        "gene_symbol":   row["gene_symbol"],
                        "chromosome":    row["chromosome"],
                        "start":         row["start"],
                        "ref":           row["ref"],
                        "alt":           row["alt"],
                        "clinsig":       row["clinsig"],
                        "review_status": row["review_status"],
                        "n_submitters":  row["n_submitters"],
                    },
                    source         = self.source_name,
                    source_version = self.source_version,
                    confidence     = REVIEW_CONFIDENCE.get(
                        row["review_status"].lower(), 0.35
                    ),
                )

    def edges(self) -> Generator[Edge, None, None]:
        if not VAR_DISEASE_TSV.exists():
            raise FileNotFoundError("Run ClinVarSource().preprocess() first.")
        with open(VAR_DISEASE_TSV, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                pheno_ids = json.loads(row["pheno_ids_json"])
                # Target priority: SNOMED -> UMLS_CUI -> OMIM
                target_id = target_system = None
                if "SNOMED" in pheno_ids:
                    target_id, target_system = pheno_ids["SNOMED"], "SNOMED"
                elif "UMLS_CUI" in pheno_ids:
                    target_id, target_system = pheno_ids["UMLS_CUI"], "UMLS_CUI"
                elif "OMIM" in pheno_ids:
                    target_id, target_system = pheno_ids["OMIM"], "OMIM"
                if not target_id:
                    continue
                try:
                    yield Edge(
                        source_id                = row["variant_primary_id"],
                        source_system            = row["variant_system"],
                        target_id                = target_id,
                        target_system            = target_system,
                        relationship_type        = row["rel_type"],
                        source_relationship_type = row["clinsig"],
                        confidence               = float(row["confidence"]),
                        primary_source           = f"ClinVar_{self.source_version}",
                        imported_via             = f"ClinVar_variant_summary_{self.source_version}",
                        study_design             = "clinical_testing",
                        source_version           = self.source_version,
                    )
                except ValueError:
                    continue

    def normalize_confidence(self, raw_value=None) -> float:
        if isinstance(raw_value, str):
            return REVIEW_CONFIDENCE.get(raw_value.lower(), 0.35)
        return 0.35
