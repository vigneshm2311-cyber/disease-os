"""
core/sources/gwas_catalog.py

GWAS Catalog adapter for Disease-OS.

Input : data/raw/gwas_catalog/associations.tsv (554MB, ~400K rows)
        data/raw/gwas_catalog/studies.tsv (70MB)

What GWAS Catalog adds that ClinVar doesn't have:
  - Common variant associations (MAF > 1%) — the polygenic risk landscape
  - Effect sizes (OR/Beta) with confidence intervals
  - Population-stratified associations
  - p-values for statistical strength

Key design decisions:
  1. Only load genome-wide significant hits: p < 5e-8 (PVALUE_MLOG >= 7.3)
     Below this threshold, associations are likely false positives.
  2. Confidence = function of p-value, NOT causality.
     GWAS = ASSOCIATION not mechanism. All edges enter as INCREASES_RISK_OF
     never as CAUSES — that distinction is architecturally enforced here.
  3. Effect size from OR or BETA column — stored as effect_size.
     OR > 1 = risk allele increases risk (direction = positive)
     OR < 1 = risk allele decreases risk (direction = negative = protective)
  4. rsID from SNPS column (cleaned), not STRONGEST SNP-RISK ALLELE
     which contains the allele appended (e.g. "rs7903146-T").

p-value -> confidence mapping (GWAS association, not causation):
  p < 5e-30  (mlog >= 29.3) -> 0.65  (extremely strong signal)
  p < 5e-15  (mlog >= 14.3) -> 0.60
  p < 5e-8   (mlog >=  7.3) -> 0.55  (genome-wide significant floor)
  p >= 5e-8                  -> skip  (not genome-wide significant)

Note: confidence ceiling is 0.65 even for the strongest GWAS hits because
association strength ≠ causal effect size. This is documented in config.py.
"""

import csv
import sys
import re
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS, PROCESSED_DIR

RAW_DIR      = Path.home() / "disease-os" / "data" / "raw" / "gwas_catalog"
ASSOC_FILE   = RAW_DIR / "associations.tsv"
STUDIES_FILE = RAW_DIR / "studies.tsv"
PROCESSED_GWAS = PROCESSED_DIR / "gwas_associations.tsv"

# Column indices (0-based)
C_PUBMED      = 1
C_TRAIT       = 7
C_MAPPED_GENE = 14
C_GENE_IDS    = 17
C_SNP_ALLELE  = 20
C_SNPS        = 21
C_MERGED      = 22
C_SNP_CURRENT = 23
C_RAF         = 26
C_PVALUE      = 27
C_PVALUE_MLOG = 28
C_OR_BETA     = 30
C_CI          = 31

# Genome-wide significance threshold
GWS_MLOG = 7.3   # corresponds to p < 5e-8

def _mlog_to_confidence(mlog: float) -> float:
    if mlog >= 29.3: return 0.65
    if mlog >= 14.3: return 0.60
    return 0.55

def _parse_rsid(snps_str: str, current_str: str) -> str | None:
    """
    Extract clean rsID from SNPS column.

    GWAS Catalog quirk: SNP_ID_CURRENT stores the numeric part only
    e.g. SNPS='rs7179075' but SNP_ID_CURRENT='7179075' (no rs prefix).

    Strategy:
      1. Try SNPS column first — usually has full rsID
      2. Fall back to 'rs' + SNP_ID_CURRENT if SNPS has no rsID
      3. Skip if neither yields a valid rsID
    """
    # Pass 1: SNPS column (col 21) — primary source
    for part in re.split(r'[;,\s]+', snps_str.strip()):
        part = part.strip()
        if part.startswith("rs") and part[2:].isdigit():
            return part

    # Pass 2: SNP_ID_CURRENT (col 23) stores numeric part only
    # Prepend "rs" to reconstruct the full rsID
    cur = current_str.strip()
    if cur and cur.isdigit():
        return f"rs{cur}"

    # Pass 3: SNP_ID_CURRENT might already have rs prefix in some releases
    if cur.startswith("rs") and cur[2:].isdigit():
        return cur

    return None

def _parse_or_beta(or_beta_str: str) -> float | None:
    """Parse OR or Beta value — returns float or None."""
    try:
        val = float(or_beta_str.strip())
        return val if val != 0 else None
    except (ValueError, AttributeError):
        return None


class GWASCatalogSource(BaseSource):
    """
    Loads GWAS Catalog associations into Disease-OS.

    All associations enter as INCREASES_RISK_OF edges — never CAUSES.
    This architectural constraint is enforced here, not downstream.

    Usage:
        source = GWASCatalogSource()
        source.preprocess()
        source.load_into(graph_store)
    """

    source_name    = "GWAS_CATALOG"
    source_version = SOURCE_VERSIONS["GWAS_CATALOG"]

    def __init__(self):
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        self._study_ancestry: dict[str, str] = {}

    def preprocess(self, force: bool = False):
        """
        One-time preprocessing:
        1. Load study ancestry from studies.tsv
        2. Stream associations.tsv, filter GWS hits, write processed TSV
        """
        if not force and PROCESSED_GWAS.exists():
            print("[GWAS] Processed file exists — skipping preprocessing.")
            return

        # Load study ancestry index
        self._load_study_ancestry()

        print(f"[GWAS] Preprocessing {ASSOC_FILE.name} "
              f"({ASSOC_FILE.stat().st_size // 1_000_000}MB)...")
        print(f"[GWAS] Filter: PVALUE_MLOG >= {GWS_MLOG} (p < 5e-8)")

        n_read = n_written = n_skipped = 0

        with open(ASSOC_FILE, encoding="utf-8") as fin, \
             open(PROCESSED_GWAS, "w", encoding="utf-8", newline="") as fout:

            writer = csv.writer(fout, delimiter="\t")
            writer.writerow([
                "rsid", "pubmed_id", "trait", "mapped_gene",
                "gene_ids", "risk_allele", "raf",
                "pvalue_mlog", "or_beta", "ci_text",
                "confidence", "ancestry",
            ])

            reader = csv.reader(fin, delimiter="\t")
            header = next(reader)

            for row in reader:
                n_read += 1
                if n_read % 100_000 == 0:
                    print(f"  {n_read:,} rows read | "
                          f"{n_written:,} GWS hits written")

                if len(row) < 32:
                    n_skipped += 1
                    continue

                # ── Filter: genome-wide significant only ───────────────
                try:
                    mlog = float(row[C_PVALUE_MLOG].strip())
                except (ValueError, IndexError):
                    n_skipped += 1
                    continue

                if mlog < GWS_MLOG:
                    n_skipped += 1
                    continue

                # ── Extract rsID ───────────────────────────────────────
                rsid = _parse_rsid(
                    row[C_SNPS].strip(),
                    row[C_SNP_CURRENT].strip() if len(row) > C_SNP_CURRENT else ""
                )
                if not rsid:
                    n_skipped += 1
                    continue

                # ── Extract trait ──────────────────────────────────────
                trait = row[C_TRAIT].strip()
                if not trait or trait == "NR":
                    n_skipped += 1
                    continue

                # ── Extract effect size ────────────────────────────────
                or_beta = _parse_or_beta(
                    row[C_OR_BETA] if len(row) > C_OR_BETA else ""
                )

                # ── Risk allele (strip rsID prefix e.g. "rs123-A" -> "A")
                risk_allele_raw = row[C_SNP_ALLELE].strip()
                risk_allele = risk_allele_raw.split("-")[-1] \
                    if "-" in risk_allele_raw else risk_allele_raw

                # ── Ancestry from study lookup ─────────────────────────
                pubmed_id = row[C_PUBMED].strip()
                ancestry  = self._study_ancestry.get(pubmed_id, "NR")

                confidence = _mlog_to_confidence(mlog)

                writer.writerow([
                    rsid,
                    pubmed_id,
                    trait,
                    row[C_MAPPED_GENE].strip(),
                    row[C_GENE_IDS].strip(),
                    risk_allele,
                    row[C_RAF].strip(),
                    mlog,
                    or_beta if or_beta is not None else "",
                    row[C_CI].strip() if len(row) > C_CI else "",
                    confidence,
                    ancestry,
                ])
                n_written += 1

        print(f"\n[GWAS] Preprocessing complete:")
        print(f"  Rows read    : {n_read:,}")
        print(f"  GWS hits     : {n_written:,}  -> {PROCESSED_GWAS.name}")
        print(f"  Skipped      : {n_skipped:,}")

    def _load_study_ancestry(self):
        """Build {pubmed_id: ancestry} from studies.tsv."""
        if not STUDIES_FILE.exists() or STUDIES_FILE.stat().st_size < 10_000:
            print("[GWAS] studies.tsv not available — ancestry will be 'NR'")
            return

        print(f"[GWAS] Loading study ancestry from {STUDIES_FILE.name}...")
        with open(STUDIES_FILE, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                pmid     = (row.get("PUBMEDID") or
                            row.get("PUBMED ID") or "").strip()
                ancestry = (row.get("BROAD ANCESTRAL CATEGORY") or
                            row.get("INITIAL SAMPLE SIZE") or "NR").strip()
                if pmid:
                    # Simplify ancestry labels
                    if "European" in ancestry:
                        ancestry = "European"
                    elif "East Asian" in ancestry:
                        ancestry = "East Asian"
                    elif "African" in ancestry:
                        ancestry = "African"
                    elif "South Asian" in ancestry:
                        ancestry = "South Asian"
                    elif "Hispanic" in ancestry or "Latin" in ancestry:
                        ancestry = "Hispanic/Latin American"
                    self._study_ancestry[pmid] = ancestry

        print(f"[GWAS] {len(self._study_ancestry):,} study ancestry entries loaded")

    def nodes(self) -> Generator[Node, None, None]:
        """
        GWAS variants are already in the graph from ClinVar (common ones)
        or get created here as new Variant nodes (novel to our graph).
        We yield new nodes only for rsIDs not already in the graph.
        The loader handles deduplication via INSERT OR IGNORE.
        """
        if not PROCESSED_GWAS.exists():
            raise FileNotFoundError("Run preprocess() first.")

        seen = set()
        with open(PROCESSED_GWAS, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rsid = row["rsid"]
                if rsid in seen:
                    continue
                seen.add(rsid)

                yield Node(
                    primary_id     = rsid,
                    primary_system = "dbSNP_rsID",
                    label          = rsid,
                    tier           = 1,
                    entity_type    = "Variant",
                    xrefs          = {},
                    properties     = {
                        "mapped_gene":  row["mapped_gene"],
                        "risk_allele":  row["risk_allele"],
                        "raf":          row["raf"],
                        "gwas_trait":   row["trait"],
                    },
                    source         = self.source_name,
                    source_version = self.source_version,
                    confidence     = float(row["confidence"]),
                )

    def edges(self) -> Generator[Edge, None, None]:
        """
        Yield INCREASES_RISK_OF edges: variant -> trait.
        Target is stored as trait label (text) — linked to disease
        nodes via label matching in the loader.
        All GWAS edges are INCREASES_RISK_OF, never CAUSES.
        """
        if not PROCESSED_GWAS.exists():
            raise FileNotFoundError("Run preprocess() first.")

        with open(PROCESSED_GWAS, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rsid      = row["rsid"]
                trait     = row["trait"]
                pubmed_id = row["pubmed_id"]
                or_beta   = row["or_beta"]
                mlog      = float(row["pvalue_mlog"])
                confidence= float(row["confidence"])
                ancestry  = row["ancestry"]

                # Parse effect size and direction
                effect_size = None
                direction   = None
                if or_beta:
                    try:
                        ef = float(or_beta)
                        effect_size = ef
                        direction = "positive" if ef >= 1.0 else "negative"
                    except ValueError:
                        pass

                try:
                    yield Edge(
                        source_id                = rsid,
                        source_system            = "dbSNP_rsID",
                        target_id                = trait,
                        target_system            = "GWAS_trait_label",
                        relationship_type        = "INCREASES_RISK_OF",
                        source_relationship_type = "GWAS_association",
                        effect_size              = effect_size,
                        effect_unit              = "OR_or_Beta",
                        direction                = direction,
                        confidence               = confidence,
                        primary_source           = f"PMID:{pubmed_id}",
                        imported_via             = f"GWAS_Catalog_{self.source_version}",
                        study_design             = "GWAS",
                        population_context       = ancestry,
                        source_version           = self.source_version,
                    )
                except ValueError:
                    continue

    def normalize_confidence(self, raw_value=None) -> float:
        if isinstance(raw_value, (int, float)):
            return _mlog_to_confidence(float(raw_value))
        return 0.55
