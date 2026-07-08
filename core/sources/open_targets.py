"""
core/sources/open_targets.py

Open Targets Platform 26.06 adapter for Disease-OS.

Files used:
  associations.parquet  74MB  — target-disease association scores
  disease.parquet        7MB  — disease metadata with synonyms/xrefs

What Open Targets adds:
  - 509,765 target→disease association edges scored across:
    genetics, somatic mutations, drugs, pathways, literature,
    animal models, RNA expression (aggregated into one score 0-1)
  - Disease nodes enriched with therapeutic area tags
  - Closes the Drug→Target→Disease chain when combined with
    DrugBank (which adds Drug→Target edges)

Edge type: Gene/Protein ASSOCIATED_WITH Disease
  - Uses Ensembl gene IDs on source side
  - Uses EFO IDs on target side
  - associationScore → confidence (filter >= 0.05 to keep meaningful hits)

ID mapping:
  targetId  ENSG00000009413 → stored as source_system=Ensembl
  diseaseId EFO_0000180     → normalised to EFO:EFO_0000180 to match our nodes
"""

import sys
import json
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS, PROCESSED_DIR

RAW_DIR      = Path.home() / "disease-os" / "data" / "raw" / "open_targets"
ASSOC_FILE   = RAW_DIR / "associations.parquet"
DISEASE_FILE = RAW_DIR / "disease.parquet"
TARGET_FILE  = RAW_DIR / "target.parquet"

# Minimum association score to load — below 0.05 is very weak
MIN_SCORE = 0.05


def _normalise_disease_id(eid: str) -> str:
    """
    Convert Open Targets EFO format to our graph format.
    EFO_0000180  ->  EFO:EFO_0000180
    MONDO_0005148 -> MONDO:MONDO_0005148
    HP_0001513   ->  HP:HP_0001513
    """
    if not eid or "_" not in eid:
        return eid
    prefix = eid.split("_")[0]
    return f"{prefix}:{eid}"


def _score_to_confidence(score: float) -> float:
    """
    Open Targets association score (0-1) → our confidence scale.
    OT scores are already well-calibrated 0-1 values.
    We cap at 0.85 because OT aggregates association evidence,
    not direct causation.
    """
    return min(round(score, 3), 0.85)


class OpenTargetsSource(BaseSource):
    """
    Loads Open Targets 26.06 associations into Disease-OS.

    Usage:
        source = OpenTargetsSource()
        source.load_into(graph_store)
    """

    source_name    = "OPEN_TARGETS"
    source_version = "26.06"

    def nodes(self) -> Generator[Node, None, None]:
        """
        Yield disease nodes enriched from disease.parquet.
        These complement our existing EFO nodes with OT-specific
        therapeutic area annotations.
        """
        import pyarrow.parquet as pq

        if not DISEASE_FILE.exists():
            print("[OT] disease.parquet not found — skipping disease nodes")
            return

        print(f"[OT] Loading disease nodes from {DISEASE_FILE.name}...")
        table = pq.read_table(DISEASE_FILE,
                              columns=["id", "name", "description",
                                       "dbXRefs", "exactSynonyms",
                                       "therapeuticAreas"])
        n = 0
        for batch in table.to_batches(max_chunksize=1000):
            for row in batch.to_pylist():
                eid  = row.get("id", "")
                name = row.get("name", "")

                if not eid or not name:
                    continue

                # Only load disease/phenotype nodes (not GO biological processes)
                prefix = eid.split("_")[0] if "_" in eid else eid.split(":")[0]
                if prefix not in ("EFO", "MONDO", "HP", "Orphanet", "DOID"):
                    continue

                primary_id = _normalise_disease_id(eid)

                # Extract xrefs from dbXRefs list
                xrefs = {}
                for xref in (row.get("dbXRefs") or []):
                    if not xref or ":" not in xref:
                        continue
                    sys_name, _, code = xref.partition(":")
                    if sys_name in ("ICD10", "ICD10CM", "ICD10WHO"):
                        xrefs["ICD-10"] = code
                    elif sys_name == "OMIM":
                        xrefs["OMIM"] = code
                    elif sys_name in ("MESH", "MeSH"):
                        xrefs["MeSH"] = code

                # Therapeutic areas as properties
                tas = row.get("therapeuticAreas") or []

                synonyms = row.get("exactSynonyms") or []

                yield Node(
                    primary_id     = primary_id,
                    primary_system = "EFO",
                    label          = name,
                    tier           = 7,
                    entity_type    = "Disease_scientific",
                    xrefs          = xrefs,
                    synonyms       = synonyms[:15],
                    definition     = (row.get("description") or "")[:500],
                    icd10_code     = xrefs.get("ICD-10"),
                    properties     = {
                        "therapeutic_areas": tas[:5],
                        "ot_id": eid,
                    },
                    source         = self.source_name,
                    source_version = self.source_version,
                    confidence     = 0.90,
                )
                n += 1

        print(f"[OT] {n:,} disease nodes yielded")

    def edges(self) -> Generator[Edge, None, None]:
        """
        Yield Gene ASSOCIATED_WITH Disease edges from associations.parquet.
        Filters to associationScore >= MIN_SCORE.
        """
        import pyarrow.parquet as pq

        if not ASSOC_FILE.exists():
            raise FileNotFoundError(f"Not found: {ASSOC_FILE}")

        print(f"[OT] Streaming {ASSOC_FILE.name} "
              f"({ASSOC_FILE.stat().st_size // 1_000_000}MB)...")
        print(f"[OT] Filter: associationScore >= {MIN_SCORE}")

        table = pq.read_table(ASSOC_FILE,
                              columns=["diseaseId", "targetId",
                                       "associationScore", "evidenceCount"])

        n_yielded = n_skipped = 0

        for batch in table.to_batches(max_chunksize=5000):
            for row in batch.to_pylist():
                disease_id = row.get("diseaseId", "")
                target_id  = row.get("targetId", "")
                score      = row.get("associationScore") or 0.0
                evidence   = row.get("evidenceCount") or 0

                if not disease_id or not target_id:
                    n_skipped += 1
                    continue

                if score < MIN_SCORE:
                    n_skipped += 1
                    continue

                # Normalise disease ID to our format
                norm_disease_id = _normalise_disease_id(disease_id)
                confidence      = _score_to_confidence(score)

                try:
                    yield Edge(
                        source_id                = target_id,
                        source_system            = "Ensembl",
                        target_id                = norm_disease_id,
                        target_system            = "EFO",
                        relationship_type        = "ASSOCIATED_WITH",
                        source_relationship_type = "ot_overall_association",
                        effect_size              = round(score, 4),
                        effect_unit              = "OT_association_score",
                        confidence               = confidence,
                        primary_source           = f"OpenTargets_{self.source_version}",
                        imported_via             = f"OpenTargets_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                    n_yielded += 1

                    if n_yielded % 50_000 == 0:
                        print(f"  {n_yielded:,} edges yielded, "
                              f"{n_skipped:,} skipped...")

                except ValueError:
                    n_skipped += 1
                    continue

        print(f"[OT] Done: {n_yielded:,} edges | {n_skipped:,} skipped")

    def normalize_confidence(self, raw_value=None) -> float:
        if isinstance(raw_value, (int, float)):
            return _score_to_confidence(float(raw_value))
        return 0.50
