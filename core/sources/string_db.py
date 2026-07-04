"""
core/sources/string_db.py

STRING v12.0 adapter for Disease-OS.

Files used (human only — taxon 9606):
  protein_info.txt   6MB   — Ensembl protein ID -> gene symbol + annotation
  protein_links.txt  602MB — protein1 | protein2 | combined_score (0-1000)

STRING combined_score -> confidence mapping:
  >= 900  -> 0.90  (highest confidence)
  >= 700  -> 0.75  (high confidence)
  >= 400  -> 0.60  (medium confidence)
  <  400  -> skip  (low confidence — too noisy for our graph)

We only load edges with combined_score >= 400 (STRING's "medium" threshold).
This cuts 602MB down to ~30% of rows.

Edge type:
  All STRING edges -> INTERACTS_WITH
  STRING does not specify direction or mechanism type in the links file.
  Directional/mechanistic types (ACTIVATES, INHIBITS) come from
  the experimental/pathway-derived files which we add in Phase 3.

Node strategy:
  STRING proteins are already in our graph as UniProt nodes (from UMLS).
  We do NOT create new nodes — we only add edges between existing nodes.
  Linking: STRING ENSP ID -> gene symbol -> match against node properties.
"""

import sys
import csv
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.edge import Edge
from core.node import Node
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS

RAW_DIR      = Path.home() / "disease-os" / "data" / "raw" / "string"
PROTEIN_INFO = RAW_DIR / "protein_info.txt"
PROTEIN_LINKS= RAW_DIR / "protein_links.txt"

MIN_SCORE    = 400    # STRING medium-confidence threshold
TAXON_PREFIX = "9606."

def _score_to_confidence(score: int) -> float:
    if score >= 900: return 0.90
    if score >= 700: return 0.75
    if score >= 400: return 0.60
    return 0.0


class StringSource(BaseSource):
    """
    Loads STRING v12.0 protein-protein interactions into Disease-OS.

    Strategy:
      1. Build ENSP -> gene_symbol index from protein_info.txt
      2. Stream protein_links.txt, filter by score >= MIN_SCORE
      3. Yield INTERACTS_WITH edges between ENSP IDs
         (stored as HGNC_Symbol system so they join to existing nodes)

    Usage:
        source = StringSource()
        source.load_into(graph_store)
    """

    source_name    = "STRING"
    source_version = SOURCE_VERSIONS["STRING"]

    def __init__(self):
        self._ensp_to_gene: dict[str, str] = {}

    def _load_protein_index(self):
        """
        Build {ENSP_id: gene_symbol} from protein_info.txt.
        Called once before streaming edges.
        """
        if self._ensp_to_gene:
            return   # already loaded

        print(f"[STRING] Loading protein index from {PROTEIN_INFO.name}...")
        with open(PROTEIN_INFO, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                raw_id   = row.get("#string_protein_id", "").strip()
                gene_sym = row.get("preferred_name", "").strip()
                if raw_id and gene_sym:
                    # Strip taxon prefix: "9606.ENSP00000000233" -> "ENSP00000000233"
                    ensp = raw_id.replace(TAXON_PREFIX, "")
                    self._ensp_to_gene[ensp] = gene_sym

        print(f"[STRING] {len(self._ensp_to_gene):,} proteins indexed")

    def nodes(self) -> Generator[Node, None, None]:
        """
        STRING does not add new node types — proteins are already in the
        graph from UMLS/UniProt. Yield nothing here.
        """
        return
        yield   # makes this a generator

    def edges(self) -> Generator[Edge, None, None]:
        """
        Stream protein_links.txt and yield INTERACTS_WITH edges.
        Filters: combined_score >= MIN_SCORE (400).
        """
        self._load_protein_index()

        print(f"[STRING] Streaming {PROTEIN_LINKS.name} "
              f"(score >= {MIN_SCORE})...")
        n_yielded = n_skipped_score = n_skipped_lookup = 0

        with open(PROTEIN_LINKS, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=" ")
            for row in reader:
                raw1  = row.get("protein1", "").strip()
                raw2  = row.get("protein2", "").strip()
                score_str = row.get("combined_score", "0").strip()

                try:
                    score = int(score_str)
                except ValueError:
                    n_skipped_score += 1
                    continue

                # Filter by confidence threshold
                if score < MIN_SCORE:
                    n_skipped_score += 1
                    continue

                # Strip taxon prefix
                ensp1 = raw1.replace(TAXON_PREFIX, "")
                ensp2 = raw2.replace(TAXON_PREFIX, "")

                # Look up gene symbols
                gene1 = self._ensp_to_gene.get(ensp1)
                gene2 = self._ensp_to_gene.get(ensp2)

                if not gene1 or not gene2:
                    n_skipped_lookup += 1
                    continue

                confidence = _score_to_confidence(score)

                try:
                    yield Edge(
                        source_id                = gene1,
                        source_system            = "HGNC_Symbol",
                        target_id                = gene2,
                        target_system            = "HGNC_Symbol",
                        relationship_type        = "INTERACTS_WITH",
                        source_relationship_type = f"STRING_combined_{score}",
                        effect_size              = score / 1000.0,
                        effect_unit              = "STRING_score",
                        confidence               = confidence,
                        primary_source           = f"STRING_{self.source_version}",
                        imported_via             = f"STRING_protein_links_{self.source_version}",
                        study_design             = "curated",
                        population_context       = "human",
                        source_version           = self.source_version,
                    )
                    n_yielded += 1
                    if n_yielded % 500_000 == 0:
                        print(f"  {n_yielded:,} edges yielded, "
                              f"{n_skipped_score:,} low-score skipped, "
                              f"{n_skipped_lookup:,} lookup misses...")
                except ValueError:
                    n_skipped_lookup += 1
                    continue

        print(f"[STRING] Done: {n_yielded:,} edges | "
              f"{n_skipped_score:,} low-score | "
              f"{n_skipped_lookup:,} lookup miss")

    def normalize_confidence(self, raw_value=None) -> float:
        if isinstance(raw_value, (int, float)):
            return _score_to_confidence(int(raw_value))
        return 0.60

    def _relationship_map(self) -> dict:
        return {}   # STRING doesn't provide typed relationships in links file
