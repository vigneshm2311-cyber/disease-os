"""
core/sources/reactome.py

Reactome adapter for Disease-OS.
Loads human pathways, protein-pathway membership,
pathway hierarchy, and functional gene interactions.

Files used:
  pathway_list.txt          -> Pathway nodes (Tier 2)
  pathway_hierarchy.txt     -> Pathway ISA/PART_OF edges
  pathways_to_proteins.txt  -> Protein PART_OF Pathway edges
  interactions.txt          -> Gene/protein functional interactions

Evidence codes in pathways_to_proteins.txt:
  TAS = Traceable Author Statement  (manually curated, confidence 0.85)
  IEA = Inferred by Electronic Annotation (computational,  confidence 0.70)
  NAS = Non-traceable Author Statement                    (confidence 0.75)
  IC  = Inferred by Curator                              (confidence 0.80)

interactions.txt annotation types -> canonical edge types:
  activation, expression    -> ACTIVATES / UPREGULATES
  inhibition                -> INHIBITS
  catalyzed by              -> CATALYZES
  complex, binding          -> BINDS
  predicted, input, others  -> ASSOCIATED_WITH (lower confidence)
"""

import sys
import csv
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS

RAW_DIR          = Path.home() / "disease-os" / "data" / "raw" / "reactome"
PATHWAY_LIST     = RAW_DIR / "pathway_list.txt"
PATHWAY_HIER     = RAW_DIR / "pathway_hierarchy.txt"
PATHWAYS_TO_PROT = RAW_DIR / "pathways_to_proteins.txt"
INTERACTIONS     = RAW_DIR / "interactions.txt"

EVIDENCE_CONFIDENCE = {
    "TAS": 0.85,
    "NAS": 0.75,
    "IEA": 0.70,
    "IC":  0.80,
}

# Annotation keyword -> canonical relationship type
# interactions.txt Annotation column can contain semicolon-separated terms
ANNOTATION_REL_MAP = {
    "activation":    "ACTIVATES",
    "expression":    "UPREGULATES",
    "inhibition":    "INHIBITS",
    "catalyzed by":  "CATALYZES",
    "catalysis":     "CATALYZES",
    "complex":       "BINDS",
    "binding":       "BINDS",
    "input":         "ASSOCIATED_WITH",
    "output":        "ASSOCIATED_WITH",
    "reaction":      "ASSOCIATED_WITH",
    "predicted":     "ASSOCIATED_WITH",
}


def _annotation_to_rel(annotation: str) -> tuple[str, float]:
    """
    Parse a semicolon-separated Reactome annotation string.
    Returns (canonical_relationship_type, confidence).
    Picks the most specific / highest-confidence type if multiple.
    """
    if not annotation:
        return "ASSOCIATED_WITH", 0.60

    terms = [t.strip().lower() for t in annotation.split(";")]

    # Priority order — more specific types win
    priority = [
        "catalyzed by", "catalysis",
        "activation", "inhibition", "expression",
        "complex", "binding",
        "input", "output", "reaction", "predicted",
    ]
    for p in priority:
        if any(p in t for t in terms):
            rel = ANNOTATION_REL_MAP.get(p, "ASSOCIATED_WITH")
            # Curated mechanistic types get higher confidence than predicted
            conf = 0.85 if p not in ("predicted", "input", "output") else 0.65
            return rel, conf

    return "ASSOCIATED_WITH", 0.60


class ReactomeSource(BaseSource):
    """
    Loads Reactome human pathways and interactions into Disease-OS.

    Usage:
        source = ReactomeSource()
        source.load_into(graph_store)
    """

    source_name    = "REACTOME"
    source_version = SOURCE_VERSIONS["REACTOME"]

    def nodes(self) -> Generator[Node, None, None]:
        yield from self._pathway_nodes()

    def edges(self) -> Generator[Edge, None, None]:
        yield from self._hierarchy_edges()
        yield from self._protein_pathway_edges()
        if INTERACTIONS.exists() and not self._is_html(INTERACTIONS):
            yield from self._interaction_edges()
        else:
            print("[REACTOME] interactions.txt missing or invalid — skipping PPI")

    def normalize_confidence(self, raw_value=None) -> float:
        if isinstance(raw_value, str):
            return EVIDENCE_CONFIDENCE.get(raw_value, 0.70)
        return 0.85

    # ── Node generators ────────────────────────────────────────────────

    def _pathway_nodes(self) -> Generator[Node, None, None]:
        """
        pathway_list.txt (no header, tab-delimited):
          col 0: Pathway stable ID  e.g. R-HSA-164843
          col 1: Pathway name       e.g. 2-LTR circle formation
          col 2: Species            e.g. Homo sapiens
        """
        n = 0
        with open(PATHWAY_LIST, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                pathway_id = parts[0].strip()
                name       = parts[1].strip()
                species    = parts[2].strip()

                if species != "Homo sapiens":
                    continue
                if not pathway_id.startswith("R-HSA-"):
                    continue

                yield Node(
                    primary_id     = pathway_id,
                    primary_system = "Reactome",
                    label          = name,
                    tier           = 2,
                    entity_type    = "Pathway",
                    xrefs          = {
                        "Reactome": pathway_id,
                        "URL": f"https://reactome.org/PathwayBrowser/#/{pathway_id}",
                    },
                    properties     = {"species": species},
                    source         = self.source_name,
                    source_version = self.source_version,
                    confidence     = 1.0,
                )
                n += 1
        print(f"[REACTOME] {n:,} human pathway nodes")

    # ── Edge generators ────────────────────────────────────────────────

    def _hierarchy_edges(self) -> Generator[Edge, None, None]:
        """
        pathway_hierarchy.txt (no header, tab-delimited):
          col 0: Parent pathway ID
          col 1: Child pathway ID

        Semantics: child is a sub-pathway of parent -> child PART_OF parent
        """
        n = 0
        with open(PATHWAY_HIER, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                parent = parts[0].strip()
                child  = parts[1].strip()

                if not parent.startswith("R-HSA-"):
                    continue
                if not child.startswith("R-HSA-"):
                    continue

                try:
                    yield Edge(
                        source_id                = child,
                        source_system            = "Reactome",
                        target_id                = parent,
                        target_system            = "Reactome",
                        relationship_type        = "PART_OF",
                        source_relationship_type = "has_parent_pathway",
                        confidence               = 1.0,
                        primary_source           = f"Reactome_{self.source_version}",
                        imported_via             = f"Reactome_hierarchy_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                    n += 1
                except ValueError:
                    continue
        print(f"[REACTOME] {n:,} pathway hierarchy edges")

    def _protein_pathway_edges(self) -> Generator[Edge, None, None]:
        """
        pathways_to_proteins.txt (no header, tab-delimited):
          col 0: UniProt accession
          col 1: Pathway stable ID
          col 2: URL
          col 3: Pathway name
          col 4: Evidence code (TAS/IEA/NAS/IC)
          col 5: Species

        Semantics: protein PART_OF pathway
        """
        n       = 0
        skipped = 0
        with open(PATHWAYS_TO_PROT, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                uniprot    = parts[0].strip()
                pathway_id = parts[1].strip()
                evidence   = parts[4].strip()
                species    = parts[5].strip()

                if species != "Homo sapiens":
                    skipped += 1
                    continue
                if not pathway_id.startswith("R-HSA-"):
                    skipped += 1
                    continue
                if not uniprot:
                    skipped += 1
                    continue

                try:
                    yield Edge(
                        source_id                = uniprot,
                        source_system            = "UniProt",
                        target_id                = pathway_id,
                        target_system            = "Reactome",
                        relationship_type        = "PART_OF",
                        source_relationship_type = f"pathway_member_{evidence}",
                        confidence               = self.normalize_confidence(evidence),
                        primary_source           = f"Reactome_{self.source_version}",
                        imported_via             = f"Reactome_UniProt2Reactome_{self.source_version}",
                        study_design             = "curated",
                        population_context       = "human",
                        source_version           = self.source_version,
                    )
                    n += 1
                except ValueError:
                    skipped += 1
                    continue

        print(f"[REACTOME] {n:,} protein->pathway edges ({skipped:,} non-human skipped)")

    def _interaction_edges(self) -> Generator[Edge, None, None]:
        """
        interactions.txt (tab-delimited with header):
          Gene1       e.g. A1CF
          Gene2       e.g. APOBEC1
          Annotation  e.g. catalyzed by; complex; input
          Direction   e.g. <- / -> / -
          Score       e.g. 1.00

        Direction semantics:
          ->   Gene1 acts on Gene2
          <-   Gene2 acts on Gene1 (swap source/target)
          -    undirected

        Gene symbols are stored as HGNC_Symbol.
        These edges will be resolved to canonical NCBI Gene IDs
        when UniProt/Ensembl data is loaded in a later step.
        """
        n       = 0
        skipped = 0

        with open(INTERACTIONS, encoding="utf-8") as f:
            # Handle header line — Reactome uses "Gene1" directly
            reader = csv.DictReader(f, delimiter="\t")

            # Normalise header — strip leading # if present
            reader.fieldnames = [
                h.lstrip("#").strip() for h in (reader.fieldnames or [])
            ]

            for row in reader:
                gene1      = row.get("Gene1", "").strip()
                gene2      = row.get("Gene2", "").strip()
                annotation = row.get("Annotation", "").strip()
                direction  = row.get("Direction", "-").strip()
                score_str  = row.get("Score", "0").strip()

                if not gene1 or not gene2:
                    skipped += 1
                    continue

                # Parse score
                try:
                    score = float(score_str)
                except ValueError:
                    score = 0.70

                # Determine canonical relationship + base confidence
                rel_type, base_conf = _annotation_to_rel(annotation)

                # Final confidence = average of annotation confidence + score
                confidence = round((base_conf + min(score, 1.0)) / 2, 3)

                # Apply direction — swap source/target if arrow points left
                if direction == "<-":
                    src, tgt = gene2, gene1
                else:
                    src, tgt = gene1, gene2

                try:
                    yield Edge(
                        source_id                = src,
                        source_system            = "HGNC_Symbol",
                        target_id                = tgt,
                        target_system            = "HGNC_Symbol",
                        relationship_type        = rel_type,
                        source_relationship_type = annotation,
                        effect_size              = score,
                        effect_unit              = "ReactomeScore",
                        confidence               = confidence,
                        primary_source           = f"Reactome_{self.source_version}",
                        imported_via             = f"Reactome_FI_{self.source_version}",
                        study_design             = "curated",
                        population_context       = "human",
                        source_version           = self.source_version,
                    )
                    n += 1
                except ValueError:
                    skipped += 1
                    continue

        print(f"[REACTOME] {n:,} functional interaction edges ({skipped:,} skipped)")

    def _is_html(self, path: Path) -> bool:
        """True if the file is an HTML error page, not real data."""
        try:
            with open(path, encoding="utf-8") as f:
                first = f.readline().strip()
                return first.startswith("<!DOCTYPE") or first.startswith("<html")
        except Exception:
            return True
