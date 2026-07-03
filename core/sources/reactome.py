"""
core/sources/reactome.py

Reactome adapter for Disease-OS.
Loads human pathways, protein-pathway membership, and pathway hierarchy.
Protein-protein interactions loaded separately once the file is available.

Files used:
  pathway_list.txt          -> Pathway nodes (Tier 2)
  pathway_hierarchy.txt     -> Pathway ISA/PART_OF edges
  pathways_to_proteins.txt  -> Protein PART_OF Pathway edges

Evidence codes in pathways_to_proteins.txt:
  TAS = Traceable Author Statement  (manually curated, high confidence)
  IEA = Inferred by Electronic Annotation (computational, lower confidence)
  NAS = Non-traceable Author Statement
  IC  = Inferred by Curator

All Reactome edges carry confidence >= 0.85 (TAS) or 0.70 (IEA).
This is significantly higher than UMLS ASSOCIATED_WITH edges (0.60).
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

# Evidence code -> confidence score
EVIDENCE_CONFIDENCE = {
    "TAS": 0.85,   # manually curated — traceable author statement
    "NAS": 0.75,   # non-traceable author statement
    "IEA": 0.70,   # inferred by electronic annotation
    "IC":  0.80,   # inferred by curator
}


class ReactomeSource(BaseSource):
    """
    Loads Reactome human pathways into Disease-OS.

    Usage:
        source = ReactomeSource()
        source.load_into(graph_store)
    """

    source_name    = "REACTOME"
    source_version = SOURCE_VERSIONS["REACTOME"]

    def nodes(self) -> Generator[Node, None, None]:
        """Yield Pathway nodes from pathway_list.txt."""
        yield from self._pathway_nodes()

    def edges(self) -> Generator[Edge, None, None]:
        """Yield all Reactome edges — hierarchy + protein membership."""
        yield from self._hierarchy_edges()
        yield from self._protein_pathway_edges()
        if INTERACTIONS.exists() and not self._is_html(INTERACTIONS):
            yield from self._interaction_edges()
        else:
            print("[REACTOME] interactions.txt not available — skipping PPI edges")
            print("           Download from reactome.org and re-run to add them")

    def normalize_confidence(self, raw_value=None) -> float:
        if isinstance(raw_value, str):
            return EVIDENCE_CONFIDENCE.get(raw_value, 0.70)
        return 0.85

    # ── Node generators ────────────────────────────────────────────────

    def _pathway_nodes(self) -> Generator[Node, None, None]:
        """
        pathway_list.txt columns (tab-delimited, no header):
          0: Pathway stable ID  (e.g. R-HSA-164843)
          1: Pathway name       (e.g. 2-LTR circle formation)
          2: Species            (e.g. Homo sapiens)
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

                # Human only
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
        pathway_hierarchy.txt columns (tab-delimited, no header):
          0: Parent pathway ID
          1: Child pathway ID

        Relationship: Child ISA Parent (child is a sub-pathway of parent)
        """
        n = 0
        with open(PATHWAY_HIER, encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                parent = parts[0].strip()
                child  = parts[1].strip()

                # Human only
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
                        imported_via             = f"Reactome_pathway_hierarchy_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                    n += 1
                except ValueError:
                    continue
        print(f"[REACTOME] {n:,} pathway hierarchy edges")

    def _protein_pathway_edges(self) -> Generator[Edge, None, None]:
        """
        pathways_to_proteins.txt columns (tab-delimited, no header):
          0: UniProt accession
          1: Pathway stable ID
          2: URL
          3: Pathway name
          4: Evidence code (TAS/IEA/NAS/IC)
          5: Species

        Relationship: Protein PART_OF Pathway
        """
        n     = 0
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

                confidence = self.normalize_confidence(evidence)

                try:
                    yield Edge(
                        source_id                = uniprot,
                        source_system            = "UniProt",
                        target_id                = pathway_id,
                        target_system            = "Reactome",
                        relationship_type        = "PART_OF",
                        source_relationship_type = f"pathway_member_{evidence}",
                        confidence               = confidence,
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
        interactions.txt — protein-protein interactions with annotation.
        Loaded only if file exists and is not an HTML error page.
        Column structure varies by Reactome release — detected dynamically.
        """
        print("[REACTOME] Loading protein-protein interactions...")
        n = 0
        with open(INTERACTIONS, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # Try to extract the two proteins and interaction type
                # Column names vary — handle common variants
                prot_a = (row.get("# Interactor 1") or
                          row.get("Interactor 1") or
                          row.get("Gene1") or "").strip()
                prot_b = (row.get("Interactor 2") or
                          row.get("Gene2") or "").strip()
                itype  = (row.get("Interaction Type") or
                          row.get("Annotation") or "").strip().upper()

                if not prot_a or not prot_b:
                    continue

                # Map interaction type to canonical relationship
                rel_type = {
                    "ACTIVATION":    "ACTIVATES",
                    "INHIBITION":    "INHIBITS",
                    "BINDING":       "BINDS",
                    "CATALYSIS":     "CATALYZES",
                    "EXPRESSION":    "UPREGULATES",
                    "INPUT":         "ASSOCIATED_WITH",
                    "REACTION":      "ASSOCIATED_WITH",
                }.get(itype, "INTERACTS_WITH")

                try:
                    yield Edge(
                        source_id                = prot_a,
                        source_system            = "UniProt",
                        target_id                = prot_b,
                        target_system            = "UniProt",
                        relationship_type        = rel_type,
                        source_relationship_type = itype,
                        confidence               = 0.85,
                        primary_source           = f"Reactome_{self.source_version}",
                        imported_via             = f"Reactome_interactions_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                    n += 1
                except ValueError:
                    continue

        print(f"[REACTOME] {n:,} protein-protein interaction edges")

    def _is_html(self, path: Path) -> bool:
        """Check if a downloaded file is actually an HTML error page."""
        try:
            with open(path, encoding="utf-8") as f:
                first_line = f.readline().strip()
                return first_line.startswith("<!DOCTYPE") or first_line.startswith("<html")
        except Exception:
            return True
