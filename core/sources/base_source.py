"""
base_source.py — Abstract contract for every source adapter.
Enforces: canonical ID, confidence normalisation,
relationship normalisation, provenance.
"""
from abc import ABC, abstractmethod
from typing import Generator
from core.node import Node
from core.edge import Edge
from core.edge import CANONICAL_RELATIONSHIP_TYPES


class BaseSource(ABC):
    source_name:    str = ""
    source_version: str = ""

    @abstractmethod
    def nodes(self) -> Generator[Node, None, None]: ...

    @abstractmethod
    def edges(self) -> Generator[Edge, None, None]: ...

    @abstractmethod
    def normalize_confidence(self, raw_value) -> float: ...

    def normalize_relationship(self, raw_type: str) -> str:
        canonical = self._relationship_map().get(raw_type.upper(), "ASSOCIATED_WITH")
        return canonical if canonical in CANONICAL_RELATIONSHIP_TYPES else "ASSOCIATED_WITH"

    def _relationship_map(self) -> dict:
        return {}

    def provenance(self) -> dict:
        return {
            "source":         self.source_name,
            "source_version": self.source_version,
            "imported_via":   f"{self.source_name}_{self.source_version}",
        }

    def load_into(self, graph_store) -> dict:
        nodes_added = edges_added = edges_skipped = 0
        print(f"[{self.source_name}] Loading nodes...")
        for node in self.nodes():
            graph_store.add_node(node)
            nodes_added += 1
            if nodes_added % 10_000 == 0:
                print(f"  {nodes_added:,} nodes...")
        print(f"[{self.source_name}] Loading edges...")
        for edge in self.edges():
            if graph_store.add_edge(edge):
                edges_added += 1
            else:
                edges_skipped += 1
            total = edges_added + edges_skipped
            if total % 10_000 == 0:
                print(f"  {edges_added:,} edges added, {edges_skipped:,} skipped...")
        summary = {
            "source":        self.source_name,
            "nodes_added":   nodes_added,
            "edges_added":   edges_added,
            "edges_skipped": edges_skipped,
        }
        print(f"[{self.source_name}] Done: {summary}")
        return summary
