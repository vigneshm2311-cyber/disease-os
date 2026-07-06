"""
core/causal_engine.py

Causal Reasoning Engine for Disease-OS.

Implements backward causal tracing — starting from a disease node,
traversing upstream along causal edges to identify contributing
mechanisms, ranked by path confidence.

Algorithm:
  Backward BFS from anchor node, following causal edges in reverse.
  Each path is scored as the product of edge confidences, with a
  depth discount to prefer short strong chains over long weak ones.

Usage:
    from core.causal_engine import CausalEngine
    engine = CausalEngine(db_path)

    # Trace root causes of T2D
    results = engine.trace(
        disease_id     = "E11.9",
        disease_system = "ICD-10-CM",
        max_depth      = 4,
        min_confidence = 0.50,
        top_n          = 20,
    )

    for path in results:
        print(path.summary())
"""

import sqlite3
import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


# ── Edge types followed upstream, in priority order ───────────────────────────
CAUSAL_EDGE_TYPES = [
    "CAUSES",
    "CONTRIBUTES_TO",
    "INCREASES_RISK_OF",
    "UPREGULATES",
    "DOWNREGULATES",
    "ACTIVATES",
    "INHIBITS",
    "ASSOCIATED_WITH",   # weakest — only if nothing stronger
]

# Edges that can be followed in either direction (undirected)
UNDIRECTED_TYPES = {"INTERACTS_WITH", "BINDS", "ASSOCIATED_WITH"}

# Depth discount factor — each hop multiplies path score by this
DEPTH_DISCOUNT = 0.85

# Maximum nodes to expand per BFS level (prevents explosion)
MAX_EXPAND_PER_LEVEL = 500


@dataclass
class CausalEdge:
    """One edge in a causal path."""
    source_id:     str
    source_system: str
    source_label:  str
    rel_type:      str
    target_id:     str
    target_system: str
    target_label:  str
    confidence:    float
    effect_size:   Optional[float]
    effect_unit:   Optional[str]
    study_design:  str
    primary_source: str
    tier_from:     int
    tier_to:       int

    def arrow(self) -> str:
        return f"--[{self.rel_type} conf={self.confidence:.2f}]-->"


@dataclass
class CausalPath:
    """A complete causal chain from a root cause to the anchor disease."""
    edges:      list[CausalEdge]
    path_score: float
    depth:      int

    def root(self) -> CausalEdge:
        return self.edges[0]

    def anchor(self) -> CausalEdge:
        return self.edges[-1]

    def nodes(self) -> list[tuple]:
        """Return list of (id, system, label, tier) along the path."""
        result = []
        if self.edges:
            e = self.edges[0]
            result.append((e.source_id, e.source_system,
                           e.source_label, e.tier_from))
        for e in self.edges:
            result.append((e.target_id, e.target_system,
                           e.target_label, e.tier_to))
        return result

    def intervention_points(self) -> list[tuple]:
        """
        Nodes in this path that could be intervened on.
        Returns list of (id, system, label, tier).
        Excludes the anchor disease itself.
        """
        return [n for n in self.nodes()[:-1]]

    def summary(self) -> str:
        lines = [f"Path score: {self.path_score:.4f}  depth: {self.depth}"]
        for e in self.edges:
            lines.append(
                f"  [{e.tier_from}] {e.source_label[:40]:<40} "
                f"{e.arrow()} "
                f"[{e.tier_to}] {e.target_label[:40]}"
            )
            if e.effect_size:
                lines.append(
                    f"       effect={e.effect_size} {e.effect_unit or ''} "
                    f"study={e.study_design} src={e.primary_source[:40]}"
                )
        return "\n".join(lines)

    def to_dict(self) -> dict:
        return {
            "path_score": round(self.path_score, 4),
            "depth":      self.depth,
            "edges": [
                {
                    "source":        e.source_label,
                    "source_id":     e.source_id,
                    "source_tier":   e.tier_from,
                    "relationship":  e.rel_type,
                    "confidence":    e.confidence,
                    "effect_size":   e.effect_size,
                    "target":        e.target_label,
                    "target_id":     e.target_id,
                    "target_tier":   e.tier_to,
                    "study_design":  e.study_design,
                    "source_ref":    e.primary_source,
                }
                for e in self.edges
            ],
            "root_cause":       self.edges[0].source_label if self.edges else "",
            "root_tier":        self.edges[0].tier_from    if self.edges else 0,
            "intervention_points": [
                {"id": n[0], "system": n[1], "label": n[2], "tier": n[3]}
                for n in self.intervention_points()
            ],
        }


class CausalEngine:
    """
    Backward causal reasoning engine over the Disease-OS graph.
    Traces from a disease node upstream to identify root causes.
    """

    def __init__(self, db_path: str):
        self.db_path = str(db_path)
        self.conn    = sqlite3.connect(self.db_path)
        self.conn.execute("PRAGMA cache_size=-256000;")
        self.conn.execute("PRAGMA temp_store=MEMORY;")
        # Label cache — avoid repeated lookups
        self._label_cache: dict[tuple, tuple] = {}

    # ── Public API ─────────────────────────────────────────────────────────────

    def trace(
        self,
        disease_id:     str,
        disease_system: str,
        max_depth:      int   = 4,
        min_confidence: float = 0.50,
        min_path_score: float = 0.05,
        top_n:          int   = 30,
        patient_data:   dict  = None,
        edge_types:     list  = None,
    ) -> list[CausalPath]:
        """
        Trace root causes of a disease using backward BFS.

        Args:
            disease_id:     Primary ID of the anchor disease node
            disease_system: Primary system of the anchor disease node
            max_depth:      Maximum hops upstream (default 4)
            min_confidence: Minimum edge confidence to follow (default 0.50)
            min_path_score: Minimum overall path score to return (default 0.05)
            top_n:          Return top N paths by path_score (default 30)
            patient_data:   Optional dict of patient-specific data for
                            personalized path scoring
            edge_types:     Optional list of edge types to follow
                            (defaults to CAUSAL_EDGE_TYPES)

        Returns:
            List of CausalPath objects, ranked by path_score descending
        """
        if edge_types is None:
            edge_types = CAUSAL_EDGE_TYPES

        # Get anchor node
        anchor = self._get_node(disease_id, disease_system)
        if not anchor:
            raise ValueError(
                f"Node not found: {disease_id} ({disease_system})"
            )

        anchor_id, anchor_sys, anchor_label, anchor_tier = anchor

        # BFS state:
        # frontier = list of (current_node_id, current_node_system,
        #                     path_so_far, path_score, visited_set)
        frontier = [(
            anchor_id,
            anchor_sys,
            [],     # path edges so far (empty at start)
            1.0,    # path score starts at 1.0
            {(anchor_id, anchor_sys)},  # visited nodes
        )]

        all_paths: list[CausalPath] = []

        for depth in range(1, max_depth + 1):
            next_frontier = []

            # Limit expansion per level
            if len(frontier) > MAX_EXPAND_PER_LEVEL:
                # Keep highest-scoring paths
                frontier.sort(key=lambda x: x[3], reverse=True)
                frontier = frontier[:MAX_EXPAND_PER_LEVEL]

            for (node_id, node_sys, path_edges,
                 path_score, visited) in frontier:

                # Get all upstream edges for this node
                upstream = self._get_upstream_edges(
                    node_id, node_sys,
                    edge_types, min_confidence
                )

                for edge in upstream:
                    src_key = (edge.source_id, edge.source_system)

                    # Skip visited nodes (prevent cycles)
                    if src_key in visited:
                        continue

                    # Score this path
                    edge_score     = edge.confidence
                    new_path_score = (path_score * edge_score
                                      * (DEPTH_DISCOUNT ** depth))

                    # Apply patient personalization if provided
                    if patient_data:
                        new_path_score = self._personalize_score(
                            new_path_score, edge, patient_data
                        )

                    if new_path_score < min_path_score:
                        continue

                    new_path  = path_edges + [edge]
                    new_visited = visited | {src_key}

                    # Record this as a complete path
                    all_paths.append(CausalPath(
                        edges      = new_path,
                        path_score = new_path_score,
                        depth      = depth,
                    ))

                    # Add to next frontier for deeper exploration
                    next_frontier.append((
                        edge.source_id,
                        edge.source_system,
                        new_path,
                        new_path_score,
                        new_visited,
                    ))

            frontier = next_frontier
            if not frontier:
                break

        # Rank by path score, deduplicate by root cause node
        all_paths.sort(key=lambda p: p.path_score, reverse=True)

        # Return top N
        return all_paths[:top_n]

    def trace_between(
        self,
        source_id:     str,
        source_system: str,
        target_id:     str,
        target_system: str,
        max_depth:     int = 5,
    ) -> list[CausalPath]:
        """
        Find causal paths between two specific nodes.
        Useful for: "How does poor sleep cause insulin resistance?"
        """
        all_paths = self.trace(
            disease_id     = target_id,
            disease_system = target_system,
            max_depth      = max_depth,
            min_confidence = 0.40,
            min_path_score = 0.01,
            top_n          = 100,
        )

        # Filter to paths that contain the source node
        filtered = [
            p for p in all_paths
            if any(
                e.source_id == source_id
                and e.source_system == source_system
                for e in p.edges
            )
        ]
        return filtered[:10]

    def intervention_impact(
        self,
        disease_id:     str,
        disease_system: str,
        intervene_on:   str,   # node primary_id to intervene on
        max_depth:      int = 4,
    ) -> dict:
        """
        Estimate the impact of removing a node from the causal graph.
        Returns: paths that pass through that node and their combined score.
        Used for: "If we fix sleep, how much of T2D risk is addressable?"
        """
        all_paths = self.trace(
            disease_id     = disease_id,
            disease_system = disease_system,
            max_depth      = max_depth,
            top_n          = 200,
        )

        affected = [
            p for p in all_paths
            if any(
                e.source_id == intervene_on or e.target_id == intervene_on
                for e in p.edges
            )
        ]

        total_score    = sum(p.path_score for p in all_paths)
        affected_score = sum(p.path_score for p in affected)

        return {
            "intervene_on":     intervene_on,
            "total_paths":      len(all_paths),
            "affected_paths":   len(affected),
            "total_score":      round(total_score, 4),
            "affected_score":   round(affected_score, 4),
            "fraction_covered": round(
                affected_score / total_score if total_score > 0 else 0, 3
            ),
            "top_affected_paths": [p.to_dict() for p in affected[:5]],
        }

    # ── Internal helpers ───────────────────────────────────────────────────────

    def _get_node(
        self, node_id: str, node_system: str
    ) -> tuple | None:
        """Get (id, system, label, tier) for a node."""
        key = (node_id, node_system)
        if key in self._label_cache:
            return self._label_cache[key]

        row = self.conn.execute("""
            SELECT primary_id, primary_system, label, tier
            FROM nodes
            WHERE primary_id = ? AND primary_system = ?
        """, (node_id, node_system)).fetchone()

        if row:
            self._label_cache[key] = row
        return row

    def _get_node_by_id_only(self, node_id: str) -> tuple | None:
        """Get node when system is unknown — tries primary_id match."""
        row = self.conn.execute("""
            SELECT primary_id, primary_system, label, tier
            FROM nodes
            WHERE primary_id = ?
            LIMIT 1
        """, (node_id,)).fetchone()
        return row

    def _get_upstream_edges(
        self,
        node_id:        str,
        node_system:    str,
        edge_types:     list,
        min_confidence: float,
    ) -> list[CausalEdge]:
        """
        Get all edges pointing TO this node (upstream causes).
        In the graph: source --[CAUSES]--> target
        We want edges where target = this node,
        then we traverse to source (upstream).
        """
        placeholders = ",".join(["?"] * len(edge_types))

        rows = self.conn.execute(f"""
            SELECT
                e.source_id, e.source_system,
                e.target_id, e.target_system,
                e.relationship_type,
                e.confidence,
                e.effect_size,
                e.effect_unit,
                e.study_design,
                e.primary_source
            FROM edges e
            WHERE e.target_id = ?
              AND e.target_system = ?
              AND e.relationship_type IN ({placeholders})
              AND e.confidence >= ?
            ORDER BY e.confidence DESC
            LIMIT 100
        """, [node_id, node_system] + edge_types + [min_confidence]).fetchall()

        causal_edges = []
        for row in rows:
            (src_id, src_sys, tgt_id, tgt_sys,
             rel, conf, eff_size, eff_unit,
             study, primary_src) = row

            # Get labels and tiers for both nodes
            src_node = (self._get_node(src_id, src_sys)
                        or self._get_node_by_id_only(src_id))
            tgt_node = (self._get_node(tgt_id, tgt_sys)
                        or self._get_node_by_id_only(tgt_id))

            src_label = src_node[2] if src_node else src_id
            tgt_label = tgt_node[2] if tgt_node else tgt_id
            src_tier  = src_node[3] if src_node else 0
            tgt_tier  = tgt_node[3] if tgt_node else 0

            causal_edges.append(CausalEdge(
                source_id      = src_id,
                source_system  = src_sys,
                source_label   = src_label[:80],
                rel_type       = rel,
                target_id      = tgt_id,
                target_system  = tgt_sys,
                target_label   = tgt_label[:80],
                confidence     = conf,
                effect_size    = eff_size,
                effect_unit    = eff_unit,
                study_design   = study or "unknown",
                primary_source = primary_src or "unknown",
                tier_from      = src_tier,
                tier_to        = tgt_tier,
            ))

        return causal_edges

    def _personalize_score(
        self,
        base_score:  float,
        edge:        CausalEdge,
        patient_data: dict,
    ) -> float:
        """
        Boost path score if patient has data confirming this edge is active.

        patient_data keys:
          "variants":    list of rsIDs the patient has
          "metabolites": list of HMDB IDs with abnormal levels
          "diseases":    list of confirmed diagnoses
          "genes":       list of NCBI Gene IDs
        """
        boost = 1.0

        # Boost if patient has this specific variant
        variants = patient_data.get("variants", [])
        if edge.source_id in variants or edge.target_id in variants:
            boost *= 1.5

        # Boost if patient has this metabolite abnormal
        metabolites = patient_data.get("metabolites", [])
        if edge.source_id in metabolites or edge.target_id in metabolites:
            boost *= 1.3

        # Boost if patient has a confirmed related diagnosis
        diseases = patient_data.get("diseases", [])
        if edge.source_id in diseases or edge.target_id in diseases:
            boost *= 1.4

        return min(base_score * boost, 0.999)  # cap at <1.0
