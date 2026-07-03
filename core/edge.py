"""
edge.py — Universal Edge for the Disease-OS knowledge graph.
All relationship types normalised to canonical vocabulary.
Provenance split into primary_source vs imported_via to prevent
circular overcounting across databases.
"""
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Optional

CANONICAL_RELATIONSHIP_TYPES = {
    "ISA", "PART_OF",
    "CAUSES", "CONTRIBUTES_TO", "PROTECTS_AGAINST",
    "UPREGULATES", "DOWNREGULATES", "ACTIVATES", "INHIBITS",
    "BINDS", "CATALYZES", "ENCODES",
    "INTERACTS_WITH",
    "INDICATES", "TREATS", "TARGETED_BY", "DIAGNOSED_BY", "HAS_SYMPTOM",
    "ASSOCIATED_WITH",
    "INCREASES_RISK_OF", "DECREASES_RISK_OF",
    "LOCATED_IN", "EXPRESSED_IN", "DETECTED_IN",
    "PRECEDES",
}


@dataclass
class Edge:
    # ── Endpoints ──────────────────────────────────────────────────────
    source_id:     str
    source_system: str
    target_id:     str
    target_system: str

    # ── Relationship ───────────────────────────────────────────────────
    relationship_type:        str
    source_relationship_type: str   = ""

    # ── Quantitative ───────────────────────────────────────────────────
    effect_size:  Optional[float] = None
    effect_unit:  Optional[str]   = None
    direction:    Optional[str]   = None

    # ── Confidence ─────────────────────────────────────────────────────
    confidence: float = 0.5

    # ── Feedback loop ──────────────────────────────────────────────────
    feedback:       bool = False
    feedback_notes: str  = ""

    # ── Provenance (split to prevent circular overcounting) ────────────
    primary_source: str = "unknown"   # original study/DB that generated claim
    imported_via:   str = "unknown"   # file/DB we actually read it from
    study_design:   str = "unknown"   # RCT / cohort / mechanistic / curated

    # ── Context ────────────────────────────────────────────────────────
    population_context: str          = "general"
    tissue_context:     str          = "systemic"
    species:            str          = "human"
    typical_latency:    Optional[str] = None

    # ── Versioning ─────────────────────────────────────────────────────
    source_version: str = "unknown"
    loaded_at:      str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    def __post_init__(self):
        if self.relationship_type not in CANONICAL_RELATIONSHIP_TYPES:
            raise ValueError(
                f"'{self.relationship_type}' not a canonical type. "
                f"Use ASSOCIATED_WITH or add to CANONICAL_RELATIONSHIP_TYPES."
            )

    def to_dict(self) -> dict:
        return {
            "source_id":                self.source_id,
            "source_system":            self.source_system,
            "target_id":                self.target_id,
            "target_system":            self.target_system,
            "relationship_type":        self.relationship_type,
            "source_relationship_type": self.source_relationship_type,
            "effect_size":              self.effect_size,
            "effect_unit":              self.effect_unit,
            "direction":                self.direction,
            "confidence":               self.confidence,
            "feedback":                 int(self.feedback),
            "feedback_notes":           self.feedback_notes,
            "primary_source":           self.primary_source,
            "imported_via":             self.imported_via,
            "study_design":             self.study_design,
            "population_context":       self.population_context,
            "tissue_context":           self.tissue_context,
            "species":                  self.species,
            "typical_latency":          self.typical_latency,
            "source_version":           self.source_version,
            "loaded_at":                self.loaded_at,
        }
