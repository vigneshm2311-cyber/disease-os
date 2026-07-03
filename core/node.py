"""
node.py — Universal Node for the Disease-OS knowledge graph.
Primary ID is entity-type specific and open.
UMLS CUI lives in xrefs as a cross-reference hub, not the spine.
Insurance codes are explicit first-class fields.
"""
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Optional
import json


@dataclass
class Node:
    # ── Core identity ──────────────────────────────────────────────────
    primary_id:     str
    primary_system: str
    label:          str
    tier:           int
    entity_type:    str

    # ── Cross-references ───────────────────────────────────────────────
    xrefs: dict = field(default_factory=dict)
    # {
    #   "UMLS_CUI": "C0011860",   ← hub, not primary
    #   "SNOMED":   "44054006",
    #   "MONDO":    "MONDO:0005148",
    #   "OMIM":     "125853",
    # }

    # ── Insurance-critical codes (explicit, not buried in xrefs) ───────
    icd10_code:  Optional[str] = None   # diagnosis billing
    icd11_code:  Optional[str] = None   # forward compatibility
    snomed_code: Optional[str] = None   # clinical finding
    loinc_code:  Optional[str] = None   # lab test
    rxnorm_cui:  Optional[str] = None   # pharmacy claims
    cpt_code:    Optional[str] = None   # procedure billing
    hcc_code:    Optional[str] = None   # risk adjustment
    ndc_code:    Optional[str] = None   # drug packaging

    # ── Biological metadata ────────────────────────────────────────────
    synonyms:   list  = field(default_factory=list)
    definition: Optional[str] = None
    properties: dict  = field(default_factory=dict)

    # ── Provenance ─────────────────────────────────────────────────────
    source:         str   = "unknown"
    source_version: str   = "unknown"
    confidence:     float = 1.0
    loaded_at:      str   = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )

    def umls_cui(self) -> Optional[str]:
        return self.xrefs.get("UMLS_CUI")

    def has_insurance_codes(self) -> bool:
        return any([self.icd10_code, self.rxnorm_cui,
                    self.loinc_code, self.cpt_code])

    def to_dict(self) -> dict:
        return {
            "primary_id":     self.primary_id,
            "primary_system": self.primary_system,
            "label":          self.label,
            "tier":           self.tier,
            "entity_type":    self.entity_type,
            "xrefs":          json.dumps(self.xrefs),
            "icd10_code":     self.icd10_code,
            "icd11_code":     self.icd11_code,
            "snomed_code":    self.snomed_code,
            "loinc_code":     self.loinc_code,
            "rxnorm_cui":     self.rxnorm_cui,
            "cpt_code":       self.cpt_code,
            "hcc_code":       self.hcc_code,
            "ndc_code":       self.ndc_code,
            "synonyms":       json.dumps(self.synonyms),
            "definition":     self.definition,
            "properties":     json.dumps(self.properties),
            "source":         self.source,
            "source_version": self.source_version,
            "confidence":     self.confidence,
            "loaded_at":      self.loaded_at,
        }
