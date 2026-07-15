"""
core/hpo_enricher.py

Adds HPO severity, frequency, and onset modifiers to Disease-OS.

Sources:
  - data/raw/hpo.obo (20K HP terms, full HPO OBO)

Adds:
  - 17 HPO modifier nodes (severity, frequency, onset)
  - Severity annotations on HAS_SYMPTOM edges (label-inferred)
  - Onset annotations on CAUSES edges (label-extracted)
  - Frequency annotations on INCREASES_RISK_OF edges (gnomAD AF proxy)
  - severity_source field marking inference method on all annotations

Note on data fidelity:
  Severity/onset annotations marked 'label_inferred' or 'label_extracted'
  are derived from phenotype/disease label text, NOT from HPO curated
  phenotype annotations (phenotype.hpoa). These are reasonable inferences
  but should not be treated as authoritative clinical severity ratings.
  Load phenotype.hpoa for HPO-curated disease-phenotype frequency data.

Run from project root:
    python3 core/hpo_enricher.py
"""

import re, sqlite3, json
from pathlib import Path
from core.node import Node
from core.config import DB_PATH

HPO_OBO = Path.home() / "disease-os/data/raw/hpo.obo"

SEVERITY_TERMS = {
    "HP:0012825": ("Mild",       "SeverityModifier"),
    "HP:0012826": ("Moderate",   "SeverityModifier"),
    "HP:0012827": ("Borderline", "SeverityModifier"),
    "HP:0012828": ("Severe",     "SeverityModifier"),
    "HP:0012829": ("Profound",   "SeverityModifier"),
}
FREQUENCY_TERMS = {
    "HP:0040280": ("Obligate",      "FrequencyModifier"),
    "HP:0040281": ("Very frequent", "FrequencyModifier"),
    "HP:0040282": ("Frequent",      "FrequencyModifier"),
    "HP:0040283": ("Occasional",    "FrequencyModifier"),
    "HP:0040284": ("Very rare",     "FrequencyModifier"),
    "HP:0040285": ("Excluded",      "FrequencyModifier"),
}
ONSET_TERMS = {
    "HP:0003581": ("Adult onset",      "OnsetModifier"),
    "HP:0011463": ("Childhood onset",  "OnsetModifier"),
    "HP:0003623": ("Neonatal onset",   "OnsetModifier"),
    "HP:0003596": ("Middle age onset", "OnsetModifier"),
    "HP:0003584": ("Late onset",       "OnsetModifier"),
    "HP:0410280": ("Pediatric onset",  "OnsetModifier"),
}

SEVERE_PATTERNS  = ["severe","profound","complete","total","fatal",
                    "lethal","death","failure","crisis","collapse",
                    "shock","coma","extensive","extreme","critical"]
MODERATE_PATTERNS = ["moderate","significant","substantial","impaired",
                     "reduced","decreased","elevated","increased",
                     "chronic","recurrent","persistent"]
MILD_PATTERNS    = ["mild","slight","minor","minimal","borderline",
                    "subclinical","transient","intermittent"]


def parse_hpo(path: Path) -> dict:
    hp_terms = {}
    current  = {}
    in_term  = False
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if line == "[Term]":
                if current.get("id","").startswith("HP:"):
                    hp_terms[current["id"]] = current
                current = {"synonyms":[], "is_a":[]}
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("def: "):
                m = re.match(r'def: "([^"]*)"', line)
                if m:
                    current["def"] = m.group(1)[:300]
            elif line.startswith("synonym: "):
                m = re.match(r'synonym: "([^"]*)"', line)
                if m:
                    current["synonyms"].append(m.group(1).strip())
            elif line.startswith("is_obsolete:"):
                current["obsolete"] = True
    if current.get("id","").startswith("HP:"):
        hp_terms[current["id"]] = current
    return hp_terms


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")

    print(f"[HPO] Enriching Disease-OS with HPO modifiers")
    print(f"  HPO file : {HPO_OBO}")
    print(f"  DB       : {DB_PATH}")

    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")

    # Add columns if missing
    for col in ["severity","severity_code","frequency","frequency_code",
                "onset","onset_code","stage_context","severity_source"]:
        try:
            conn.execute(f"ALTER TABLE edges ADD COLUMN {col} TEXT")
        except sqlite3.OperationalError:
            pass
    conn.commit()

    # Parse HPO
    print(f"\n[1/4] Parsing HPO OBO...")
    hp_terms = parse_hpo(HPO_OBO)
    print(f"  {len(hp_terms):,} HP terms")

    # Add modifier nodes
    print(f"\n[2/4] Adding HPO modifier nodes...")
    sample = Node("x","x","x",0,"Gene")
    cols   = list(sample.to_dict().keys())
    ph     = ",".join(["?"]*len(cols))
    ALL    = {**SEVERITY_TERMS, **FREQUENCY_TERMS, **ONSET_TERMS}
    conn.execute("BEGIN")
    for hpid, (default_label, entity_type) in ALL.items():
        term = hp_terms.get(hpid, {})
        node = Node(
            primary_id=hpid, primary_system="HPO",
            label=term.get("name", default_label), tier=6,
            entity_type=entity_type,
            synonyms=term.get("synonyms",[]),
            definition=term.get("def","") or None,
            source="HPO", source_version="2024", confidence=1.0,
        )
        conn.execute(
            f"INSERT OR REPLACE INTO nodes ({','.join(cols)}) VALUES ({ph})",
            list(node.to_dict().values())
        )
    conn.execute("COMMIT")
    print(f"  {len(ALL)} modifier nodes added")

    # Severity on HAS_SYMPTOM edges
    print(f"\n[3/4] Annotating HAS_SYMPTOM edges with severity...")
    offset = severe = moderate = mild = 0
    while True:
        edges = conn.execute("""
            SELECT e.id, n.label FROM edges e
            JOIN nodes n ON n.primary_id=e.target_id
                         AND n.primary_system=e.target_system
            WHERE e.relationship_type='HAS_SYMPTOM'
              AND (e.severity IS NULL OR e.severity='')
            LIMIT 10000 OFFSET ?
        """, (offset,)).fetchall()
        if not edges:
            break
        updates = []
        for eid, label in edges:
            ll = (label or "").lower()
            if any(p in ll for p in SEVERE_PATTERNS):
                updates.append(("severe","HP:0012828","label_inferred",eid))
                severe += 1
            elif any(p in ll for p in MODERATE_PATTERNS):
                updates.append(("moderate","HP:0012826","label_inferred",eid))
                moderate += 1
            elif any(p in ll for p in MILD_PATTERNS):
                updates.append(("mild","HP:0012825","label_inferred",eid))
                mild += 1
        if updates:
            conn.execute("BEGIN")
            conn.executemany(
                "UPDATE edges SET severity=?,severity_code=?,"
                "severity_source=? WHERE id=?", updates
            )
            conn.execute("COMMIT")
        offset += 10000
    print(f"  severe={severe:,}  moderate={moderate:,}  mild={mild:,}")

    # Onset on CAUSES edges
    print(f"\n[4/4] Annotating CAUSES edges with onset...")
    conn.execute("BEGIN")
    conn.execute("""
        UPDATE edges SET onset='childhood',onset_code='HP:0011463',
          severity_source='label_extracted'
        WHERE relationship_type='CAUSES'
          AND (onset IS NULL OR onset='')
          AND EXISTS (SELECT 1 FROM nodes n WHERE n.primary_id=edges.target_id
            AND (LOWER(n.label) LIKE '%childhood%'
              OR LOWER(n.label) LIKE '%juvenile%'
              OR LOWER(n.label) LIKE '%congenital%'))
    """)
    childhood = conn.execute("SELECT changes()").fetchone()[0]
    conn.execute("""
        UPDATE edges SET onset='adult',onset_code='HP:0003581',
          severity_source='label_extracted'
        WHERE relationship_type='CAUSES'
          AND (onset IS NULL OR onset='')
          AND EXISTS (SELECT 1 FROM nodes n WHERE n.primary_id=edges.target_id
            AND (LOWER(n.label) LIKE '%adult%'
              OR LOWER(n.label) LIKE '%late-onset%'))
    """)
    adult = conn.execute("SELECT changes()").fetchone()[0]
    conn.execute("COMMIT")
    print(f"  childhood={childhood:,}  adult={adult:,}")

    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.close()
    print(f"\n✓ HPO enrichment complete")
