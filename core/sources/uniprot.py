"""
core/sources/uniprot.py

UniProt Swiss-Prot human proteome adapter for Disease-OS.

Input : data/raw/uniprot/uniprot_sprot_human.dat (453MB)
Format: UniProt flat file (.dat) with two-letter field codes

What UniProt adds:
  - ~20K manually curated human protein nodes (primary_id = UniProt accession)
  - Resolves 393K staged HMDB→UniProt_pending edges
  - Gene→Protein ENCODES edges (links NCBI Gene nodes to proteins)
  - Subcellular localization, function, active sites
  - Disease associations from Swiss-Prot curation
  - Cross-references: Ensembl, NCBI Gene, PDB, OMIM, Reactome

Key field codes parsed:
  ID  entry name, review status, sequence length
  AC  primary accession + secondary accessions
  DE  protein names (RecName = canonical, AltName = alternatives)
  GN  gene name(s)
  DR  database cross-references
  CC  comments: function, subcellular location, disease
  FT  feature annotations: active site, binding site, domain
"""

import sys
import re
from pathlib import Path
from typing import Generator

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from core.node import Node
from core.edge import Edge
from core.sources.base_source import BaseSource
from core.config import SOURCE_VERSIONS

RAW_DIR     = Path.home() / "disease-os/data/raw/uniprot"
UNIPROT_DAT = RAW_DIR / "uniprot_sprot_human.dat"


def _parse_dat(path: Path) -> Generator[dict, None, None]:
    """
    Stream UniProt flat file one entry at a time.
    Yields one dict per protein entry.
    """
    entry = {}
    lines = {"AC": [], "DE": [], "GN": [], "DR": [],
             "CC": [], "FT": [], "OX": []}
    current = None

    with open(path, encoding="utf-8") as f:
        for line in f:
            code = line[:2].strip()
            data = line[5:].rstrip()

            if code == "//":
                # End of entry — yield and reset
                if entry:
                    entry["_lines"] = dict(lines)
                    yield entry
                entry = {}
                lines = {"AC": [], "DE": [], "GN": [],
                         "DR": [], "CC": [], "FT": [], "OX": []}
                current = None
                continue

            if code == "ID":
                parts = data.split()
                entry["_id_name"] = parts[0] if parts else ""
                entry["_reviewed"] = "Reviewed" in data
                entry["_length"] = int(parts[-2]) if parts[-1] == "AA." else 0

            elif code in lines:
                lines[code].append(data)

    # Yield last entry if file doesn't end with //
    if entry:
        entry["_lines"] = dict(lines)
        yield entry


def _extract_entry(raw: dict) -> dict | None:
    """
    Parse a raw entry dict into a structured protein record.
    """
    ls = raw.get("_lines", {})

    # ── Accessions ──────────────────────────────────────────────────────
    ac_text = " ".join(ls.get("AC", []))
    accessions = [a.strip().rstrip(";") for a in ac_text.split(";")
                  if a.strip()]
    if not accessions:
        return None
    primary_ac = accessions[0]
    secondary_acs = accessions[1:]

    # ── Protein name ────────────────────────────────────────────────────
    label = ""
    alt_names = []
    for line in ls.get("DE", []):
        if "RecName: Full=" in line:
            m = re.search(r'RecName: Full=([^;{]+)', line)
            if m and not label:
                label = m.group(1).strip()
        elif "AltName: Full=" in line:
            m = re.search(r'AltName: Full=([^;{]+)', line)
            if m:
                alt_names.append(m.group(1).strip())
        elif "Short=" in line and not any("RecName" in l
                                          for l in ls.get("DE", [])):
            m = re.search(r'Short=([^;{]+)', line)
            if m:
                alt_names.append(m.group(1).strip())

    if not label:
        label = raw.get("_id_name", primary_ac)

    # ── Gene name ────────────────────────────────────────────────────────
    gene_name = ""
    gene_synonyms = []
    for line in ls.get("GN", []):
        if "Name=" in line and not gene_name:
            m = re.search(r'Name=([^;{,]+)', line)
            if m:
                gene_name = m.group(1).strip()
        syns = re.findall(r'Synonyms=([^;{]+)', line)
        for s in syns:
            gene_synonyms.extend([x.strip() for x in s.split(",")])

    # ── Cross-references ────────────────────────────────────────────────
    xrefs = {}
    ncbi_gene_ids = []
    ensembl_ids   = []
    omim_ids      = []
    pdb_ids       = []
    reactome_ids  = []

    for line in ls.get("DR", []):
        parts = [p.strip().rstrip(".") for p in line.split(";")]
        if not parts:
            continue
        db = parts[0]

        if db == "GeneID" and len(parts) > 1:
            ncbi_gene_ids.append(parts[1])
        elif db == "Ensembl" and len(parts) > 1:
            ensembl_ids.append(parts[1])
        elif db == "MIM" and len(parts) > 1:
            omim_ids.append(parts[1])
        elif db == "PDB" and len(parts) > 1:
            pdb_ids.append(parts[1])
        elif db == "Reactome" and len(parts) > 1:
            reactome_ids.append(parts[1])
        elif db == "HGNC" and len(parts) > 1:
            xrefs["HGNC"] = parts[1]

    if ncbi_gene_ids:
        xrefs["NCBI_Gene"] = ncbi_gene_ids[0]
        if len(ncbi_gene_ids) > 1:
            xrefs["NCBI_Gene_alt"] = ",".join(ncbi_gene_ids[1:])
    if ensembl_ids:
        xrefs["Ensembl"] = ensembl_ids[0]
    if omim_ids:
        xrefs["OMIM"] = omim_ids[0]
    if reactome_ids:
        xrefs["Reactome"] = reactome_ids[0]

    # ── Function and subcellular location (from CC) ──────────────────────
    function    = ""
    subcell     = ""
    disease_txt = ""
    in_function = in_subcell = in_disease = False

    for line in ls.get("CC", []):
        if "-!- FUNCTION:" in line:
            function    = line.replace("-!- FUNCTION:", "").strip()
            in_function = True
            in_subcell  = in_disease = False
        elif "-!- SUBCELLULAR LOCATION:" in line:
            subcell    = line.replace("-!- SUBCELLULAR LOCATION:", "").strip()
            in_subcell = True
            in_function = in_disease = False
        elif "-!- DISEASE:" in line:
            disease_txt = line.replace("-!- DISEASE:", "").strip()
            in_disease  = True
            in_function = in_subcell = False
        elif line.startswith("    ") and not "-!-" in line:
            if in_function:
                function    += " " + line.strip()
            elif in_subcell:
                subcell     += " " + line.strip()
            elif in_disease:
                disease_txt += " " + line.strip()
        elif "-!-" in line:
            in_function = in_subcell = in_disease = False

    # ── Active sites and binding sites (from FT) ─────────────────────────
    active_sites = []
    binding_sites = []
    for line in ls.get("FT", []):
        if line.startswith("ACTIVE_SITE"):
            m = re.search(r'/note="([^"]+)"', line)
            if m:
                active_sites.append(m.group(1))
        elif line.startswith("BINDING"):
            m = re.search(r'/note="([^"]+)"', line)
            if m:
                binding_sites.append(m.group(1))

    return {
        "primary_ac":     primary_ac,
        "secondary_acs":  secondary_acs,
        "label":          label[:200],
        "alt_names":      alt_names[:10],
        "gene_name":      gene_name,
        "gene_synonyms":  gene_synonyms[:10],
        "ncbi_gene_ids":  ncbi_gene_ids,
        "ensembl_ids":    ensembl_ids,
        "xrefs":          xrefs,
        "function":       function[:500],
        "subcell":        subcell[:200],
        "disease_txt":    disease_txt[:300],
        "active_sites":   active_sites[:5],
        "binding_sites":  binding_sites[:5],
        "length":         raw.get("_length", 0),
        "reviewed":       raw.get("_reviewed", True),
        "pdb_ids":        pdb_ids[:5],
        "reactome_ids":   reactome_ids[:5],
    }


class UniProtSource(BaseSource):
    """
    Loads UniProt Swiss-Prot human proteins into Disease-OS.

    Creates:
      - Protein nodes (primary_system = "UniProt")
      - Gene→Protein ENCODES edges (NCBI_Gene → UniProt)
      - Protein→Pathway PART_OF edges (UniProt → Reactome)

    Also resolves:
      - 393K staged HMDB→UniProt_pending edges
        (done separately via resolve_hmdb_pending)

    Usage:
        source = UniProtSource()
        source.load_into(graph_store)
        source.resolve_hmdb_pending(conn)
    """

    source_name    = "UNIPROT_SWISSPROT"
    source_version = SOURCE_VERSIONS.get("UNIPROT", "2026_06")

    def nodes(self) -> Generator[Node, None, None]:
        print(f"[UniProt] Streaming {UNIPROT_DAT.name}...")
        n = 0
        for raw in _parse_dat(UNIPROT_DAT):
            entry = _extract_entry(raw)
            if not entry:
                continue

            synonyms = list(set(
                entry["alt_names"] +
                entry["gene_synonyms"] +
                entry["secondary_acs"][:5]
            ))[:20]

            xrefs = dict(entry["xrefs"])
            for ac in entry["secondary_acs"][:5]:
                xrefs[f"UniProt_secondary_{ac}"] = ac

            properties = {}
            if entry["function"]:
                properties["function"] = entry["function"]
            if entry["subcell"]:
                properties["subcellular_location"] = entry["subcell"]
            if entry["active_sites"]:
                properties["active_sites"] = entry["active_sites"]
            if entry["binding_sites"]:
                properties["binding_sites"] = entry["binding_sites"]
            if entry["length"]:
                properties["sequence_length"] = entry["length"]
            if entry["gene_name"]:
                properties["gene_name"] = entry["gene_name"]
            if entry["pdb_ids"]:
                properties["pdb_ids"] = entry["pdb_ids"]

            yield Node(
                primary_id     = entry["primary_ac"],
                primary_system = "UniProt",
                label          = entry["label"],
                tier           = 1,
                entity_type    = "Protein",
                xrefs          = xrefs,
                synonyms       = synonyms,
                definition     = entry["function"][:500] or None,
                properties     = properties,
                source         = self.source_name,
                source_version = self.source_version,
                confidence     = 1.0,  # Swiss-Prot is manually curated
            )
            n += 1
            if n % 5000 == 0:
                print(f"  {n:,} proteins yielded...")

        print(f"[UniProt] {n:,} protein nodes total")

    def edges(self) -> Generator[Edge, None, None]:
        """
        Yield Gene→Protein ENCODES edges and Protein→Pathway PART_OF edges.
        """
        print(f"[UniProt] Generating edges...")
        n_encodes = n_pathway = 0

        for raw in _parse_dat(UNIPROT_DAT):
            entry = _extract_entry(raw)
            if not entry:
                continue

            primary_ac  = entry["primary_ac"]
            ncbi_gene_id = (entry["ncbi_gene_ids"][0]
                            if entry["ncbi_gene_ids"] else None)

            # Gene → ENCODES → Protein
            if ncbi_gene_id:
                try:
                    yield Edge(
                        source_id                = ncbi_gene_id,
                        source_system            = "NCBI_Gene",
                        target_id                = primary_ac,
                        target_system            = "UniProt",
                        relationship_type        = "ENCODES",
                        source_relationship_type = "uniprot_gene_protein",
                        confidence               = 1.0,
                        primary_source           = f"UniProt_{self.source_version}",
                        imported_via             = f"UniProt_SwissProt_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                    n_encodes += 1
                except ValueError:
                    pass

            # Protein → PART_OF → Reactome Pathway
            for reactome_id in entry["reactome_ids"]:
                try:
                    yield Edge(
                        source_id                = primary_ac,
                        source_system            = "UniProt",
                        target_id                = reactome_id,
                        target_system            = "Reactome",
                        relationship_type        = "PART_OF",
                        source_relationship_type = "uniprot_reactome_pathway",
                        confidence               = 0.90,
                        primary_source           = f"UniProt_{self.source_version}",
                        imported_via             = f"UniProt_SwissProt_{self.source_version}",
                        study_design             = "curated",
                        source_version           = self.source_version,
                    )
                    n_pathway += 1
                except ValueError:
                    pass

        print(f"[UniProt] {n_encodes:,} ENCODES edges  |  "
              f"{n_pathway:,} PART_OF edges")

    def normalize_confidence(self, raw_value=None) -> float:
        return 1.0  # Swiss-Prot is manually curated — highest confidence

    def resolve_hmdb_pending(self, conn) -> dict:
        """
        Resolve the 393K staged HMDB→UniProt_pending edges
        now that UniProt nodes exist.
        Simply updates target_system from 'UniProt_pending' to 'UniProt'.
        Only updates edges where the target UniProt ID now exists as a node.
        """
        print(f"\n[UniProt] Resolving staged HMDB→UniProt_pending edges...")

        # Count how many can be resolved
        resolvable = conn.execute("""
            SELECT COUNT(*) FROM edges e
            WHERE e.target_system = 'UniProt_pending'
              AND EXISTS (
                SELECT 1 FROM nodes n
                WHERE n.primary_id = e.target_id
                  AND n.primary_system = 'UniProt'
              )
        """).fetchone()[0]

        total = conn.execute(
            "SELECT COUNT(*) FROM edges WHERE target_system='UniProt_pending'"
        ).fetchone()[0]

        print(f"  Total staged    : {total:,}")
        print(f"  Resolvable now  : {resolvable:,}")
        print(f"  Still pending   : {total - resolvable:,} "
              f"(protein not in Swiss-Prot)")

        # Resolve what we can
        conn.execute("BEGIN")
        conn.execute("""
            UPDATE edges
            SET target_system = 'UniProt'
            WHERE target_system = 'UniProt_pending'
              AND EXISTS (
                SELECT 1 FROM nodes n
                WHERE n.primary_id = edges.target_id
                  AND n.primary_system = 'UniProt'
              )
        """)
        conn.execute("COMMIT")

        still_pending = conn.execute(
            "SELECT COUNT(*) FROM edges WHERE target_system='UniProt_pending'"
        ).fetchone()[0]

        print(f"  Resolved        : {total - still_pending:,}")
        print(f"  Remaining staged: {still_pending:,}")

        return {
            "total":     total,
            "resolved":  total - still_pending,
            "pending":   still_pending,
        }
