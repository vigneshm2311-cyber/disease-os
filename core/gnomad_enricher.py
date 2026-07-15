"""
core/gnomad_enricher.py

Enriches variant nodes with gnomAD r4 population allele frequencies.

Strategy:
  1. Extract GRCh38 coordinates from ClinVar variant_summary.txt
  2. Backfill GWAS Catalog GRCh37 coordinates (note: these need liftover
     to GRCh38 before gnomAD queries will work — currently stored as-is)
  3. Query gnomAD GraphQL API per variant using region endpoint
     (region query matches rsID in response, bypassing allele mismatch)

Known limitation:
  GWAS Catalog positions are GRCh37. gnomAD r4 uses GRCh38.
  Liftover required before GWAS variants can be queried.
  ClinVar GRCh38 positions work correctly.

Current state:
  - 1,045,484 variants have GRCh38 coords from ClinVar
  - 34 variants have gnomAD AF data (pilot batch)
  - GWAS liftover pending (requires UCSC liftOver tool)

Run from project root:
    python3 core/gnomad_enricher.py
"""

import csv, sqlite3, json, time, urllib.request
from pathlib import Path
from core.config import DB_PATH

CLINVAR  = Path.home() / "disease-os/data/raw/clinvar/variant_summary.txt"
GWAS_TSV = Path.home() / "disease-os/data/raw/gwas_catalog/associations.tsv"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"


def extract_clinvar_coords(clinvar_path: Path) -> dict:
    """Extract GRCh38 VCF coordinates from ClinVar variant_summary.txt."""
    coord_map = {}
    with open(clinvar_path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 34 or row[16] != "GRCh38":
                continue
            rs_raw  = row[9].strip()
            if not rs_raw or rs_raw in ("-1","0",""):
                continue
            rsid    = f"rs{rs_raw}"
            chrom   = row[18].strip()
            pos_vcf = row[31].strip()
            ref_vcf = row[32].strip()
            alt_vcf = row[33].strip()
            if (chrom and pos_vcf and ref_vcf and alt_vcf
                    and ref_vcf != "na" and alt_vcf != "na"
                    and pos_vcf != "0"):
                coord_map[rsid] = (chrom, pos_vcf, ref_vcf, alt_vcf)
    return coord_map


def extract_gwas_coords(gwas_path: Path) -> dict:
    """
    Extract GRCh37 positions from GWAS Catalog.
    NOTE: These are GRCh37 — need liftover before gnomAD query.
    """
    coord_map = {}
    with open(gwas_path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 25:
                continue
            chrom    = row[11].strip()
            pos      = row[12].strip()
            snp_raw  = row[23].strip()
            snps_col = row[21].strip()
            if not chrom or not pos or pos == "NA":
                continue
            rsid = None
            if snps_col.startswith("rs"):
                rsid = snps_col.split(";")[0].split(",")[0].strip()
            elif snp_raw and snp_raw.isdigit():
                rsid = f"rs{snp_raw}"
            if rsid and rsid not in coord_map:
                coord_map[rsid] = (chrom, pos)
    return coord_map


def query_gnomad_region(chrom: str, pos: int) -> list:
    """Query gnomAD for all variants at a genomic position (GRCh38)."""
    query = f"""
    {{
      region(chrom: "{chrom}", start: {pos-1}, stop: {pos+1},
             reference_genome: GRCh38) {{
        variants(dataset: gnomad_r4) {{
          variant_id rsids
          genome {{ ac an af
            populations {{ id ac an af }}
          }}
        }}
      }}
    }}
    """
    req = urllib.request.Request(
        GNOMAD_API,
        data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json"},
    )
    try:
        resp = urllib.request.urlopen(req, timeout=10)
        data = json.loads(resp.read())
        return (data.get("data",{})
                    .get("region",{})
                    .get("variants",[]))
    except Exception:
        return []


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    print(f"[gnomAD] Variant coordinate extractor + gnomAD enricher")
    print(f"  DB: {DB_PATH}")

    conn = sqlite3.connect(str(DB_PATH))
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA cache_size=-256000;")

    # Step 1: Extract ClinVar GRCh38 coords
    print(f"\n[1/3] Extracting ClinVar GRCh38 coordinates...")
    clinvar_coords = extract_clinvar_coords(CLINVAR)
    print(f"  {len(clinvar_coords):,} rsIDs with GRCh38 coords")

    # Step 2: Backfill variant nodes
    print(f"\n[2/3] Updating variant nodes with coordinates...")
    updated = 0
    conn.execute("BEGIN")
    for rsid, (chrom, pos, ref, alt) in clinvar_coords.items():
        row = conn.execute(
            "SELECT id, properties FROM nodes "
            "WHERE primary_id=? AND primary_system='dbSNP_rsID'",
            (rsid,)
        ).fetchone()
        if row:
            props = json.loads(row[1] or "{}")
            if not props.get("pos_vcf") or props.get("pos_vcf") in ("0","None"):
                props.update({
                    "chromosome": chrom,
                    "pos_vcf":    pos,
                    "ref":        ref,
                    "alt":        alt,
                    "gnomad_id":  f"{chrom}-{pos}-{ref}-{alt}",
                })
                conn.execute(
                    "UPDATE nodes SET properties=? WHERE id=?",
                    (json.dumps(props), row[0])
                )
                updated += 1
            if updated % 50_000 == 0 and updated > 0:
                conn.execute("COMMIT")
                conn.execute("BEGIN")
                print(f"  {updated:,} nodes updated...")
    conn.execute("COMMIT")
    print(f"  {updated:,} variant nodes updated")

    # Step 3: gnomAD enrichment for clinical variants
    print(f"\n[3/3] Querying gnomAD for high-confidence clinical variants...")
    print(f"  NOTE: Only ClinVar GRCh38 variants are queried.")
    print(f"  GWAS variants need GRCh37→GRCh38 liftover first.")

    clinical = conn.execute("""
        SELECT DISTINCT n.id, n.primary_id,
               json_extract(n.properties,'$.chromosome') as chrom,
               json_extract(n.properties,'$.pos_vcf') as pos,
               n.properties
        FROM nodes n
        JOIN edges e ON e.source_id=n.primary_id
                     AND e.source_system=n.primary_system
        WHERE n.entity_type='Variant'
          AND n.primary_system='dbSNP_rsID'
          AND e.relationship_type IN ('CAUSES','CONTRIBUTES_TO',
                                      'INCREASES_RISK_OF')
          AND json_extract(n.properties,'$.gnomad_af_global') IS NULL
          AND json_extract(n.properties,'$.chromosome') IS NOT NULL
          AND json_extract(n.properties,'$.pos_vcf') NOT IN ('0','None')
          AND json_extract(n.properties,'$.ref') != 'na'
        LIMIT 5000
    """).fetchall()

    print(f"  Variants to query: {len(clinical):,}")
    enriched = not_found = 0
    batch    = []

    for i, (node_id, rsid, chrom, pos, props_json) in enumerate(clinical):
        try:
            pos_int = int(pos)
        except:
            not_found += 1
            continue

        variants = query_gnomad_region(str(chrom), pos_int)
        for v in variants:
            if rsid in (v.get("rsids") or []):
                genome = v.get("genome") or {}
                af     = genome.get("af")
                if af is not None:
                    pop_afs = {
                        p["id"].lower(): round(p["af"], 6)
                        for p in (genome.get("populations") or [])
                        if p.get("af", 0) > 0
                    }
                    props = json.loads(props_json or "{}")
                    props.update({
                        "gnomad_id":        v["variant_id"],
                        "gnomad_af_global": round(af, 6),
                        "gnomad_ac":        genome.get("ac"),
                        "gnomad_an":        genome.get("an"),
                        "gnomad_pop_afs":   pop_afs,
                        "gnomad_dataset":   "gnomad_r4",
                    })
                    batch.append((json.dumps(props), node_id))
                    enriched += 1
                break
        else:
            not_found += 1

        if len(batch) >= 100:
            conn.execute("BEGIN")
            conn.executemany(
                "UPDATE nodes SET properties=? WHERE id=?", batch
            )
            conn.execute("COMMIT")
            batch = []

        if (i + 1) % 500 == 0:
            print(f"  {i+1:,}/{len(clinical)} | "
                  f"enriched={enriched} | not_found={not_found}")
        time.sleep(0.25)

    if batch:
        conn.execute("BEGIN")
        conn.executemany(
            "UPDATE nodes SET properties=? WHERE id=?", batch
        )
        conn.execute("COMMIT")

    total = conn.execute("""
        SELECT COUNT(*) FROM nodes
        WHERE json_extract(properties,'$.gnomad_af_global') IS NOT NULL
    """).fetchone()[0]

    print(f"\n  Enriched this run : {enriched:,}")
    print(f"  Total with gnomAD : {total:,}")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.close()
    print(f"\n✓ gnomAD enrichment complete")
