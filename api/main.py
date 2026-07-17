"""
api/main.py — Disease-OS REST API v2

Changes from v1:
  - API key authentication (X-API-Key header)
  - In-memory LRU cache for causal traces (1 hour TTL)
  - Pagination on /causes and /search
  - TREATS edges ranked first in /drugs
  - New: GET /drug/{rxnorm}/targets  (ChEMBL TARGETS)
  - New: GET /drug/{rxnorm}/mechanism (MOA text)
  - New: GET /gene/{id}/drugs  (drugs targeting this gene's protein)

Run:
    cd ~/disease-os
    DISEASE_OS_API_KEY=your-secret-key uvicorn api.main:app --port 8000 --reload

Docs: http://localhost:8000/docs
"""

import sys, json, os, time, sqlite3
from pathlib import Path
from typing import Optional
from functools import lru_cache
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from fastapi import FastAPI, HTTPException, Query, Depends, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.security.api_key import APIKeyHeader
from pydantic import BaseModel

from core.config import DB_PATH
from core.causal_engine import CausalEngine

# ── App ────────────────────────────────────────────────────────
app = FastAPI(
    title       = "Disease-OS API",
    description = """
Universal causal disease modeling API — v2.

Requires X-API-Key header for all endpoints.
Contact the Disease-OS team for an API key.

**Coverage:** 3,057,056 nodes · 10,734,822 edges · 14 databases
**GitHub:** github.com/vigneshm2311-cyber/disease-os
    """,
    version  = "2.0.0",
    docs_url = "/docs",
    redoc_url= "/redoc",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_methods=["*"], allow_headers=["*"],
)

# ── Auth ───────────────────────────────────────────────────────
API_KEY        = os.environ.get("DISEASE_OS_API_KEY", "dev-key-change-in-production")
api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)

async def require_api_key(key: str = Depends(api_key_header)):
    if not key or key != API_KEY:
        raise HTTPException(
            status_code=403,
            detail="Invalid or missing API key. Pass X-API-Key header."
        )
    return key

# Health endpoint is public — all others require auth
PUBLIC_PATHS = {"/health", "/docs", "/redoc", "/openapi.json"}

@app.middleware("http")
async def auth_middleware(request: Request, call_next):
    if request.url.path in PUBLIC_PATHS:
        return await call_next(request)
    key = request.headers.get("X-API-Key","")
    if key != API_KEY:
        return JSONResponse(
            status_code=403,
            content={"detail": "Invalid or missing API key. Pass X-API-Key header."}
        )
    return await call_next(request)

# ── Rate limiting (simple in-memory) ──────────────────────────
_request_counts: dict = {}
_request_window = 60   # seconds
_request_limit  = 120  # requests per window per key

def check_rate_limit(key: str) -> None:
    now  = time.time()
    data = _request_counts.get(key, {"count": 0, "reset": now + _request_window})
    if now > data["reset"]:
        data = {"count": 0, "reset": now + _request_window}
    data["count"] += 1
    _request_counts[key] = data
    if data["count"] > _request_limit:
        raise HTTPException(
            status_code=429,
            detail=f"Rate limit exceeded. {_request_limit} requests per minute."
        )

# ── Shared resources ───────────────────────────────────────────
_conn   = None
_engine = None
_cache: dict = {}   # simple TTL cache
CACHE_TTL = 3600    # 1 hour

def get_conn() -> sqlite3.Connection:
    global _conn
    if _conn is None:
        _conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
        _conn.execute("PRAGMA cache_size=-256000;")
        _conn.execute("PRAGMA temp_store=MEMORY;")
    return _conn

def get_engine() -> CausalEngine:
    global _engine
    if _engine is None:
        _engine = CausalEngine(str(DB_PATH))
    return _engine

def cache_get(key: str):
    entry = _cache.get(key)
    if entry and time.time() < entry["expires"]:
        return entry["data"]
    return None

def cache_set(key: str, data):
    _cache[key] = {"data": data, "expires": time.time() + CACHE_TTL}

def jload(s):
    try:    return json.loads(s) if s else None
    except: return s

def find_disease(icd10: str):
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, primary_system, label, icd10_code,
               hcc_code, icd11_code, xrefs, synonyms, definition, confidence
        FROM nodes WHERE icd10_code=? ORDER BY confidence DESC LIMIT 1
    """, (icd10,)).fetchone()
    if not row:
        row = conn.execute("""
            SELECT primary_id, primary_system, label, icd10_code,
                   hcc_code, icd11_code, xrefs, synonyms, definition, confidence
            FROM nodes WHERE primary_id=? AND primary_system='ICD-10-CM' LIMIT 1
        """, (icd10,)).fetchone()
    if not row: return None
    cols = ["primary_id","primary_system","label","icd10_code",
            "hcc_code","icd11_code","xrefs","synonyms","definition","confidence"]
    d = dict(zip(cols, row))
    for f in ["xrefs","synonyms"]: d[f] = jload(d[f])
    return d

# ── Request models ─────────────────────────────────────────────
class PatientData(BaseModel):
    variants:    list[str] = []
    metabolites: list[str] = []
    diseases:    list[str] = []

class TraceRequest(BaseModel):
    disease_id:     str
    disease_system: str   = "ICD-10-CM"
    max_depth:      int   = 4
    min_confidence: float = 0.50
    top_n:          int   = 20
    page:           int   = 1
    page_size:      int   = 20
    patient_data:   Optional[PatientData] = None

# ══════════════════════════════════════════════════════════════
# ENDPOINTS
# ══════════════════════════════════════════════════════════════

@app.get("/health", tags=["System"])
def health():
    """System status — public endpoint, no auth required."""
    conn  = get_conn()
    nodes = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    edges = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    treats= conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='TREATS'"
    ).fetchone()[0]
    targets=conn.execute(
        "SELECT COUNT(*) FROM edges WHERE relationship_type='TARGETS'"
    ).fetchone()[0]
    tiers = dict(conn.execute(
        "SELECT tier, COUNT(*) FROM nodes GROUP BY tier ORDER BY tier"
    ).fetchall())
    return {
        "status":  "healthy",
        "version": "2.0.0",
        "graph": {
            "total_nodes":  nodes,
            "total_edges":  edges,
            "treats_edges": treats,
            "targets_edges":targets,
            "nodes_by_tier":tiers,
        },
        "sources": [
            "UMLS 2025AB","Reactome v88","ClinVar 2024",
            "STRING v12","GWAS Catalog 2024","HMDB v5.0",
            "EFO v3.91","Open Targets 26.06","NCBI Gene 2026",
            "UniProt Swiss-Prot 2026","HPO 2024",
            "Orphanet 2026","MONDO 2026","DrugBank 5.1.10","ChEMBL 37",
        ],
        "cache_entries": len(_cache),
    }


@app.get("/disease/{icd10}", tags=["Disease"])
def get_disease(icd10: str):
    """Get disease node by ICD-10 code."""
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}. Try /search?q={icd10}")
    return node


@app.get("/disease/{icd10}/causes", tags=["Disease","Causal"])
def get_causes(
    icd10:          str,
    max_depth:      int   = Query(4,   ge=1, le=6),
    min_confidence: float = Query(0.5, ge=0.0, le=1.0),
    top_n:          int   = Query(20,  ge=1, le=100),
    page:           int   = Query(1,   ge=1),
    page_size:      int   = Query(20,  ge=1, le=100),
):
    """
    Trace root causes of a disease using backward causal BFS.

    Results are paginated. Use page/page_size to navigate.
    Results are cached for 1 hour — repeated queries are instant.
    """
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")

    # Check cache
    cache_key = f"causes:{icd10}:{max_depth}:{min_confidence}"
    cached = cache_get(cache_key)

    if cached:
        paths_data = cached
        from_cache = True
    else:
        engine = get_engine()
        paths  = engine.trace(
            disease_id     = node["primary_id"],
            disease_system = node["primary_system"],
            max_depth      = max_depth,
            min_confidence = min_confidence,
            min_path_score = 0.01,
            top_n          = 500,
        )
        # Deduplicate by root
        best = {}
        for p in paths:
            k = (p.edges[0].source_id, p.edges[0].source_system)
            if k not in best or p.path_score > best[k].path_score:
                best[k] = p
        ranked = sorted(best.values(), key=lambda p: p.path_score, reverse=True)

        # Tier breakdown
        tier_names = {1:"Molecular",2:"Networks",3:"Cellular",4:"Tissue",
                      5:"Systemic",6:"Phenotype",7:"Disease",8:"Behavior",
                      9:"Social",10:"Healthcare"}
        tier_scores = {}
        for p in paths:
            t = p.edges[0].tier_from
            tier_scores[t] = tier_scores.get(t,0) + p.path_score
        total = sum(tier_scores.values()) or 1

        paths_data = {
            "total_paths":    len(paths),
            "unique_roots":   len(best),
            "tier_breakdown": {
                tier_names.get(t, f"Tier {t}"): round(s/total*100,1)
                for t,s in sorted(tier_scores.items(), key=lambda x:-x[1])
            },
            "all_paths": [p.to_dict() for p in ranked],
        }
        cache_set(cache_key, paths_data)
        from_cache = False

    # Paginate
    all_paths = paths_data["all_paths"]
    start     = (page - 1) * page_size
    end       = start + page_size
    page_paths= all_paths[start:end]
    total_p   = len(all_paths)

    return {
        "disease":        node["label"],
        "icd10":          icd10,
        "hcc_code":       node.get("hcc_code"),
        "total_paths":    paths_data["total_paths"],
        "unique_roots":   paths_data["unique_roots"],
        "tier_breakdown": paths_data["tier_breakdown"],
        "pagination": {
            "page":        page,
            "page_size":   page_size,
            "total_results": total_p,
            "total_pages": (total_p + page_size - 1) // page_size,
            "has_next":    end < total_p,
            "from_cache":  from_cache,
        },
        "causal_paths": page_paths,
    }


@app.get("/disease/{icd10}/drugs", tags=["Disease"])
def get_drugs(
    icd10:  str,
    top_n:  int = Query(20, ge=1, le=100),
):
    """
    Get drugs for a disease. TREATS edges ranked first, then ASSOCIATED_WITH.
    """
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")
    conn = get_conn()
    rows = conn.execute("""
        SELECT DISTINCT
            n.primary_id, n.label, n.rxnorm_cui,
            e.relationship_type, e.confidence, e.primary_source,
            json_extract(n.properties,'$.moa') as moa,
            json_extract(n.properties,'$.drugbank_id') as drugbank_id
        FROM edges e
        JOIN nodes n ON n.primary_id=e.source_id AND n.primary_system=e.source_system
        WHERE e.target_id=? AND e.target_system=?
          AND n.entity_type IN ('Drug_clinical','Pharmacologic Substance')
        ORDER BY
            CASE e.relationship_type WHEN 'TREATS' THEN 0 ELSE 1 END,
            e.confidence DESC
        LIMIT ?
    """, (node["primary_id"], node["primary_system"], top_n)).fetchall()

    return {
        "disease": node["label"],
        "icd10":   icd10,
        "count":   len(rows),
        "drugs": [
            {
                "primary_id":   r[0],
                "label":        r[1],
                "rxnorm_cui":   r[2],
                "relationship": r[3],
                "confidence":   r[4],
                "source":       r[5],
                "moa":          r[6],
                "drugbank_id":  r[7],
            }
            for r in rows
        ],
    }


@app.get("/disease/{icd10}/genes", tags=["Disease"])
def get_genes(
    icd10:          str,
    min_confidence: float = Query(0.1, ge=0.0, le=1.0),
    top_n:          int   = Query(20,  ge=1, le=100),
):
    """Gene targets associated with a disease from Open Targets scores."""
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")
    conn = get_conn()
    rows = conn.execute("""
        SELECT DISTINCT n.primary_id,
               json_extract(n.properties,'$.hgnc_symbol') as symbol,
               n.label, e.confidence, e.source_relationship_type
        FROM edges e
        JOIN nodes n ON n.primary_id=e.source_id AND n.primary_system=e.source_system
        WHERE e.target_id=? AND n.primary_system='NCBI_Gene' AND e.confidence>=?
        ORDER BY e.confidence DESC LIMIT ?
    """, (node["primary_id"], min_confidence, top_n)).fetchall()
    return {"disease": node["label"], "icd10": icd10, "count": len(rows),
            "genes": [{"ncbi_gene_id":r[0],"hgnc_symbol":r[1],
                       "label":r[2],"score":r[3],"evidence":r[4]} for r in rows]}


@app.get("/variant/{rsid}", tags=["Genomics"])
def get_variant(rsid: str):
    """Variant node with clinical significance, gnomAD AF, and associated diseases."""
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, label, xrefs, properties, confidence, source
        FROM nodes WHERE primary_id=? AND primary_system='dbSNP_rsID' LIMIT 1
    """, (rsid,)).fetchone()
    if not row:
        raise HTTPException(404, f"Variant not found: {rsid}")
    props = jload(row[3]) or {}
    diseases_raw = conn.execute("""
        SELECT n.label, n.icd10_code, e.relationship_type,
               e.confidence, e.primary_source
        FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=?
          AND e.relationship_type IN ('CAUSES','CONTRIBUTES_TO','INCREASES_RISK_OF')
        ORDER BY e.confidence DESC LIMIT 100
    """, (rsid,)).fetchall()
    # Deduplicate by label
    seen = {}
    for label, icd10, rel, conf, src in diseases_raw:
        if label not in seen or conf > seen[label]["confidence"]:
            seen[label] = {"label":label,"icd10_code":icd10,
                           "relationship":rel,"confidence":conf,
                           "sources":[src] if src else []}
        elif src and src not in seen[label]["sources"]:
            seen[label]["sources"].append(src)
    diseases = sorted(seen.values(), key=lambda x: -x["confidence"])[:20]
    return {
        "rsid":             rsid,
        "label":            row[1],
        "gene_symbol":      props.get("gene_symbol") or props.get("gene"),
        "chromosome":       props.get("chromosome"),
        "position":         props.get("pos_vcf"),
        "clinical_significance": props.get("clinsig"),
        "review_status":    props.get("review_status"),
        "confidence":       row[4],
        "gnomad":           {"af_global": props.get("gnomad_af_global"),
                             "pop_afs": props.get("gnomad_pop_afs")
                             } if props.get("gnomad_af_global") else None,
        "associated_diseases": diseases,
    }


@app.get("/gene/{ncbi_id}", tags=["Genomics"])
def get_gene(ncbi_id: str):
    """Gene node with encoded protein, pathways, and disease associations."""
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, label, xrefs, properties, definition, confidence
        FROM nodes WHERE primary_id=? AND primary_system='NCBI_Gene' LIMIT 1
    """, (ncbi_id,)).fetchone()
    if not row:
        row = conn.execute("""
            SELECT primary_id, label, xrefs, properties, definition, confidence
            FROM nodes WHERE json_extract(properties,'$.hgnc_symbol')=?
              AND primary_system='NCBI_Gene' LIMIT 1
        """, (ncbi_id.upper(),)).fetchone()
    if not row:
        raise HTTPException(404, f"Gene not found: {ncbi_id}")
    props = jload(row[3]) or {}
    xrefs = jload(row[2]) or {}
    protein = conn.execute("""
        SELECT n.primary_id, n.label, n.definition FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type='ENCODES' LIMIT 1
    """, (row[0],)).fetchone()
    pathways = conn.execute("""
        SELECT DISTINCT n.primary_id, n.label FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type='PART_OF' LIMIT 10
    """, (row[0],)).fetchall()
    # Drugs targeting this gene's protein
    drugs = []
    if protein:
        drugs = conn.execute("""
            SELECT DISTINCT nd.label, nd.primary_id,
                   e.source_relationship_type, e.confidence
            FROM edges e
            JOIN nodes nd ON nd.primary_id=e.source_id AND nd.primary_system=e.source_system
            WHERE e.target_id=? AND e.relationship_type='TARGETS'
            ORDER BY e.confidence DESC LIMIT 10
        """, (protein[0],)).fetchall()
    return {
        "ncbi_gene_id": row[0],
        "hgnc_symbol":  props.get("hgnc_symbol"),
        "label":        row[1],
        "chromosome":   props.get("chromosome"),
        "ensembl_id":   xrefs.get("Ensembl"),
        "description":  row[4],
        "protein": {"uniprot_id":protein[0],"label":protein[1],
                    "function":(protein[2] or "")[:200]} if protein else None,
        "pathways": [{"id":r[0],"name":r[1]} for r in pathways],
        "drugs_targeting_protein": [
            {"label":r[0],"primary_id":r[1],
             "action":r[2],"confidence":r[3]} for r in drugs
        ],
    }


@app.get("/protein/{uniprot_id}", tags=["Molecular"])
def get_protein(uniprot_id: str):
    """Protein node with function, gene, interactors, and drugs targeting it."""
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, label, xrefs, properties, definition, confidence
        FROM nodes WHERE primary_id=? AND primary_system='UniProt' LIMIT 1
    """, (uniprot_id,)).fetchone()
    if not row:
        raise HTTPException(404, f"Protein not found: {uniprot_id}")
    props = jload(row[3]) or {}
    xrefs = jload(row[2]) or {}
    gene  = conn.execute("""
        SELECT n.primary_id, n.label, json_extract(n.properties,'$.hgnc_symbol')
        FROM edges e JOIN nodes n ON n.primary_id=e.source_id AND n.primary_system=e.source_system
        WHERE e.target_id=? AND e.relationship_type='ENCODES'
          AND n.primary_system='NCBI_Gene' LIMIT 1
    """, (uniprot_id,)).fetchone()
    drugs = conn.execute("""
        SELECT nd.label, nd.primary_id, e.source_relationship_type, e.confidence
        FROM edges e
        JOIN nodes nd ON nd.primary_id=e.source_id AND nd.primary_system=e.source_system
        WHERE e.target_id=? AND e.relationship_type='TARGETS'
        ORDER BY e.confidence DESC LIMIT 10
    """, (uniprot_id,)).fetchall()
    return {
        "uniprot_id":           uniprot_id,
        "label":                row[1],
        "gene_name":            props.get("gene_name"),
        "function":             row[4],
        "subcellular_location": props.get("subcellular_location"),
        "ncbi_gene_xref":       xrefs.get("NCBI_Gene"),
        "gene": {"ncbi_gene_id":gene[0],"label":gene[1],
                 "hgnc_symbol":gene[2]} if gene else None,
        "drugs_targeting": [
            {"label":r[0],"primary_id":r[1],
             "action":r[2],"confidence":r[3]} for r in drugs
        ],
    }


@app.get("/metabolite/{hmdb_id}", tags=["Molecular"])
def get_metabolite(hmdb_id: str):
    """Metabolite node with formula, biofluid locations, and disease associations."""
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, label, xrefs, properties, definition, confidence
        FROM nodes WHERE primary_id=? AND primary_system='HMDB' LIMIT 1
    """, (hmdb_id,)).fetchone()
    if not row:
        raise HTTPException(404, f"Metabolite not found: {hmdb_id}")
    props = jload(row[3]) or {}
    diseases = conn.execute("""
        SELECT n.label, n.icd10_code, e.relationship_type, e.confidence
        FROM edges e JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.source_system='HMDB'
          AND e.relationship_type='ASSOCIATED_WITH'
          AND n.entity_type='Disease_clinical'
        ORDER BY e.confidence DESC LIMIT 10
    """, (hmdb_id,)).fetchall()
    return {
        "hmdb_id":   hmdb_id,
        "label":     row[1],
        "formula":   props.get("formula"),
        "mass_da":   props.get("mass_da"),
        "biofluids": props.get("biofluids",[]),
        "function":  row[4],
        "associated_diseases": [{"label":r[0],"icd10_code":r[1],
                                  "relationship":r[2],"confidence":r[3]}
                                 for r in diseases],
    }


@app.get("/drug/{rxnorm}", tags=["Drug"])
def get_drug(rxnorm: str):
    """Drug node with MOA, DrugBank ID, ChEMBL targets, and disease indications."""
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, primary_system, label, rxnorm_cui,
               definition, confidence, properties
        FROM nodes WHERE rxnorm_cui=?
          OR (primary_id=? AND primary_system='RxNorm')
        ORDER BY confidence DESC LIMIT 1
    """, (rxnorm, rxnorm)).fetchone()
    if not row:
        raise HTTPException(404, f"Drug not found: {rxnorm}")
    props = jload(row[6]) or {}
    indications = conn.execute("""
        SELECT n.label, n.icd10_code, n.hcc_code,
               e.relationship_type, e.confidence
        FROM edges e JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type IN ('TREATS','ASSOCIATED_WITH')
          AND n.entity_type='Disease_clinical'
        ORDER BY CASE e.relationship_type WHEN 'TREATS' THEN 0 ELSE 1 END,
                 e.confidence DESC LIMIT 20
    """, (row[0],)).fetchall()
    targets = conn.execute("""
        SELECT n.primary_id, n.label, e.source_relationship_type,
               e.confidence, e.primary_source
        FROM edges e JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type='TARGETS'
        ORDER BY e.confidence DESC LIMIT 10
    """, (row[0],)).fetchall()
    return {
        "primary_id":   row[0],
        "label":        row[2],
        "rxnorm_cui":   row[3],
        "drugbank_id":  props.get("drugbank_id"),
        "groups":       props.get("drugbank_groups"),
        "moa":          props.get("moa"),
        "confidence":   row[5],
        "indications": [{"disease":r[0],"icd10_code":r[1],"hcc_code":r[2],
                          "relationship":r[3],"confidence":r[4]}
                         for r in indications],
        "protein_targets": [{"uniprot_id":r[0],"label":r[1],"action":r[2],
                              "confidence":r[3],"source":r[4]}
                             for r in targets],
    }


@app.get("/drug/{rxnorm}/targets", tags=["Drug"])
def get_drug_targets(rxnorm: str):
    """
    Protein targets for a drug with action types from ChEMBL.

    Returns the full chain: Drug→Protein→Gene→Pathway
    """
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, primary_system, label FROM nodes
        WHERE rxnorm_cui=? OR (primary_id=? AND primary_system='RxNorm')
        ORDER BY confidence DESC LIMIT 1
    """, (rxnorm, rxnorm)).fetchone()
    if not row:
        raise HTTPException(404, f"Drug not found: {rxnorm}")

    targets = conn.execute("""
        SELECT np.primary_id, np.label, np.definition,
               e.source_relationship_type, e.confidence, e.primary_source
        FROM edges e
        JOIN nodes np ON np.primary_id=e.target_id AND np.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type='TARGETS'
        ORDER BY e.confidence DESC
    """, (row[0],)).fetchall()

    result = []
    for uniprot, prot_label, prot_func, action, conf, src in targets:
        # Get encoding gene
        gene = conn.execute("""
            SELECT n.primary_id, json_extract(n.properties,'$.hgnc_symbol'), n.label
            FROM edges e JOIN nodes n ON n.primary_id=e.source_id
            WHERE e.target_id=? AND e.relationship_type='ENCODES'
              AND n.primary_system='NCBI_Gene' LIMIT 1
        """, (uniprot,)).fetchone()
        # Get top pathways
        pathways = conn.execute("""
            SELECT n.primary_id, n.label FROM edges e
            JOIN nodes n ON n.primary_id=e.target_id
            WHERE e.source_id=? AND e.relationship_type='PART_OF' LIMIT 5
        """, (uniprot,)).fetchall()
        result.append({
            "uniprot_id":     uniprot,
            "protein_label":  prot_label,
            "function":       (prot_func or "")[:200],
            "action_type":    action,
            "confidence":     conf,
            "pmid_source":    src,
            "encoding_gene":  {"ncbi_id":gene[0],"symbol":gene[1],
                               "label":gene[2]} if gene else None,
            "pathways":       [{"id":r[0],"name":r[1]} for r in pathways],
        })

    return {
        "drug":    row[2],
        "rxnorm":  rxnorm,
        "targets": result,
        "target_count": len(result),
    }


@app.get("/drug/{rxnorm}/mechanism", tags=["Drug"])
def get_drug_mechanism(rxnorm: str):
    """
    Mechanism of action for a drug — MOA text, action types, and target proteins.
    """
    conn = get_conn()
    row  = conn.execute("""
        SELECT primary_id, primary_system, label, properties FROM nodes
        WHERE rxnorm_cui=? OR (primary_id=? AND primary_system='RxNorm')
        ORDER BY confidence DESC LIMIT 1
    """, (rxnorm, rxnorm)).fetchone()
    if not row:
        raise HTTPException(404, f"Drug not found: {rxnorm}")
    props = jload(row[3]) or {}
    targets = conn.execute("""
        SELECT np.label, e.source_relationship_type, e.confidence
        FROM edges e JOIN nodes np ON np.primary_id=e.target_id
        WHERE e.source_id=? AND e.relationship_type='TARGETS'
        ORDER BY e.confidence DESC LIMIT 10
    """, (row[0],)).fetchall()
    action_types = list(set(r[1] for r in targets))
    return {
        "drug":           row[2],
        "rxnorm":         rxnorm,
        "drugbank_id":    props.get("drugbank_id"),
        "moa_text":       props.get("moa"),
        "action_types":   action_types,
        "target_proteins":[{"protein":r[0],"action":r[1],"confidence":r[2]}
                            for r in targets],
        "be_target_ids":  props.get("drugbank_targets"),
    }


@app.post("/patient/trace", tags=["Personalized"])
def patient_trace(request: TraceRequest):
    """
    Personalized causal trace with patient variant/metabolite/comorbidity boosting.

    Results are paginated and cached.
    """
    conn = get_conn()
    if request.disease_system == "ICD-10-CM":
        node = find_disease(request.disease_id)
    else:
        row = conn.execute("""
            SELECT primary_id, primary_system, label FROM nodes
            WHERE primary_id=? AND primary_system=? LIMIT 1
        """, (request.disease_id, request.disease_system)).fetchone()
        node = {"primary_id":row[0],"primary_system":row[1],
                "label":row[2]} if row else None
    if not node:
        raise HTTPException(404, f"Disease not found: {request.disease_id}")

    patient_data = None
    if request.patient_data:
        patient_data = {"variants": request.patient_data.variants,
                        "metabolites": request.patient_data.metabolites,
                        "diseases": request.patient_data.diseases}

    # Cache key includes patient data
    cache_key = (f"trace:{request.disease_id}:{request.disease_system}:"
                 f"{request.max_depth}:{request.min_confidence}:"
                 f"{json.dumps(patient_data, sort_keys=True)}")
    cached = cache_get(cache_key)

    if cached:
        all_paths  = cached["all_paths"]
        total_paths= cached["total_paths"]
        from_cache = True
    else:
        engine = get_engine()
        paths  = engine.trace(
            disease_id     = node["primary_id"],
            disease_system = node["primary_system"],
            max_depth      = request.max_depth,
            min_confidence = request.min_confidence,
            min_path_score = 0.01,
            top_n          = 500,
            patient_data   = patient_data,
        )
        best = {}
        for p in paths:
            k = (p.edges[0].source_id, p.edges[0].source_system)
            if k not in best or p.path_score > best[k].path_score:
                best[k] = p
        ranked = sorted(best.values(), key=lambda p: p.path_score, reverse=True)
        all_paths   = [p.to_dict() for p in ranked]
        total_paths = len(paths)
        cache_set(cache_key, {"all_paths":all_paths,"total_paths":total_paths})
        from_cache  = False

    # Paginate
    start  = (request.page - 1) * request.page_size
    end    = start + request.page_size
    total_r= len(all_paths)

    return {
        "disease":      node["label"],
        "disease_id":   node["primary_id"],
        "personalized": patient_data is not None,
        "patient_data": patient_data,
        "total_paths":  total_paths,
        "unique_roots": len(all_paths),
        "pagination": {
            "page":          request.page,
            "page_size":     request.page_size,
            "total_results": total_r,
            "total_pages":   (total_r + request.page_size - 1) // request.page_size,
            "has_next":      end < total_r,
            "from_cache":    from_cache,
        },
        "top_paths": all_paths[start:end],
    }


@app.get("/search", tags=["Search"])
def search(
    q:      str           = Query(..., min_length=2),
    entity: Optional[str] = Query(None),
    limit:  int           = Query(10, ge=1, le=50),
    page:   int           = Query(1,  ge=1),
):
    """Cross-entity full-text search with pagination."""
    conn         = get_conn()
    page_size    = limit
    offset       = (page - 1) * page_size
    entity_filter= f'AND entity_type = "{entity}"' if entity else ""

    # Total count
    total = conn.execute(f"""
        SELECT COUNT(*) FROM nodes
        WHERE label LIKE ?  {entity_filter}
    """, (f"%{q}%",)).fetchone()[0]

    rows = conn.execute(f"""
        SELECT primary_id, primary_system, label, entity_type,
               tier, icd10_code, hcc_code, rxnorm_cui, confidence
        FROM nodes WHERE label LIKE ? {entity_filter}
        ORDER BY
            CASE WHEN LOWER(label)=LOWER(?) THEN 0
                 WHEN LOWER(label) LIKE LOWER(?)||'%' THEN 1
                 ELSE 2 END,
            confidence DESC
        LIMIT ? OFFSET ?
    """, (f"%{q}%", q, q, page_size, offset)).fetchall()

    return {
        "query": q,
        "pagination": {
            "page": page, "page_size": page_size,
            "total_results": total,
            "total_pages": (total + page_size - 1) // page_size,
            "has_next": offset + page_size < total,
        },
        "results": [{"primary_id":r[0],"primary_system":r[1],"label":r[2],
                     "entity_type":r[3],"tier":r[4],"icd10_code":r[5],
                     "hcc_code":r[6],"rxnorm_cui":r[7],"confidence":r[8]}
                    for r in rows],
    }


@app.delete("/cache", tags=["System"])
def clear_cache():
    """Clear the in-memory causal trace cache. Admin use only."""
    _cache.clear()
    return {"message": "Cache cleared", "entries_cleared": len(_cache)}


@app.get("/cache/stats", tags=["System"])
def cache_stats():
    """View cache statistics."""
    now = time.time()
    active = sum(1 for v in _cache.values() if now < v["expires"])
    return {
        "total_entries": len(_cache),
        "active_entries": active,
        "expired_entries": len(_cache) - active,
        "cache_ttl_seconds": CACHE_TTL,
    }
