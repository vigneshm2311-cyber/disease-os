"""
api/main.py - Disease-OS REST API
Run: uvicorn api.main:app --reload --port 8000
"""
import sys, json, sqlite3
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).parent.parent))

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from core.config import DB_PATH
from core.causal_engine import CausalEngine

app = FastAPI(
    title="Disease-OS API",
    description="Universal causal disease modeling — 3.2M nodes, 10.6M edges",
    version="1.0.0",
)
app.add_middleware(CORSMiddleware, allow_origins=["*"],
                   allow_methods=["*"], allow_headers=["*"])

_conn   = None
_engine = None

def get_conn():
    global _conn
    if _conn is None:
        _conn = sqlite3.connect(str(DB_PATH))
        _conn.execute("PRAGMA cache_size=-256000;")
        _conn.execute("PRAGMA temp_store=MEMORY;")
    return _conn

def get_engine():
    global _engine
    if _engine is None:
        _engine = CausalEngine(str(DB_PATH))
    return _engine

def jload(s):
    try: return json.loads(s) if s else None
    except: return s

def find_disease(icd10):
    conn = get_conn()
    row = conn.execute("""
        SELECT primary_id, primary_system, label, icd10_code,
               hcc_code, icd11_code, xrefs, synonyms, definition, confidence
        FROM nodes WHERE icd10_code=? ORDER BY confidence DESC LIMIT 1
    """, (icd10,)).fetchone()
    if not row:
        row = conn.execute("""
            SELECT primary_id, primary_system, label, icd10_code,
                   hcc_code, icd11_code, xrefs, synonyms, definition, confidence
            FROM nodes WHERE primary_id=? AND primary_system="ICD-10-CM" LIMIT 1
        """, (icd10,)).fetchone()
    if not row: return None
    cols = ["primary_id","primary_system","label","icd10_code",
            "hcc_code","icd11_code","xrefs","synonyms","definition","confidence"]
    d = dict(zip(cols, row))
    for f in ["xrefs","synonyms"]: d[f] = jload(d[f])
    return d


class PatientData(BaseModel):
    variants:    list[str] = []
    metabolites: list[str] = []
    diseases:    list[str] = []

class TraceRequest(BaseModel):
    disease_id:     str
    disease_system: str = "ICD-10-CM"
    max_depth:      int = 4
    min_confidence: float = 0.50
    top_n:          int = 20
    patient_data:   Optional[PatientData] = None


@app.get("/health", tags=["System"])
def health():
    conn = get_conn()
    nodes = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    edges = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    tiers = dict(conn.execute(
        "SELECT tier, COUNT(*) FROM nodes GROUP BY tier ORDER BY tier"
    ).fetchall())
    rels = dict(conn.execute("""
        SELECT relationship_type, COUNT(*) FROM edges
        GROUP BY relationship_type ORDER BY COUNT(*) DESC
    """).fetchall())
    return {
        "status": "healthy", "version": "1.0.0",
        "graph": {"total_nodes": nodes, "total_edges": edges,
                  "nodes_by_tier": tiers, "edges_by_relationship": rels},
        "sources": ["UMLS 2025AB","Reactome v88","ClinVar 2024",
                    "STRING v12","GWAS Catalog","HMDB v5.0",
                    "EFO v3.91","Open Targets 26.06",
                    "NCBI Gene","UniProt Swiss-Prot","HPO 2024"],
    }


@app.get("/disease/{icd10}", tags=["Disease"])
def get_disease(icd10: str):
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")
    return node


@app.get("/disease/{icd10}/causes", tags=["Disease","Causal"])
def get_causes(
    icd10: str,
    max_depth: int = Query(4, ge=1, le=6),
    min_confidence: float = Query(0.5, ge=0.0, le=1.0),
    top_n: int = Query(20, ge=1, le=100),
):
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")
    engine = get_engine()
    paths = engine.trace(
        disease_id=node["primary_id"],
        disease_system=node["primary_system"],
        max_depth=max_depth, min_confidence=min_confidence,
        min_path_score=0.01, top_n=500,
    )
    best = {}
    for p in paths:
        k = (p.edges[0].source_id, p.edges[0].source_system)
        if k not in best or p.path_score > best[k].path_score:
            best[k] = p
    ranked = sorted(best.values(), key=lambda p: p.path_score, reverse=True)[:top_n]
    tier_names = {1:"Molecular",2:"Networks",3:"Cellular",4:"Tissue",
                  5:"Systemic",6:"Phenotype",7:"Disease",8:"Behavior",9:"Social",10:"Healthcare"}
    tier_scores = {}
    for p in paths:
        t = p.edges[0].tier_from
        tier_scores[t] = tier_scores.get(t,0) + p.path_score
    total = sum(tier_scores.values()) or 1
    return {
        "disease": node["label"], "icd10": icd10, "hcc_code": node.get("hcc_code"),
        "total_paths": len(paths), "unique_roots": len(best),
        "tier_breakdown": {tier_names.get(t,f"Tier {t}"): round(s/total*100,1)
                           for t,s in sorted(tier_scores.items(), key=lambda x:-x[1])},
        "causal_paths": [p.to_dict() for p in ranked],
    }


@app.get("/disease/{icd10}/drugs", tags=["Disease"])
def get_drugs(icd10: str, top_n: int = Query(20, ge=1, le=100)):
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")
    conn = get_conn()
    rows = conn.execute("""
        SELECT DISTINCT n.primary_id, n.label, n.rxnorm_cui,
               e.relationship_type, e.confidence, e.primary_source
        FROM edges e
        JOIN nodes n ON n.primary_id=e.source_id AND n.primary_system=e.source_system
        WHERE e.target_id=? AND e.target_system=?
          AND n.entity_type IN ("Drug_clinical","Pharmacologic Substance")
        ORDER BY e.confidence DESC LIMIT ?
    """, (node["primary_id"], node["primary_system"], top_n)).fetchall()
    return {"disease": node["label"], "icd10": icd10, "count": len(rows),
            "drugs": [{"primary_id":r[0],"label":r[1],"rxnorm_cui":r[2],
                       "relationship":r[3],"confidence":r[4],"source":r[5]} for r in rows]}


@app.get("/disease/{icd10}/genes", tags=["Disease"])
def get_genes(icd10: str,
              min_confidence: float = Query(0.1, ge=0.0, le=1.0),
              top_n: int = Query(20, ge=1, le=100)):
    node = find_disease(icd10)
    if not node:
        raise HTTPException(404, f"Disease not found: {icd10}")
    conn = get_conn()
    rows = conn.execute("""
        SELECT DISTINCT n.primary_id,
               json_extract(n.properties,"$.hgnc_symbol") as symbol,
               n.label, e.confidence, e.source_relationship_type
        FROM edges e
        JOIN nodes n ON n.primary_id=e.source_id AND n.primary_system=e.source_system
        WHERE e.target_id=? AND n.primary_system="NCBI_Gene" AND e.confidence>=?
        ORDER BY e.confidence DESC LIMIT ?
    """, (node["primary_id"], min_confidence, top_n)).fetchall()
    return {"disease": node["label"], "icd10": icd10, "count": len(rows),
            "genes": [{"ncbi_gene_id":r[0],"hgnc_symbol":r[1],
                       "label":r[2],"score":r[3],"evidence":r[4]} for r in rows]}


@app.get("/variant/{rsid}", tags=["Genomics"])
def get_variant(rsid: str):
    conn = get_conn()
    row = conn.execute("""
        SELECT primary_id, label, xrefs, properties, confidence, source
        FROM nodes WHERE primary_id=? AND primary_system="dbSNP_rsID" LIMIT 1
    """, (rsid,)).fetchone()
    if not row:
        raise HTTPException(404, f"Variant not found: {rsid}")
    props = jload(row[3]) or {}
    diseases = conn.execute("""
        SELECT n.label, n.icd10_code, e.relationship_type, e.confidence, e.primary_source
        FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=?
          AND e.relationship_type IN ("CAUSES","CONTRIBUTES_TO","INCREASES_RISK_OF")
        ORDER BY e.confidence DESC LIMIT 20
    """, (rsid,)).fetchall()
    return {
        "rsid": rsid, "label": row[1],
        "gene_symbol": props.get("gene_symbol") or props.get("gene"),
        "chromosome": props.get("chromosome"),
        "position": props.get("pos_vcf"),
        "clinical_significance": props.get("clinsig"),
        "review_status": props.get("review_status"),
        "confidence": row[4],
        "gnomad": {"af_global": props.get("gnomad_af_global"),
                   "pop_afs": props.get("gnomad_pop_afs")} if props.get("gnomad_af_global") else None,
        "associated_diseases": [{"label":r[0],"icd10_code":r[1],
                                  "relationship":r[2],"confidence":r[3],"source":r[4]}
                                 for r in diseases],
    }


@app.get("/gene/{ncbi_id}", tags=["Genomics"])
def get_gene(ncbi_id: str):
    conn = get_conn()
    row = conn.execute("""
        SELECT primary_id, label, xrefs, properties, definition, confidence
        FROM nodes WHERE primary_id=? AND primary_system="NCBI_Gene"
        LIMIT 1
    """, (ncbi_id,)).fetchone()
    if not row:
        row = conn.execute("""
            SELECT primary_id, label, xrefs, properties, definition, confidence
            FROM nodes WHERE json_extract(properties,"$.hgnc_symbol")=?
              AND primary_system="NCBI_Gene" LIMIT 1
        """, (ncbi_id.upper(),)).fetchone()
    if not row:
        raise HTTPException(404, f"Gene not found: {ncbi_id}")
    props = jload(row[3]) or {}
    xrefs = jload(row[2]) or {}
    protein = conn.execute("""
        SELECT n.primary_id, n.label, n.definition FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type="ENCODES" LIMIT 1
    """, (row[0],)).fetchone()
    pathways = conn.execute("""
        SELECT DISTINCT n.primary_id, n.label FROM edges e
        JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type="PART_OF" LIMIT 10
    """, (row[0],)).fetchall()
    return {
        "ncbi_gene_id": row[0], "label": row[1],
        "hgnc_symbol": props.get("hgnc_symbol"),
        "chromosome": props.get("chromosome"),
        "ensembl_id": xrefs.get("Ensembl"),
        "description": row[4],
        "protein": {"uniprot_id":protein[0],"label":protein[1],
                    "function":(protein[2] or "")[:200]} if protein else None,
        "pathways": [{"id":r[0],"name":r[1]} for r in pathways],
    }


@app.get("/protein/{uniprot_id}", tags=["Molecular"])
def get_protein(uniprot_id: str):
    conn = get_conn()
    row = conn.execute("""
        SELECT primary_id, label, xrefs, properties, definition, confidence
        FROM nodes WHERE primary_id=? AND primary_system="UniProt" LIMIT 1
    """, (uniprot_id,)).fetchone()
    if not row:
        raise HTTPException(404, f"Protein not found: {uniprot_id}")
    props = jload(row[3]) or {}
    xrefs = jload(row[2]) or {}
    gene = conn.execute("""
        SELECT n.primary_id, n.label, json_extract(n.properties,"$.hgnc_symbol")
        FROM edges e JOIN nodes n ON n.primary_id=e.source_id AND n.primary_system=e.source_system
        WHERE e.target_id=? AND e.relationship_type="ENCODES"
          AND n.primary_system="NCBI_Gene" LIMIT 1
    """, (uniprot_id,)).fetchone()
    return {
        "uniprot_id": uniprot_id, "label": row[1],
        "gene_name": props.get("gene_name"),
        "function": row[4],
        "subcellular_location": props.get("subcellular_location"),
        "active_sites": props.get("active_sites",[]),
        "ncbi_gene_xref": xrefs.get("NCBI_Gene"),
        "ensembl_xref": xrefs.get("Ensembl"),
        "gene": {"ncbi_gene_id":gene[0],"label":gene[1],"hgnc_symbol":gene[2]} if gene else None,
    }


@app.get("/metabolite/{hmdb_id}", tags=["Molecular"])
def get_metabolite(hmdb_id: str):
    conn = get_conn()
    row = conn.execute("""
        SELECT primary_id, label, xrefs, properties, definition, confidence
        FROM nodes WHERE primary_id=? AND primary_system="HMDB" LIMIT 1
    """, (hmdb_id,)).fetchone()
    if not row:
        raise HTTPException(404, f"Metabolite not found: {hmdb_id}")
    props = jload(row[3]) or {}
    diseases = conn.execute("""
        SELECT n.label, n.icd10_code, e.relationship_type, e.confidence
        FROM edges e JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.source_system="HMDB"
          AND e.relationship_type="ASSOCIATED_WITH"
          AND n.entity_type="Disease_clinical"
        ORDER BY e.confidence DESC LIMIT 10
    """, (hmdb_id,)).fetchall()
    return {
        "hmdb_id": hmdb_id, "label": row[1],
        "formula": props.get("formula"),
        "mass_da": props.get("mass_da"),
        "biofluids": props.get("biofluids",[]),
        "function": row[4],
        "associated_diseases": [{"label":r[0],"icd10_code":r[1],
                                  "relationship":r[2],"confidence":r[3]} for r in diseases],
    }


@app.get("/drug/{rxnorm}", tags=["Drug"])
def get_drug(rxnorm: str):
    conn = get_conn()
    row = conn.execute("""
        SELECT primary_id, primary_system, label, rxnorm_cui, definition, confidence
        FROM nodes WHERE rxnorm_cui=? OR (primary_id=? AND primary_system="RxNorm")
        ORDER BY confidence DESC LIMIT 1
    """, (rxnorm, rxnorm)).fetchone()
    if not row:
        raise HTTPException(404, f"Drug not found: {rxnorm}")
    indications = conn.execute("""
        SELECT n.label, n.icd10_code, n.hcc_code, e.relationship_type, e.confidence
        FROM edges e JOIN nodes n ON n.primary_id=e.target_id AND n.primary_system=e.target_system
        WHERE e.source_id=? AND e.relationship_type IN ("TREATS","ASSOCIATED_WITH")
          AND n.entity_type="Disease_clinical"
        ORDER BY e.confidence DESC LIMIT 20
    """, (row[0],)).fetchall()
    return {
        "primary_id": row[0], "label": row[2], "rxnorm_cui": row[3],
        "definition": row[4], "confidence": row[5],
        "indications": [{"disease":r[0],"icd10_code":r[1],"hcc_code":r[2],
                          "relationship":r[3],"confidence":r[4]} for r in indications],
    }


@app.post("/patient/trace", tags=["Personalized"])
def patient_trace(request: TraceRequest):
    conn = get_conn()
    if request.disease_system == "ICD-10-CM":
        node = find_disease(request.disease_id)
    else:
        row = conn.execute("""
            SELECT primary_id, primary_system, label FROM nodes
            WHERE primary_id=? AND primary_system=? LIMIT 1
        """, (request.disease_id, request.disease_system)).fetchone()
        node = {"primary_id":row[0],"primary_system":row[1],"label":row[2]} if row else None
    if not node:
        raise HTTPException(404, f"Disease not found: {request.disease_id}")
    patient_data = None
    if request.patient_data:
        patient_data = {"variants": request.patient_data.variants,
                        "metabolites": request.patient_data.metabolites,
                        "diseases": request.patient_data.diseases}
    engine = get_engine()
    paths = engine.trace(
        disease_id=node["primary_id"], disease_system=node["primary_system"],
        max_depth=request.max_depth, min_confidence=request.min_confidence,
        min_path_score=0.01, top_n=500, patient_data=patient_data,
    )
    best = {}
    for p in paths:
        k = (p.edges[0].source_id, p.edges[0].source_system)
        if k not in best or p.path_score > best[k].path_score:
            best[k] = p
    ranked = sorted(best.values(), key=lambda p: p.path_score, reverse=True)[:request.top_n]
    return {
        "disease": node["label"], "personalized": patient_data is not None,
        "total_paths": len(paths), "unique_roots": len(best),
        "top_paths": [p.to_dict() for p in ranked],
    }


@app.get("/search", tags=["Search"])
def search(
    q: str = Query(..., min_length=2),
    entity: Optional[str] = Query(None),
    limit: int = Query(10, ge=1, le=50),
):
    conn = get_conn()
    entity_filter = f'AND entity_type = "{entity}"' if entity else ""
    rows = conn.execute(f"""
        SELECT primary_id, primary_system, label, entity_type,
               tier, icd10_code, hcc_code, rxnorm_cui, confidence
        FROM nodes
        WHERE label LIKE ?
        {entity_filter}
        ORDER BY
            CASE WHEN LOWER(label)=LOWER(?) THEN 0
                 WHEN LOWER(label) LIKE LOWER(?)||"%" THEN 1
                 ELSE 2 END,
            confidence DESC
        LIMIT ?
    """, (f"%{q}%", q, q, limit)).fetchall()
    return {
        "query": q, "count": len(rows),
        "results": [{"primary_id":r[0],"primary_system":r[1],"label":r[2],
                     "entity_type":r[3],"tier":r[4],"icd10_code":r[5],
                     "hcc_code":r[6],"rxnorm_cui":r[7],"confidence":r[8]} for r in rows],
    }
