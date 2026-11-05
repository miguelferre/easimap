"""Conjunto actual de estructuras de una familia, con los filtros del plan.

Filtros, en este orden (cada paso deja su cuenta en el informe):
1. Entidades de polímero con la anotación de la familia en RCSB (InterPro IPR001360 para GH1).
2. Estructura experimental (se excluyen los modelos computados).
3. Resolución de 3 Å o mejor. Sin resolución (RMN) no pasa.
4. Sin mutante catalítico. Lo decide `build`: la entidad tiene mutación anotada y el residuo de una columna catalítica no es Glu.
5. Agrupado al 90 % de identidad de secuencia. Cada grupo aporta la entidad de mejor resolución.
6. Cadena del dominio: la que PDBe mapea a la familia, no siempre la A.

No hay scraping de CAZy. Todo sale de RCSB y PDBe.
"""
from __future__ import annotations

import datetime
import json
import re
from pathlib import Path

import requests

from . import fetch

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"

ENTITY_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity { pdbx_mutation pdbx_description }
    rcsb_polymer_entity_container_identifiers { entry_id entity_id auth_asym_ids asym_ids }
    rcsb_polymer_entity_group_membership { aggregation_method group_id similarity_cutoff }
    rcsb_entity_source_organism { ncbi_scientific_name taxonomy_lineage { id name } }
    entry { rcsb_entry_info { resolution_combined experimental_method } exptl { method } }
  }
}
"""


def annotated_entities(accession: str, kind: str = "Pfam") -> list[str]:
    body = {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id", "operator": "exact_match", "value": accession}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.type", "operator": "exact_match", "value": kind}}]},
        "return_type": "polymer_entity", "request_options": {"return_all_hits": True}}
    r = requests.post(SEARCH, json=body, timeout=60)
    r.raise_for_status()
    return [x["identifier"] for x in r.json()["result_set"]]


def entity_metadata(ids: list[str], batch: int = 40) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(ids), batch):
        r = requests.post(GRAPHQL, json={"query": ENTITY_QUERY, "variables": {"ids": ids[i:i + batch]}}, timeout=120)
        r.raise_for_status()
        out += r.json()["data"]["polymer_entities"]
    return out


def lineage_domain(entity: dict) -> str:
    """Eukaryota, Bacteria, Archaea o Unknown, según el linaje taxonómico de la entidad.

    Una quimera de dos dominios (4GXP: Thermotoga y Trichoderma) da Unknown: no tiene un solo dominio.
    """
    found = set()
    for org in entity.get("rcsb_entity_source_organism") or []:
        names = {t["name"] for t in (org.get("taxonomy_lineage") or [])}
        found |= names & {"Eukaryota", "Bacteria", "Archaea"}
    return found.pop() if len(found) == 1 else "Unknown"


def _mutated_positions(entity: dict) -> list[int]:
    txt = (entity["rcsb_polymer_entity"] or {}).get("pdbx_mutation") or ""
    return [int(m) for m in re.findall(r"[A-Z](\d+)[A-Z]", txt)]


def select(accession: str, kind: str = "InterPro", max_res: float = 3.0, identity: int = 90) -> dict:
    """Aplica los filtros y devuelve la lista final y las cuentas de cada paso.

    La anotación de RCSB por Pfam (PF00232) deja fuera casi todo GH1 con el modelo vigente: solo 22 de las 66
    estructuras de 2020. La entrada de InterPro de la familia (IPR001360) recoge 60 de las 66.
    """
    counts: dict[str, int] = {}
    ids = sorted(annotated_entities(accession, kind))     # RCSB no garantiza el orden: se fija para que la salida sea estable
    counts[f"entidades con anotación {kind} {accession}"] = len(ids)
    ents = entity_metadata(ids)
    rows = []
    for e in ents:
        info = e["entry"]["rcsb_entry_info"]
        method = (info.get("experimental_method") or "")
        res = (info.get("resolution_combined") or [None])[0]
        ident = e["rcsb_polymer_entity_container_identifiers"]
        grp = next((g for g in (e["rcsb_polymer_entity_group_membership"] or [])
                    if g["aggregation_method"] == "sequence_identity" and g["similarity_cutoff"] == float(identity)), None)
        rows.append({
            "entity": e["rcsb_id"], "pdb": ident["entry_id"], "entity_id": ident["entity_id"],
            "method": method, "resolution": res, "auth_asym_ids": ident["auth_asym_ids"],
            "group": grp["group_id"] if grp else None, "domain": lineage_domain(e),
            "organism": ((e.get("rcsb_entity_source_organism") or [{}])[0] or {}).get("ncbi_scientific_name"),
            "mutation": (e["rcsb_polymer_entity"] or {}).get("pdbx_mutation") or "",
            "description": (e["rcsb_polymer_entity"] or {}).get("pdbx_description") or "",
        })
    rows.sort(key=lambda r: r["entity"])                  # orden fijo: la salida no depende de la respuesta de RCSB
    rows = [r for r in rows if r["method"] and "MODEL" not in r["method"].upper()]
    counts["experimentales"] = len(rows)
    rows = [r for r in rows if r["resolution"] is not None and r["resolution"] <= max_res]
    counts[f"resolución de {max_res} Å o mejor"] = len(rows)
    return {"rows": rows, "counts": counts, "accession": accession, "max_resolution": max_res, "identity": identity,
            "date": datetime.date.today().isoformat()}


def cluster_best(rows: list[dict]) -> list[dict]:
    """Una entidad por grupo de identidad: la de mejor resolución. Sin grupo, la entidad queda sola."""
    best: dict[str, dict] = {}
    for r in sorted(rows, key=lambda r: (r["resolution"], r["entity"])):
        best.setdefault(r["group"] or r["entity"], r)
    return list(best.values())


def domain_chain(code: str, pfam: str, interpro: list[str], allowed: list[str]) -> tuple[str, str] | None:
    """Cadena (auth) donde PDBe mapea la familia. Devuelve (cadena, origen del dato).

    Primero /mappings/pfam. Si PDBe no mapea Pfam, /mappings/interpro. Solo vale una cadena de la entidad.
    """
    for kind, accs in (("pfam", [pfam]), ("interpro", interpro)):
        try:
            block = fetch.pdbe_mapping(code, kind).get(code.lower(), {}).get("Pfam" if kind == "pfam" else "InterPro", {})
        except Exception:
            block = {}
        for acc in accs:
            data = block.get(acc)
            if data:
                chains = sorted({m["chain_id"] for m in data["mappings"]} & set(allowed))
                if chains:
                    return chains[0], f"PDBe {kind}"
    return (sorted(allowed)[0], "primera cadena de la entidad") if allowed else None
