"""Construcción de web/data/<familia>.json con los métodos nuevos.

Respecto a la réplica de 2020:
* Cadena del dominio con PDBe, no siempre la A. Coordenadas del mmCIF, con la numeración de autor.
* Solo el primer altloc, sin hidrógenos. MSE cuenta como M.
* Una entrada por proteína, cadena y posición HMM. Solo residuos con posición.
* Dos superposiciones: `catalytic` (residuos catalíticos) y `core` (núcleo de la familia).
* RMSF en Å por posición, que sustituye al SCV (D-004). Mapa de horizonte.
"""
from __future__ import annotations

import collections
import datetime
import json
import math
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import numpy as np

from . import cache, fetch, select_set
from .geometry import kabsch, spherical, substrate_frame, to_frame
from .hmmalign import assign_positions_local, load_hmm
from .pfam_archive import archived_hmm
from .structure import frame_atoms

ROOT = Path(__file__).resolve().parents[2]
BACKBONE4 = ("N", "CA", "C", "O")
# Átomos que se tratan como esqueleto al medir la cadena lateral (mismo criterio que la réplica, sin hidrógenos).
BACKBONE_ALL = {"N", "CA", "C", "O", "OXT", "OT1", "OT2"}
ONE = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
       "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TYR": "Y", "TRP": "W",
       "MSE": "M"}
HORIZON_GRID = 26
HORIZON_CONE = 14.0     # grados
HORIZON_FAR = 40.0      # Å, si ningún residuo cae en el cono
MIN_POSITION_N = 5      # estructuras mínimas para publicar el RMSF de una posición


@dataclass
class Prot:
    pdb: str
    chain: str
    domain: str
    organism: str
    resolution: float
    resid: list[int]
    resname: list[str]
    pos: list[int]                    # posición HMM por residuo; 0 si no tiene
    backbone: dict[int, np.ndarray]   # posición HMM -> matriz 4x3 (N, CA, C, O)
    ca: dict[int, np.ndarray]         # posición HMM -> CA
    point: dict[int, np.ndarray]      # posición HMM -> centro de cadena lateral más CA
    skip: str = ""


def read_chain(code: str, chain: str):
    """Residuos polímero de una cadena (nombre de autor). Primer altloc, sin hidrógenos."""
    st = gemmi.read_structure(str(fetch.cif_file(code)))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    ch = next((c for c in st[0] if c.name == chain), None)
    if ch is None:
        return []
    out = []
    for res in ch:
        if res.entity_type != gemmi.EntityType.Polymer:
            continue
        atoms = {a.name: np.array(a.pos.tolist()) for a in res}
        if all(n in atoms for n in BACKBONE4):
            out.append((res.seqid.num, res.seqid.icode.strip(), res.name, atoms))
    return out


def prepare(code: str, chain: str, hmm_path: Path, meta: dict, min_positions: int = 250) -> Prot:
    """Alinea una cadena al HMM y guarda lo que hacen falta las dos superposiciones. Con caché."""
    params = {"pdb": code, "chain": chain, "hmm": cache.file_hash(hmm_path), "altloc": "first", "hydrogens": "removed",
              "mse": "M", "unknown": "X", "alignment": "local-unihit", "source": "mmcif"}

    def compute() -> Prot:
        res = read_chain(code, chain)
        p = Prot(code, chain, meta.get("domain", "Unknown"), meta.get("organism") or "", meta.get("resolution") or 0.0,
                 [], [], [], {}, {}, {})
        if not res:
            p.skip = "la cadena no tiene residuos de polímero con esqueleto completo"
            return p
        seq = "".join(ONE.get(r[2], "X") for r in res)
        pos = assign_positions_local(load_hmm(hmm_path), seq)
        for (num, icode, name, atoms), q in zip(res, pos):
            if q <= 0:
                continue
            p.resid.append(num)
            p.resname.append("MET" if name == "MSE" else name)
            p.pos.append(q)
            p.backbone[q] = np.array([atoms[n] for n in BACKBONE4])
            p.ca[q] = atoms["CA"]
            side = [v for n, v in atoms.items() if n not in BACKBONE_ALL]
            p.point[q] = np.mean(side + [atoms["CA"]], axis=0)
        return p

    prot = cache.memo("prot", params, compute)
    if meta:   # los metadatos no cambian la medida: no entran en la clave de caché
        prot.domain = meta.get("domain", prot.domain)
        prot.organism = meta.get("organism") or prot.organism
        prot.resolution = meta.get("resolution") or prot.resolution
    if not prot.skip and len(prot.pos) < min_positions:    # el umbral es de la familia, no entra en la clave de caché
        prot.skip = f"solo {len(prot.pos)} residuos con posición HMM"
    return prot


def reference_substrate(cfg: dict):
    """Átomos que fijan el marco en la estructura de referencia: origen, eje y y plano x (ver `[substrate]` de la familia)."""
    return frame_atoms(fetch.pdb_file(cfg["reference_pdb"]), cfg["substrate"])


def catalytic_fit(prots: list[Prot], ref: Prot, cols: list[int]):
    """Rotación y traslación de cada proteína sobre la referencia, con N, CA, C, O de las columnas catalíticas."""
    target = np.vstack([ref.backbone[c] for c in cols])
    fits = {}
    for p in prots:
        if not all(c in p.backbone for c in cols):
            p.skip = p.skip or "falta un residuo catalítico en las columnas del HMM"
            continue
        fits[(p.pdb, p.chain)] = kabsch(np.vstack([p.backbone[c] for c in cols]), target)
    return fits


def core_fit(prots: list[Prot], base: dict, rmsf_in: float = 3.0, rmsf_out: float = 2.0, rounds: int = 8):
    """Superposición sobre el núcleo de la familia: posiciones comunes cuyo CA apenas se mueve.

    Parte de las posiciones presentes en el 90 % de las estructuras con RMSF menor de 3 Å tras el ajuste
    catalítico. Ajusta todas a la media del núcleo y repite, con umbral de 2 Å. Devuelve (transformaciones
    acumuladas, posiciones del núcleo).
    """
    keys = [(p.pdb, p.chain) for p in prots if (p.pdb, p.chain) in base]
    byk = {(p.pdb, p.chain): p for p in prots}
    xyz = {k: {q: base[k][0] @ v + base[k][1] for q, v in byk[k].ca.items()} for k in keys}
    cnt = collections.Counter(q for k in keys for q in xyz[k])
    full = sorted(q for q, c in cnt.items() if c >= 0.9 * len(keys))

    def rmsf(coords):
        out = {}
        for q in cnt:
            X = np.array([coords[k][q] for k in keys if q in coords[k]])
            if len(X) >= 5:
                out[q] = float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))
        return out

    acc = {k: (np.eye(3), np.zeros(3)) for k in keys}
    cur = rmsf(xyz)
    core = [q for q in full if cur.get(q, 99) < rmsf_in]
    for _ in range(rounds):
        mean = {q: np.mean([xyz[k][q] for k in keys if q in xyz[k]], axis=0) for q in core}
        for k in keys:
            use = [q for q in core if q in xyz[k]]
            if len(use) < 20:
                continue
            r, t = kabsch(np.array([xyz[k][q] for q in use]), np.array([mean[q] for q in use]))
            xyz[k] = {q: r @ v + t for q, v in xyz[k].items()}
            acc[k] = (r @ acc[k][0], r @ acc[k][1] + t)
        cur = rmsf(xyz)
        core = [q for q in full if cur.get(q, 99) < rmsf_out]
    return acc, core


def measure(prots: list[Prot], transform: dict, ligand: tuple, ligand_xf: tuple):
    """Coordenadas esféricas de cada punto medido.

    `transform` lleva de las coordenadas originales al marco común. `ligand_xf` es la transformación
    de la referencia (el sustrato está en sus coordenadas), con la que se construyen los ejes.
    """
    c1, o, c4 = ligand
    r_ref, t_ref = ligand_xf
    origin, axes = substrate_frame(r_ref @ c1 + t_ref, r_ref @ o + t_ref, r_ref @ c4 + t_ref)
    rows = []
    for p in prots:
        k = (p.pdb, p.chain)
        if k not in transform:
            continue
        r, t = transform[k]
        for i, q in enumerate(p.pos):
            pt = to_frame(r @ p.point[q] + t, origin, axes)[0]
            d, th, ph = spherical(pt)
            rows.append((p, i, q, d, th, ph, pt))
    return rows


def columns(rows, domain_of: dict) -> dict:
    out = {k: [] for k in ("pdb", "chain", "auth_seq_id", "hmm_position", "residue", "domain", "theta", "phi", "distance", "x", "y", "z")}
    for p, i, q, d, th, ph, pt in rows:
        out["pdb"].append(p.pdb); out["chain"].append(p.chain); out["auth_seq_id"].append(p.resid[i])
        out["hmm_position"].append(q); out["residue"].append(p.resname[i]); out["domain"].append(p.domain)
        out["theta"].append(round(th, 4)); out["phi"].append(round(ph, 4)); out["distance"].append(round(d, 3))
        out["x"].append(round(float(pt[0]), 3)); out["y"].append(round(float(pt[1]), 3)); out["z"].append(round(float(pt[2]), 3))
    return out


def position_stats(rows) -> list[dict]:
    by = collections.defaultdict(list)
    names = collections.defaultdict(collections.Counter)
    for p, i, q, d, th, ph, pt in rows:
        by[q].append(pt)
        names[q][p.resname[i]] += 1
    out = []
    for q in sorted(by):
        X = np.array(by[q])
        if len(X) < MIN_POSITION_N:
            continue
        m = X.mean(0)
        out.append({"hmm_position": int(q), "n": int(len(X)), "rmsf": round(float(np.sqrt(((X - m) ** 2).sum(1).mean())), 2),
                    "distance": round(float(np.linalg.norm(m)), 2), "top_residue": names[q].most_common(1)[0][0]})
    return out


def horizon(rows) -> dict:
    """Mediana entre estructuras de la distancia al residuo más cercano en cada dirección (proyección de área igual).

    Malla de GRID x GRID por hemisferio, en el orden `j * GRID + i`. Fuera del disco vale 0.
    """
    by = collections.defaultdict(list)
    for p, i, q, d, th, ph, pt in rows:
        by[(p.pdb, p.chain)].append(pt)
    mols = []
    for X in by.values():
        X = np.unique(np.round(np.array(X), 4), axis=0)
        D = np.linalg.norm(X, axis=1)
        mols.append((X / D[:, None], D))
    cosc = math.cos(math.radians(HORIZON_CONE))
    G = HORIZON_GRID
    out = {}
    for name, south in (("north", False), ("south", True)):
        cells = [0.0] * (G * G)
        for i in range(G):
            for j in range(G):
                u, v = (i + .5) / G * 2 - 1, (j + .5) / G * 2 - 1
                rr = math.hypot(u, v)
                if rr > 1:
                    continue
                t = 2 * math.asin(rr / math.sqrt(2))
                ph = math.atan2(u * (1 if south else -1), -v)
                th = math.pi - t if south else t
                dv = np.array([math.sin(th) * math.cos(ph), math.cos(th), math.sin(th) * math.sin(ph)])
                vals = []
                for U, D in mols:
                    near = D[(U @ dv) > cosc]
                    vals.append(near.min() if len(near) else HORIZON_FAR)
                cells[j * G + i] = round(float(np.median(vals)), 2)
        out[name] = cells
    return {"grid": G, "north": out["north"], "south": out["south"]}


def build(family: str = "GH1", hmm_path: Path | None = None, out: Path | None = None, report_path: Path | None = None,
          limit: int | None = None) -> dict:
    cfg = tomllib.loads((ROOT / "familias" / f"{family}.toml").read_text(encoding="utf-8"))
    rcfg = cfg["replica2020"]
    hmm_path = hmm_path or archived_hmm(cfg["pfam"], cfg["hmm_release"])
    ref_code = cfg["reference_pdb"]
    ref_chain = cfg.get("reference_chain", "A")
    min_pos = cfg.get("min_positions", 250)
    expected = cfg["catalytic_residues"]
    log: dict = {"family": family, "date": datetime.date.today().isoformat(), "hmm": hmm_path.name, "hmm_sha": cache.file_hash(hmm_path),
                 "skipped": {}, "counts": {}}

    # 1. Candidatos y filtros de metadatos
    sel = select_set.select(cfg["interpro"], "InterPro", max_res=3.0, identity=90)
    log["counts"].update(sel["counts"])
    rows = sel["rows"][:limit] if limit else sel["rows"]

    # 2. Cadena del dominio y alineamiento de cada candidato
    ref = prepare(ref_code, ref_chain, hmm_path, {"domain": cfg.get("reference_domain", "Unknown")}, min_pos)
    if rcfg.get("numbering_2020", True):
        ref_cols = list(rcfg["reference_positions_2020"])
    else:
        # El HMM disponible no reproduce la numeración de 2020 (ver el hallazgo de la familia): las columnas de
        # referencia son las que el modelo da a los residuos de referencia.
        ref_cols = []
        for rid in rcfg["reference_resids"]:
            num = int(str(rid).rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
            ref_cols.append(next(q for q, r in zip(ref.pos, ref.resid) if r == num))
    # Comprobación: cada columna cae sobre el residuo de referencia. Un resid con código de inserción (GH16: 4ATF)
    # puede salir dos veces en la cadena; basta que la columna sea una de las suyas.
    for col, rid in zip(ref_cols, rcfg["reference_resids"]):
        num = int(str(rid).rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
        mine = [q for q, r in zip(ref.pos, ref.resid) if r == num]
        assert col in mine, f"la columna {col} no cae sobre el resid {rid} de la referencia (columnas {mine})"
        assert col in ref.backbone, f"la referencia no tiene la columna {col}"
    cand: list[tuple[dict, Prot, str]] = []
    for r in rows:
        ch = select_set.domain_chain(r["pdb"], cfg["pfam"], cfg["interpro_pdbe"], r["auth_asym_ids"])
        if ch is None:
            log["skipped"][r["entity"]] = "sin cadena"
            continue
        p = prepare(r["pdb"], ch[0], hmm_path, {"domain": r["domain"], "organism": r["organism"], "resolution": r["resolution"]}, min_pos)
        if p.skip:
            log["skipped"][r["entity"]] = p.skip
            continue
        cand.append((r, p, ch[1]))
    log["counts"]["entidades alineadas"] = len(cand)
    log["counts"]["cadena por PDBe pfam"] = sum(1 for _, _, s in cand if s == "PDBe pfam")
    log["counts"]["cadena por PDBe interpro"] = sum(1 for _, _, s in cand if s == "PDBe interpro")
    log["counts"]["cadena = primera de la entidad"] = sum(1 for _, _, s in cand if s.startswith("primera"))

    # 3. Sin mutantes catalíticos: mutación anotada y residuo distinto del esperado en una columna catalítica
    keep = []
    for r, p, s in cand:
        names = [p.resname[p.pos.index(c)] if c in p.pos else None for c in ref_cols]
        if r["mutation"] and any(n != e for n, e in zip(names, expected)):
            log["skipped"][r["entity"]] = f"mutante catalítico ({r['mutation']}; residuos {names})"
            continue
        keep.append((r, p))
    log["counts"]["sin mutante catalítico"] = len(keep)

    # 4. Agrupado al 90 % de identidad: la mejor resolución de cada grupo
    best = {r["entity"]: r for r in select_set.cluster_best([r for r, _ in keep])}
    chosen = [p for r, p in keep if r["entity"] in best]
    log["counts"]["grupos al 90 % de identidad"] = len(chosen)

    # 5. Ajustes
    ref_full = prepare(ref_code, ref_chain, hmm_path, {"domain": cfg.get("reference_domain", "Unknown")}, min_pos)
    allp = chosen + ([ref_full] if not any(p.pdb == ref_code for p in chosen) else [])
    cat = catalytic_fit(allp, ref_full, ref_cols)
    for p in allp:
        if p.skip:
            log["skipped"][f"{p.pdb}:{p.chain}"] = p.skip
    ligand = reference_substrate(cfg)
    emit = [p for p in chosen if (p.pdb, p.chain) in cat]
    rows_cat = measure(emit, cat, ligand, (np.eye(3), np.zeros(3))) if emit else []
    cat_all = {k: v for k, v in cat.items()}
    acc, core = core_fit(allp, cat_all)
    core_tr = {k: (a[0] @ cat_all[k][0], a[0] @ cat_all[k][1] + a[1]) for k, a in acc.items()}
    rows_core = measure(emit, core_tr, ligand, core_tr[(ref_code, ref_chain)]) if emit else []
    log["counts"]["estructuras en la salida"] = len(emit)
    log["counts"]["posiciones del núcleo"] = len(core)
    log["core_positions"] = core

    data = {
        "family": family, "source": "pipeline", "generated": datetime.date.today().isoformat(),
        "reference_pdb": ref_code, "ligand": cfg["ligand"], "reference_positions": ref_cols,
        "fits": {"catalytic": columns(rows_cat, {}), "core": columns(rows_core, {})},
        "positions": position_stats(rows_cat),
        "horizon": horizon(rows_cat),
    }
    log["entries"] = [{"pdb": p.pdb, "chain": p.chain, "domain": p.domain, "organism": p.organism, "resolution": p.resolution,
                       "residues": len(p.pos)} for p in emit]
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(data, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
    return {"data": data, "log": log}
