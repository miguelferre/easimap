"""Prueba del desfase en las posiciones 334 a 336 de GH1.

Hipótesis del plan: la posición 333 (Trp) separa eucariotas y procariotas porque un lazo cercano tiene un residuo
de diferencia, y el HMM coloca el hueco en 334-336.

Dos pruebas:
1. Pareja: alinear un eucariota y un procariota entre 330 y 340 con el HMM y ver dónde cae el hueco. Después,
   con las dos estructuras superpuestas por el núcleo, mirar qué desplazamiento de residuos las acerca más.
2. Familia: para cada posición p de 328 a 342, distancia entre el punto medio de los eucariotas en p y el de los
   procariotas en p + s, con s de -2 a +2.
"""
from __future__ import annotations

import collections
import json
from pathlib import Path

import numpy as np

from . import build, fetch
from .geometry import kabsch

LO, HI = 330, 340


def segment(p: build.Prot, lo: int = LO, hi: int = HI) -> list[dict]:
    """Posiciones HMM de lo a hi: residuo y número de autor, o hueco si el HMM no puso residuo."""
    out = []
    for q in range(lo, hi + 1):
        if q in p.ca:
            i = p.pos.index(q)
            out.append({"hmm": q, "resid": p.resid[i], "res": p.resname[i]})
        else:
            out.append({"hmm": q, "resid": None, "res": "-"})
    return out


def pair_test(euk: str, prok: str, hmm_path: Path, cols=(165, 360)) -> dict:
    """Alinea una pareja entre 330 y 340 y busca el desplazamiento de residuos que más acerca las dos estructuras."""
    a = build.prepare(euk, "A", hmm_path, {"domain": "Eukaryota"})
    b = build.prepare(prok, "A", hmm_path, {"domain": "Bacteria"})
    ref = build.prepare("3SCU", "A", hmm_path, {"domain": "Eukaryota"})
    fit = build.catalytic_fit([a, b], ref, list(cols))
    ka, kb = (a.pdb, a.chain), (b.pdb, b.chain)
    ra, ta = fit[ka]
    rb, tb = fit[kb]
    # Cadena completa en orden de secuencia: residuos con posición, para medir con desplazamiento real de residuos
    seg_a, seg_b = segment(a), segment(b)
    ca_a = {q: ra @ v + ta for q, v in a.ca.items()}
    ca_b = {q: rb @ v + tb for q, v in b.ca.items()}
    # Desplazamiento en residuos consecutivos del procariota (no en columnas HMM)
    result = []
    for s in range(-3, 4):
        d = []
        for q in range(LO, HI + 1):
            if q not in ca_a:
                continue
            # residuo del procariota s posiciones de secuencia después del que el HMM puso en q
            if q not in b.pos:
                continue
            j = b.pos.index(q) + s
            if 0 <= j < len(b.pos):
                d.append(float(np.linalg.norm(ca_a[q] - (rb @ b.ca[b.pos[j]] + tb))))
        result.append({"shift": s, "n": len(d), "mean_CA_distance": round(float(np.mean(d)), 2) if d else None})
    per_pos = []
    for q in range(LO - 2, HI + 3):
        if q in ca_a and q in ca_b:
            per_pos.append({"hmm": q, "ca_distance": round(float(np.linalg.norm(ca_a[q] - ca_b[q])), 2)})
    return {"eukaryote": euk, "prokaryote": prok, "segment_eukaryote": seg_a, "segment_prokaryote": seg_b,
            "shift_scan": result, "ca_distance_by_position": per_pos}


def family_test(data: dict, lo: int = 326, hi: int = 344, shifts=range(-2, 3)) -> list[dict]:
    """Para cada posición, qué desplazamiento acerca más los puntos medios de eucariotas y procariotas."""
    f = data["fits"]["catalytic"]
    groups: dict[tuple[str, int], list] = collections.defaultdict(list)
    for q, dom, x, y, z in zip(f["hmm_position"], f["domain"], f["x"], f["y"], f["z"]):
        g = "euk" if dom == "Eukaryota" else "prok" if dom in ("Bacteria", "Archaea") else None
        if g:
            groups[(g, q)].append((x, y, z))
    mean = {k: np.mean(v, axis=0) for k, v in groups.items() if len(v) >= 3}
    out = []
    for q in range(lo, hi + 1):
        if ("euk", q) not in mean:
            continue
        row = {"hmm": q, "n_euk": len(groups[("euk", q)])}
        best = None
        for s in shifts:
            if ("prok", q + s) in mean:
                d = float(np.linalg.norm(mean[("euk", q)] - mean[("prok", q + s)]))
                row[f"d_shift_{s:+d}"] = round(d, 2)
                if best is None or d < best[1]:
                    best = (s, d)
        row["best_shift"] = best[0] if best else None
        out.append(row)
    return out


def structural_map(euk: str, prok: str, hmm_path: Path, lo: int = 331, hi: int = 338) -> list[dict]:
    """Para cada residuo del eucariota en la ventana, el residuo más cercano (CA) del procariota.

    Las dos estructuras van superpuestas por los residuos catalíticos. Incluye los residuos sin columna HMM
    que caen dentro de la ventana, que es donde el HMM pone el hueco.
    """
    ref = build.prepare("3SCU", "A", hmm_path, {"domain": "Eukaryota"})
    a = build.prepare(euk, "A", hmm_path, {})
    b = build.prepare(prok, "A", hmm_path, {})
    fit = build.catalytic_fit([a, b], ref, [165, 360])
    if (a.pdb, a.chain) not in fit or (b.pdb, b.chain) not in fit:
        return []
    (ra, ta), (rb, tb) = fit[(a.pdb, a.chain)], fit[(b.pdb, b.chain)]
    chain_a, chain_b = build.read_chain(euk, "A"), build.read_chain(prok, "A")
    hm_a, hm_b = dict(zip(a.resid, a.pos)), dict(zip(b.resid, b.pos))
    xb = np.array([rb @ r[3]["CA"] + tb for r in chain_b])
    out = []
    for num, _, name, atoms in chain_a:
        q = hm_a.get(num, 0)
        near = any(lo <= hm_a.get(m, 0) <= hi for m in (num - 1, num + 1))
        if lo <= q <= hi or (q == 0 and near):
            d = np.linalg.norm(xb - (ra @ atoms["CA"] + ta), axis=1)
            j = int(d.argmin())
            out.append({"euk": f"{name}{num}", "euk_hmm": q or None, "prok": f"{chain_b[j][2]}{chain_b[j][0]}",
                        "prok_hmm": hm_b.get(chain_b[j][0]) or None, "ca_distance": round(float(d[j]), 2)})
    return out


def run(data_path: Path, out: Path, pairs=(("1CBG", "1BGA"), ("1E1E", "1GNX"), ("3AHZ", "1OD0"), ("2RGL", "1QOX"))) -> dict:
    hmm = fetch.CACHE / "hmm" / "PF00232_pfam32.0.hmm"
    data = json.loads(data_path.read_text(encoding="utf-8"))
    res = {"hmm": hmm.name, "family_test": family_test(data), "pairs": []}
    for e, p in pairs:
        res["pairs"].append({"eukaryote": e, "prokaryote": p, "structural_map": structural_map(e, p, hmm),
                             **{k: v for k, v in pair_test(e, p, hmm).items() if k in ("segment_eukaryote", "segment_prokaryote", "shift_scan")}})
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    return res
