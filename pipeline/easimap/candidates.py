"""Candidatas a historia: posiciones HMM cuyo residuo se coloca distinto en dos grupos de dominios.

Para cada posición compara el punto medio del grupo A (eucariotas) con el del grupo B (procariotas: bacterias y
arqueas), con la superposición catalítica de `web/data/<familia>.json`. El estadístico es la distancia entre los
dos puntos medios, en Å. El valor p sale de una prueba de permutación: se reparten al azar las etiquetas de
dominio entre las ESTRUCTURAS (no entre filas), porque todas las posiciones de una estructura comparten etiqueta.
Se corrige por comparaciones múltiples con Benjamini-Hochberg (FDR) y, aparte, con el máximo de las permutaciones
(Westfall-Young), que controla el error de familia.

Son hipótesis. No hay corrección filogenética: las estructuras de un mismo clado no son puntos independientes,
y el agrupado al 90 % de identidad solo reduce el problema.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

EUK = ("Eukaryota",)
PROK = ("Bacteria", "Archaea")


def _structures(fit: dict):
    """Estructuras (pdb, cadena) en orden de aparición y, por fila, el índice de su estructura."""
    keys, idx = {}, []
    for p, c in zip(fit["pdb"], fit["chain"]):
        idx.append(keys.setdefault((p, c), len(keys)))
    return list(keys), np.array(idx)


def scan(fit: dict, group_a=EUK, group_b=PROK, n_perm: int = 5000, seed: int = 20201104, min_per_group: int = 4) -> dict:
    """Prueba de permutación por posición. Devuelve {resumen, posiciones} ordenadas por valor q."""
    keys, sidx = _structures(fit)
    domain_of = {}
    for k, i, d in zip(zip(fit["pdb"], fit["chain"]), sidx, fit["domain"]):
        domain_of[int(i)] = d
    label = np.array([1 if domain_of[i] in group_a else (0 if domain_of[i] in group_b else -1) for i in range(len(keys))])
    use = label >= 0
    n_a, n_b = int((label == 1).sum()), int((label == 0).sum())
    summary = {"grupo_a": list(group_a), "grupo_b": list(group_b), "estructuras_a": n_a, "estructuras_b": n_b,
               "estructuras_sin_grupo": int((label < 0).sum()), "permutaciones": n_perm, "semilla": seed}
    if min(n_a, n_b) < min_per_group:
        summary["aviso"] = "grupo con menos de %d estructuras: no se hace la prueba" % min_per_group
        return {"resumen": summary, "posiciones": []}

    rng = np.random.default_rng(seed)
    n_struct = len(keys)
    # Etiquetas permutadas solo entre las estructuras con grupo; las demás quedan fuera (−1) en todas las permutaciones.
    base = label.copy()
    members = np.flatnonzero(use)
    perms = np.empty((n_perm, n_struct), dtype=np.int8)
    perms[:] = -1
    for b in range(n_perm):
        perms[b, members] = rng.permutation(base[members])
    is_a = (perms == 1).astype(np.float32)
    is_b = (perms == 0).astype(np.float32)

    pts = np.column_stack([fit["x"], fit["y"], fit["z"]]).astype(np.float64)
    pos = np.array(fit["hmm_position"])
    out = []
    for q in np.unique(pos):
        rows = np.flatnonzero((pos == q) & use[sidx])
        if len(rows) == 0:
            continue
        s = sidx[rows]
        X = pts[rows]
        na0, nb0 = int((label[s] == 1).sum()), int((label[s] == 0).sum())
        if na0 < min_per_group or nb0 < min_per_group:
            continue
        a, b = is_a[:, s], is_b[:, s]               # permutaciones x estructuras con esa posición
        na, nb = a.sum(1), b.sum(1)
        ok = (na >= 1) & (nb >= 1)
        sa, sb = a @ X, b @ X
        with np.errstate(invalid="ignore", divide="ignore"):
            stat = np.linalg.norm(sa / na[:, None] - sb / nb[:, None], axis=1)
        obs_a = X[label[s] == 1].mean(0)
        obs_b = X[label[s] == 0].mean(0)
        obs = float(np.linalg.norm(obs_a - obs_b))
        # Las permutaciones usan los mismos grupos de estructuras; la observada es la de las etiquetas reales.
        null = np.where(ok, stat, 0.0)
        out.append({"hmm_position": int(q), "n_a": na0, "n_b": nb0, "separacion": obs,
                    "p": float((1 + (null >= obs).sum()) / (n_perm + 1)), "_null": null})

    if not out:
        return {"resumen": summary, "posiciones": []}
    # Benjamini-Hochberg
    p = np.array([r["p"] for r in out])
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for rank, i in zip(range(m, 0, -1), order[::-1]):
        prev = min(prev, p[i] * m / rank)
        q[i] = prev
    # Westfall-Young: p ajustado con el máximo de la separación entre posiciones en cada permutación.
    # La separación se pondera por posición para que el máximo no lo gobiernen las posiciones con pocas estructuras:
    # se usa la propia distancia en Å, que es lo que se interpreta.
    nulls = np.vstack([r["_null"] for r in out])
    maxnull = nulls.max(axis=0)
    for r, qi in zip(out, q):
        r["q_bh"] = float(qi)
        r["p_familia"] = float((1 + (maxnull >= r["separacion"]).sum()) / (n_perm + 1))
        del r["_null"]
        r["separacion"] = round(r["separacion"], 3)
    out.sort(key=lambda r: (r["q_bh"], -r["separacion"]))
    summary["posiciones_probadas"] = m
    summary["q_bh_menor_de_0_05"] = int(sum(r["q_bh"] < 0.05 for r in out))
    summary["p_familia_menor_de_0_05"] = int(sum(r["p_familia"] < 0.05 for r in out))
    return {"resumen": summary, "posiciones": out}


def run(data_path: Path, out_path: Path, **kw) -> dict:
    data = json.loads(Path(data_path).read_text(encoding="utf-8"))
    res = scan(data["fits"]["catalytic"], **kw)
    res["familia"] = data["family"]
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
    return res


if __name__ == "__main__":
    # python -m easimap.candidates GH5
    import sys
    fam = sys.argv[1] if len(sys.argv) > 1 else "GH1"
    root = Path(__file__).resolve().parents[2]
    res = run(root / "web" / "data" / f"{fam.lower()}.json", root / "data" / f"{fam.lower()}.candidatas.json")
    print(res["resumen"])
