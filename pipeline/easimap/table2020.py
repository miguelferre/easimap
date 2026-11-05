"""Lectura de la tabla de 2020 y emparejamiento con las tablas del pipeline.

La tabla de la memoria no guarda el número de residuo del PDB (resid). Se recupera del orden de filas:
2020 iteraba `lsort -uniq` sobre los resid, que ordena como texto ("10" antes que "5"). Con ese orden
la fila i de una proteína es el resid i-ésimo. Cada emparejamiento se valida con el nombre del residuo.
"""
from __future__ import annotations

import csv
import collections
import os
from pathlib import Path

# Tabla de la memoria (referencias 165 y 360) y tabla con los códigos PDB completos, mismo orden de filas.
# Son de solo lectura y no se copian al repo.
DEFAULT_ROOT = Path(os.environ.get("EASIMAP_TFM_ROOT", "datos-2020-no-publicados"))
MEMORIA_TABLE = DEFAULT_ROOT / "Tfm" / "EASiMap" / "angulosydistancias" / "GH1.tsv"
CODES_TABLE = DEFAULT_ROOT / "easimap" / "familias" / "GH1.tsv"


def single_table(family: str) -> Path:
    """Tabla de 2020 de una familia que no tiene copia en la memoria (GH5, GH16, GH18, GH20, CE4)."""
    return DEFAULT_ROOT / "easimap" / "familias" / f"{family}.tsv"


def available(family: str = "GH1") -> bool:
    if family == "GH1":
        return MEMORIA_TABLE.exists() and CODES_TABLE.exists()
    return single_table(family).exists()


def _rows(path: Path) -> list[list[str]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        return [[c.strip() for c in r] for r in list(csv.reader(fh))[1:]]


def load(family: str = "GH1") -> dict[str, list[dict]]:
    """Filas de la tabla de 2020 agrupadas por código PDB de 4 letras, en el orden del fichero.

    GH1 tiene dos tablas (la de la memoria y la de los códigos completos). Las demás familias tienen una sola,
    que lleva a la vez el código y las medidas.
    """
    if family == "GH1":
        memoria, codes = _rows(MEMORIA_TABLE), _rows(CODES_TABLE)
    else:
        memoria = codes = _rows(single_table(family))
    if len(memoria) != len(codes):
        raise ValueError("las dos tablas de 2020 no tienen el mismo número de filas")
    out: dict[str, list[dict]] = collections.OrderedDict()
    for m, c in zip(memoria, codes):
        pos = m[1]
        out.setdefault(c[0][:4].upper(), []).append({
            "position": int(float(pos)) if " " not in pos else None,   # "0.0 170.0": resid ambiguo
            "residue": m[3], "distance": float(m[4]), "theta": float(m[5]), "phi": float(m[6]),
            "xyz": (float(m[7]), float(m[8]), float(m[9])),
        })
    return out


def pair(table, rows: list[dict], max_mismatch: float = 0.01):
    """Empareja las filas de 2020 con una ChainTable del pipeline.

    Devuelve (pares, estado). Cada par es (fila 2020, índice en `table`). Estado: "ok", "doble" si la
    proteína figura dos veces en la tabla de 2020 (se usa la primera copia), "nombres" si no coinciden los
    residuos, o "ok-sin-N" / "doble-sin-N" si N filas tienen otro nombre de residuo y se descartan.
    Una entrada del PDB revisada desde 2020 puede cambiar un residuo (5HOS: Gly 25 es hoy Pro). Se tolera
    hasta `max_mismatch` de las filas, y las filas distintas no entran en la comparación.
    """
    n = len(table.resid)
    state = "ok"
    if len(rows) == 2 * n:
        state = "doble"
    elif len(rows) != n:
        return [], "longitud"
    wrong = {i for i, r in enumerate(rows[:n]) if r["residue"] != table.resname[i]}
    if len(wrong) > max_mismatch * n:
        return [], "nombres"
    if wrong:
        state += f"-sin-{len(wrong)}"
    return [(rows[i], i) for i in range(n) if i not in wrong], state
