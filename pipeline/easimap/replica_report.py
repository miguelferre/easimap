"""Cifras de la réplica de 2020: error por proteína y global, con tres variantes de HMM y alineamiento."""
from __future__ import annotations

import json
import tomllib
from pathlib import Path

import numpy as np

from . import fetch, replica2020, table2020
from .pfam_archive import archived_hmm

ROOT = Path(__file__).resolve().parents[2]

def variants(pfam: str, release: str) -> dict:
    """Las tres variantes de HMM y alineamiento que se comparan: (fichero, local, descripción)."""
    return {
        f"pfam{release}_local": (f"{pfam}_pfam{release}.hmm", True, f"{pfam} de Pfam {release}, alineamiento local"),
        "tabla_2020": (f"{pfam}_pfam{release}.hmm", True, "posiciones de la propia tabla de 2020 (sin HMM): prueba el marco y la superposición"),
        "actual_local": (f"{pfam}.hmm", True, f"{pfam} vigente, alineamiento local"),
        "actual_glocal": (f"{pfam}.hmm", False, f"{pfam} vigente, alineamiento glocal (hmmalign de pyhmmer)"),
    }


def default_out(family: str) -> Path:
    """GH1 conserva su nombre original; las demás familias llevan el suyo."""
    name = "resumen.json" if family == "GH1" else f"resumen-{family.lower()}.json"
    return ROOT / "data" / "replica2020" / name


def _stats(tables, rows2020):
    per, allerr, unpaired = {}, [], {}
    for code, t in tables.items():
        if code not in rows2020:
            continue
        pairs, state = table2020.pair(t, rows2020[code])
        if not pairs:
            unpaired[code] = state
            continue
        e = np.abs(np.array([t.distance[i] for _, i in pairs]) - np.array([r["distance"] for r, _ in pairs]))
        per[code] = e
        allerr.extend(e)
    return per, np.array(allerr), unpaired


def report(out: Path | None = None, family: str = "GH1") -> dict:
    cfg = tomllib.loads((ROOT / "familias" / f"{family}.toml").read_text(encoding="utf-8"))
    codes = list(cfg["replica2020"]["proteins"])
    for c in codes:
        fetch.pdb_file(c)
    hmm_dir = fetch.CACHE / "hmm"
    release = cfg["hmm_release"]
    archived_hmm(cfg["pfam"], release)
    fetch.current_hmm(cfg["pfam"])
    rows2020 = table2020.load(family)
    total_rows = sum(len(v) for v in rows2020.values())
    res: dict = {"familia": family, "filas_2020": total_rows, "proteinas_2020": len(rows2020), "variantes": {}}
    for key, (hmm_name, local, label) in variants(cfg["pfam"], release).items():
        use_resids = (not key.startswith("pfam") or not cfg["replica2020"].get("numbering_2020", True)) and key != "tabla_2020"
        try:
            tables, skipped, cols = replica2020.run(
                cfg, fetch.CACHE / "pdb", hmm_dir / hmm_name, codes, local=local, use_resids=use_resids,
                table_rows=rows2020 if key == "tabla_2020" else None)
        except ValueError as exc:     # p. ej. un residuo de referencia fuera del modelo (CE4: His295)
            res["variantes"][key] = {"descripcion": label, "error": str(exc)}
            continue
        per, allerr, unpaired = _stats(tables, rows2020)
        in_table = [c for c in rows2020 if c not in tables]
        res["variantes"][key] = {
            "descripcion": label,
            "columnas_referencia": cols,
            "proteinas_reproducidas": len(per),
            "omitidas": {s.pdb: s.reason for s in skipped},
            "en_tabla_2020_y_no_reproducidas": in_table,
            "sin_emparejar": unpaired,
            "pares": int(len(allerr)),
            "error_mediana_A": float(np.median(allerr)) if len(allerr) else None,
            "error_p90_A": float(np.percentile(allerr, 90)) if len(allerr) else None,
            "error_p99_A": float(np.percentile(allerr, 99)) if len(allerr) else None,
            "error_max_A": float(allerr.max()) if len(allerr) else None,
            "mediana_por_proteina": {k: float(np.median(v)) for k, v in per.items()},
        }
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    return res



if __name__ == "__main__":
    # python -m easimap.replica_report GH5   (hasta que `easimap replica` acepte la familia)
    import sys
    fam = sys.argv[1] if len(sys.argv) > 1 else "GH1"
    res = report(default_out(fam), fam)
    for k, v in res["variantes"].items():
        print(k, v["proteinas_reproducidas"], "proteínas", v["pares"], "pares, mediana", f"{v['error_mediana_A']:.5f}", "Å")
