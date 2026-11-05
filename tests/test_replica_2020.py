"""Réplica de 2020: el pipeline en Python contra la tabla de la memoria.

Empareja por proteína y residuo y exige mediana del error absoluto de la distancia menor de 0,05 Å.
La tolerancia no se toca. Si la tabla de 2020 no está en el disco, las pruebas se omiten.
"""
from __future__ import annotations

from pathlib import Path
import tomllib

import numpy as np
import pytest

from easimap import fetch, replica2020, table2020
from easimap.pfam_archive import archived_hmm

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 0.05  # Å

pytestmark = pytest.mark.skipif(not table2020.available(), reason="falta la tabla de 2020 (carpeta TFM de solo lectura)")


@pytest.fixture(scope="module")
def cfg():
    return tomllib.loads((ROOT / "familias" / "GH1.toml").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def rows2020():
    return table2020.load()


@pytest.fixture(scope="module")
def hmm_2020():
    """PF00232.18 de Pfam 32.0: el modelo vigente cuando se hizo la tabla (453 columnas)."""
    return archived_hmm("PF00232", "32.0")


def _codes(cfg):
    codes = list(cfg["replica2020"]["proteins"])
    for c in codes:
        fetch.pdb_file(c)
    return codes


@pytest.fixture(scope="module")
def run(cfg):
    """Réplica principal: PF00232 vigente, alineamiento local, columnas de referencia por residuo del 3SCU."""
    hmm = fetch.current_hmm("PF00232")
    return replica2020.run(cfg, fetch.CACHE / "pdb", hmm, _codes(cfg), local=True, use_resids=True)


@pytest.fixture(scope="module")
def run_pfam32(cfg, hmm_2020):
    """Control: el HMM de 2020 (Pfam 32.0) con las columnas 165 y 360 tal cual."""
    return replica2020.run(cfg, fetch.CACHE / "pdb", hmm_2020, _codes(cfg), local=True, use_resids=False)


def _errors(tables, rows2020):
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


def test_la_referencia_cierra_con_la_tabla(run, rows2020):
    """3SCU no se ajusta: solo prueba el marco del sustrato y el punto medido."""
    tables, _, _ = run
    per, _, _ = _errors({"3SCU": tables["3SCU"]}, rows2020)
    assert len(per["3SCU"]) == 472
    assert float(np.max(per["3SCU"])) < 1e-3


def test_replica_mediana_menor_de_0_05(run, rows2020):
    tables, skipped, cols = run
    per, allerr, unpaired = _errors(tables, rows2020)
    assert cols[0] < cols[1]
    assert not unpaired, f"proteínas sin emparejar: {unpaired}"
    assert len(allerr) > 20000
    assert float(np.median(allerr)) < TOLERANCE, (
        f"mediana {np.median(allerr):.4f} Å; peores: "
        + ", ".join(f"{k} {np.median(v):.3f}" for k, v in sorted(per.items(), key=lambda kv: -np.median(kv[1]))[:5]))


def test_cobertura_de_proteinas(run, rows2020):
    """Al menos 55 de las 64 proteínas de la tabla se reproducen. Las demás se documentan en el hallazgo."""
    tables, _, _ = run
    per, _, _ = _errors(tables, rows2020)
    assert len(per) >= 55
    assert all(float(np.median(e)) < TOLERANCE for e in per.values())


def test_control_con_el_hmm_de_2020(run_pfam32, rows2020):
    """Con Pfam 32.0 las columnas 165 y 360 son las de 2020 y el resultado es el mismo."""
    tables, _, cols = run_pfam32
    per, allerr, unpaired = _errors(tables, rows2020)
    assert cols == [165, 360]
    assert not unpaired
    assert float(np.median(allerr)) < TOLERANCE


def test_el_alineamiento_glocal_no_cierra(cfg, rows2020):
    """Documenta el fallo: el hmmalign de pyhmmer es glocal y mueve las columnas catalíticas."""
    hmm = fetch.current_hmm("PF00232")
    tables, _, _ = replica2020.run(cfg, fetch.CACHE / "pdb", hmm, _codes(cfg), local=False, use_resids=True)
    _, allerr, _ = _errors(tables, rows2020)
    assert float(np.median(allerr)) > 1.0
