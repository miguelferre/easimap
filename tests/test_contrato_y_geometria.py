"""Geometría básica y forma de web/data/gh1.json frente a docs/contrato-datos.md."""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from easimap.geometry import kabsch, spherical, substrate_frame, to_frame

ROOT = Path(__file__).resolve().parents[1]
SALIDA = ROOT / "web" / "data" / "gh1.json"
COLUMNAS = ["pdb", "chain", "auth_seq_id", "hmm_position", "residue", "domain", "theta", "phi", "distance", "x", "y", "z"]


def test_marco_del_sustrato_es_una_rotacion():
    rng = np.random.default_rng(0)
    c1, o, c4 = rng.normal(size=(3, 3)) * 3
    origin, axes = substrate_frame(c1, o, c4)
    assert np.allclose(axes @ axes.T, np.eye(3), atol=1e-9)
    assert np.isclose(np.linalg.det(axes), 1.0)
    assert np.allclose(to_frame(c1, origin, axes), 0)
    p = to_frame(o, origin, axes)[0]
    assert abs(p[0]) < 1e-9 and abs(p[2]) < 1e-9 and p[1] > 0          # el oxígeno glicosídico va sobre +y
    q = to_frame(c4, origin, axes)[0]
    assert abs(q[2]) < 1e-9 and q[0] > 0                                # el C4 queda en el plano z = 0, con x > 0


def test_coordenadas_esfericas():
    d, th, ph = spherical(np.array([0.0, 2.0, 0.0]))
    assert (d, th) == (2.0, 0.0)
    d, th, ph = spherical(np.array([0.0, 0.0, -1.0]))                  # z < 0: phi cerca de 3 pi / 2
    assert math.isclose(th, math.pi / 2) and math.isclose(ph, 1.5 * math.pi)
    assert 0 <= spherical(np.array([-1.0, 0.0, -1e-9]))[2] < 2 * math.pi


def test_quimera_de_dos_dominios_no_tiene_dominio():
    from easimap.select_set import lineage_domain

    def org(*names):
        return {"taxonomy_lineage": [{"name": n} for n in names]}

    assert lineage_domain({"rcsb_entity_source_organism": [org("Bacteria"), org("Bacteria")]}) == "Bacteria"
    assert lineage_domain({"rcsb_entity_source_organism": [org("Bacteria"), org("Eukaryota")]}) == "Unknown"
    assert lineage_domain({"rcsb_entity_source_organism": []}) == "Unknown"


def test_kabsch_recupera_una_rotacion():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(12, 3))
    a = 0.7
    r = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
    y = x @ r.T + np.array([1.0, -2.0, 3.0])
    rot, t = kabsch(x, y)
    assert np.allclose(x @ rot.T + t, y, atol=1e-9)
    assert np.isclose(np.linalg.det(rot), 1.0)


@pytest.fixture(scope="module")
def d():
    return json.loads(SALIDA.read_text(encoding="utf-8"))


@pytest.mark.skipif(not SALIDA.exists(), reason="falta web/data/gh1.json (ejecuta `python -m easimap build GH1`)")
class TestSalida:

    def test_campos_del_contrato(self, d):
        assert d["source"] == "pipeline" and d["family"] == "GH1"
        assert d["reference_pdb"] == "3SCU" and d["reference_positions"] == [165, 360]
        assert set(d["fits"]) == {"catalytic", "core"}
        assert len(d["generated"]) == 10

    @pytest.mark.parametrize("modo", ["catalytic", "core"])
    def test_columnas_paralelas(self, d, modo):
        f = d["fits"][modo]
        assert list(f) == COLUMNAS
        n = len(f["pdb"])
        assert n > 20000 and all(len(f[k]) == n for k in COLUMNAS)
        assert all(isinstance(p, int) and p > 0 for p in f["hmm_position"])
        assert all(0 <= p < 2 * math.pi for p in f["phi"]) and all(0 <= t <= math.pi for t in f["theta"])
        assert set(f["domain"]) <= {"Eukaryota", "Bacteria", "Archaea", "Unknown"}

    @pytest.mark.parametrize("modo", ["catalytic", "core"])
    def test_una_entrada_por_proteina_cadena_y_posicion(self, d, modo):
        f = d["fits"][modo]
        repetidos = [k for k, v in Counter(zip(f["pdb"], f["chain"], f["hmm_position"])).items() if v > 1]
        assert not repetidos

    def test_distancia_coincide_con_xyz(self, d):
        f = d["fits"]["catalytic"]
        dist = np.linalg.norm(np.array([f["x"], f["y"], f["z"]]), axis=0)
        assert np.allclose(dist, f["distance"], atol=2e-3)

    def test_residuos_cataliticos_cerca_del_sustrato(self, d):
        f = d["fits"]["catalytic"]
        for pos in (165, 360):
            dist = [x for x, p in zip(f["distance"], f["hmm_position"]) if p == pos]
            assert np.median(dist) < 6.0

    def test_posiciones_y_horizonte(self, d):
        assert len(d["positions"]) > 400
        assert all({"hmm_position", "n", "rmsf", "distance", "top_residue"} <= set(p) for p in d["positions"])
        h = d["horizon"]
        assert h["grid"] == 26 and len(h["north"]) == len(h["south"]) == 26 * 26
        assert sum(v > 0 for v in h["north"]) > 400
