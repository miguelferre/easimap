"""Pruebas por familia: marco del sustrato, réplica de 2020 y forma de la salida.

Se añade una entrada a FAMILIAS por cada familia con la réplica cerrada. La tolerancia (0,05 Å, mediana) es la
misma para todas y no se toca. Si la tabla de 2020 no está en el disco, las pruebas de réplica se omiten.
"""
from __future__ import annotations

import json
import math
import tomllib
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from easimap import fetch, replica2020, table2020
from easimap.geometry import substrate_frame, to_frame
from easimap.pfam_archive import archived_hmm
from easimap.structure import atom_from_spec, frame_atoms

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 0.05  # Å
COLUMNAS = ["pdb", "chain", "auth_seq_id", "hmm_position", "residue", "domain", "theta", "phi", "distance", "x", "y", "z"]

# familia: columnas catalíticas de la salida, proteínas que se reproducen como mínimo y estructuras mínimas en la salida.
# modelo "vigente": réplica con el HMM actual (GH1, GH5). modelo "pfam32": réplica con el HMM de Pfam 32.0 y columnas por resid
# (el modelo de 2020 de GH18 y GH20 no se identificó; GH16 sí cierra con Pfam 32.0 y las columnas de 2020).
FAMILIAS = {
    "GH5": dict(cols=[130, 243], min_proteinas=70, min_salida=50, modelo="vigente", salida=True, max_dist=6.0),
    "GH16": dict(cols=[73, 75, 77], min_proteinas=38, min_salida=40, modelo="pfam32", salida=True, tabla=False, max_dist=6.0),
    "GH18": dict(cols=[109, 111, 113], min_proteinas=75, min_salida=60, modelo="pfam32", salida=True, max_dist=12.5,
                 fuera_de_tolerancia=["3CZ8"]),   # 3CZ8: mediana de 0,9 Å
    "GH20": dict(cols=[169, 170], min_proteinas=17, min_salida=20, modelo="pfam32", salida=True, max_dist=8.0),
    "CE4": dict(cols=[8, 64, 163], min_proteinas=19, min_salida=0, modelo="ninguno", salida=False),
}
CON_REPLICA_POR_MODELO = [f for f, v in FAMILIAS.items() if v["modelo"] != "ninguno"]


def _cfg(fam):
    return tomllib.loads((ROOT / "familias" / f"{fam}.toml").read_text(encoding="utf-8"))


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


def _codes(cfg):
    codes = list(cfg["replica2020"]["proteins"])
    for c in codes:
        fetch.pdb_file(c)
    return codes


# ---------- Marco del sustrato (generalizado) ----------

@pytest.mark.parametrize("fam", ["GH1", *FAMILIAS])
def test_el_marco_de_la_configuracion_es_una_rotacion_valida(fam):
    """Los tres átomos del TOML existen en la referencia y dan un marco ortonormal con x > 0 para el tercer átomo."""
    cfg = _cfg(fam)
    c0, cy, cx = frame_atoms(fetch.pdb_file(cfg["reference_pdb"]), cfg["substrate"])
    origin, axes = substrate_frame(c0, cy, cx)
    assert np.allclose(axes @ axes.T, np.eye(3), atol=1e-9) and np.isclose(np.linalg.det(axes), 1.0)
    assert to_frame(cy, origin, axes)[0][1] > 0
    assert to_frame(cx, origin, axes)[0][0] > 0


def test_atomo_con_altloc_distingue_conformaciones():
    """En GH5 el origen es el C1 de la conformación A y el eje y sale del C1 de la B, a 0,495 Å."""
    cfg = _cfg("GH5")
    pdb = fetch.pdb_file(cfg["reference_pdb"])
    a = atom_from_spec(pdb, cfg["substrate"]["origin"])
    b = atom_from_spec(pdb, cfg["substrate"]["y_axis"])
    assert np.linalg.norm(b - a) == pytest.approx(0.495, abs=0.002)


@pytest.mark.skipif(not table2020.DEFAULT_ROOT.exists(), reason="falta la carpeta del TFM")
@pytest.mark.parametrize("fam", ["GH1", *FAMILIAS])
def test_el_marco_reproduce_el_ligando_de_2020(fam):
    """Llevar el ligando original de 2020 al marco de la configuración da su fichero `_1` (salida de VMD): menos de 0,002 Å."""
    cfg = _cfg(fam)
    lig = table2020.DEFAULT_ROOT / "Tfm" / "EASiMap" / "Ligandos"
    src, dst = lig / f"{cfg['ligand_file']}.pdb", lig / f"{cfg['ligand_file']}_1.pdb"
    if not (src.exists() and dst.exists()):
        pytest.skip("faltan los ficheros del ligando")

    def coords(p):
        return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(p, encoding="latin-1")
                         if l[:6] in ("ATOM  ", "HETATM")])

    origin, axes = replica2020.substrate_from_reference(fetch.pdb_file(cfg["reference_pdb"]), cfg)
    err = np.abs(to_frame(coords(src), origin, axes) - coords(dst)).max()
    assert err < 0.002


# ---------- Réplica de 2020 ----------

@pytest.mark.parametrize("fam", CON_REPLICA_POR_MODELO)
class TestReplica:

    @pytest.fixture(autouse=True)
    def _tabla(self, fam):
        if not table2020.available(fam):
            pytest.skip("falta la tabla de 2020 (carpeta TFM de solo lectura)")

    @pytest.fixture
    def rows2020(self, fam):
        return table2020.load(fam)

    @pytest.fixture
    def run_actual(self, fam):
        """HMM vigente, alineamiento local, columnas de referencia por residuo de la estructura de referencia."""
        cfg = _cfg(fam)
        if FAMILIAS[fam]["modelo"] == "pfam32":
            hmm = archived_hmm(cfg["pfam"], cfg["hmm_release"])
            return replica2020.run(cfg, fetch.CACHE / "pdb", hmm, _codes(cfg), local=True,
                                   use_resids=not cfg["replica2020"].get("numbering_2020", True))
        return replica2020.run(cfg, fetch.CACHE / "pdb", fetch.current_hmm(cfg["pfam"]), _codes(cfg), local=True, use_resids=True)

    @pytest.fixture
    def run_2020(self, fam):
        """Control: el HMM de la versión de Pfam de 2020 con las columnas de 2020 tal cual."""
        cfg = _cfg(fam)
        hmm = archived_hmm(cfg["pfam"], cfg["hmm_release"])
        return replica2020.run(cfg, fetch.CACHE / "pdb", hmm, _codes(cfg), local=True, use_resids=False)

    def test_mediana_menor_de_la_tolerancia(self, run_actual, rows2020, fam):
        tables, _, cols = run_actual
        per, allerr, unpaired = _errors(tables, rows2020)
        assert not unpaired, f"proteínas sin emparejar: {unpaired}"
        assert len(allerr) > 100 * FAMILIAS[fam]["min_proteinas"], "demasiados pocos pares"
        assert float(np.median(allerr)) < TOLERANCE, (
            f"mediana {np.median(allerr):.4f} Å; peores: "
            + ", ".join(f"{k} {np.median(v):.3f}" for k, v in sorted(per.items(), key=lambda kv: -np.median(kv[1]))[:5]))

    def test_cobertura_de_proteinas(self, run_actual, rows2020, fam):
        tables, _, _ = run_actual
        per, _, _ = _errors(tables, rows2020)
        assert len(per) >= FAMILIAS[fam]["min_proteinas"]
        fuera = {k for k, e in per.items() if float(np.median(e)) >= TOLERANCE}
        assert fuera <= set(FAMILIAS[fam].get("fuera_de_tolerancia", [])), f"proteínas fuera de tolerancia: {fuera}"

    def test_control_con_el_hmm_de_2020(self, run_2020, rows2020, fam):
        """Con la versión de Pfam de 2020 las columnas catalíticas son las de la tabla y el resultado es el mismo."""
        if not _cfg(fam)["replica2020"].get("numbering_2020", True):
            pytest.skip("el HMM de Pfam 32.0 no reproduce la numeración de 2020 en esta familia")
        tables, _, cols = run_2020
        per, allerr, unpaired = _errors(tables, rows2020)
        assert cols == FAMILIAS[fam]["cols"]
        assert not unpaired
        assert float(np.median(allerr)) < TOLERANCE

    def test_el_alineamiento_glocal_no_cierra(self, fam, rows2020):
        """Documenta el fallo: el hmmalign de pyhmmer es glocal y mueve las columnas catalíticas."""
        cfg = _cfg(fam)
        hmm = (archived_hmm(cfg["pfam"], cfg["hmm_release"]) if FAMILIAS[fam]["modelo"] == "pfam32"
               else fetch.current_hmm(cfg["pfam"]))
        tables, _, _ = replica2020.run(cfg, fetch.CACHE / "pdb", hmm, _codes(cfg), local=False, use_resids=True)
        _, allerr, _ = _errors(tables, rows2020)
        assert float(np.median(allerr)) > 1.0


@pytest.mark.parametrize("fam", list(FAMILIAS))
def test_replica_con_las_posiciones_de_la_tabla(fam):
    """Sin modelo HMM: las posiciones de la propia tabla de 2020 fijan los residuos de referencia. Prueba el marco, la
    superposición y el punto medido. CE4 solo se comprueba así."""
    if not table2020.available(fam):
        pytest.skip("falta la tabla de 2020 (carpeta TFM de solo lectura)")
    if not FAMILIAS[fam].get("tabla", True):
        pytest.skip("en esta familia la fila del residuo de referencia es ambigua en la tabla (GH16: 4ATF, '72.0 73.0')")
    cfg = _cfg(fam)
    rows = table2020.load(fam)
    hmm = archived_hmm(cfg["pfam"], cfg["hmm_release"])
    tables, _, cols = replica2020.run(cfg, fetch.CACHE / "pdb", hmm, _codes(cfg), local=True, table_rows=rows)
    per, allerr, unpaired = _errors(tables, rows)
    assert cols == cfg["replica2020"]["reference_positions_2020"]
    assert not unpaired and len(per) >= FAMILIAS[fam]["min_proteinas"] - 3
    assert float(np.median(allerr)) < TOLERANCE


def test_pair_tolera_un_residuo_revisado_y_rechaza_muchos():
    """5HOS: un Pro actual donde 2020 tenía Gly. Se descarta esa fila. Más de un 1 % de nombres distintos no se empareja."""
    t = replica2020.ChainTable("XXXX")
    t.resid = list(range(200))
    t.resname = ["ALA"] * 200
    rows = [{"residue": "ALA", "distance": 1.0} for _ in range(200)]
    rows[7]["residue"] = "GLY"
    pares, estado = table2020.pair(t, rows)
    assert len(pares) == 199 and estado == "ok-sin-1"
    for i in (1, 2, 3):
        rows[i]["residue"] = "GLY"
    assert table2020.pair(t, rows) == ([], "nombres")


# ---------- Salida ----------

@pytest.mark.parametrize("fam", [f for f, v in FAMILIAS.items() if v["salida"]])
class TestSalida:

    @pytest.fixture
    def d(self, fam):
        path = ROOT / "web" / "data" / f"{fam.lower()}.json"
        if not path.exists():
            pytest.skip(f"falta {path.name} (ejecuta `python -m easimap build {fam}`)")
        return json.loads(path.read_text(encoding="utf-8"))

    def test_campos_del_contrato(self, d, fam):
        cfg = _cfg(fam)
        assert d["source"] == "pipeline" and d["family"] == fam
        assert d["reference_pdb"] == cfg["reference_pdb"] and d["reference_positions"] == FAMILIAS[fam]["cols"]
        assert set(d["fits"]) == {"catalytic", "core"} and len(d["generated"]) == 10

    @pytest.mark.parametrize("modo", ["catalytic", "core"])
    def test_columnas_paralelas(self, d, modo, fam):
        f = d["fits"][modo]
        assert list(f) == COLUMNAS
        n = len(f["pdb"])
        assert n > 100 * FAMILIAS[fam]["min_salida"] and all(len(f[k]) == n for k in COLUMNAS)
        assert all(isinstance(p, int) and p > 0 for p in f["hmm_position"])
        assert all(0 <= p < 2 * math.pi for p in f["phi"]) and all(0 <= t <= math.pi for t in f["theta"])
        assert set(f["domain"]) <= {"Eukaryota", "Bacteria", "Archaea", "Unknown"}

    @pytest.mark.parametrize("modo", ["catalytic", "core"])
    def test_una_entrada_por_proteina_cadena_y_posicion(self, d, modo):
        f = d["fits"][modo]
        assert not [k for k, v in Counter(zip(f["pdb"], f["chain"], f["hmm_position"])).items() if v > 1]

    def test_distancia_coincide_con_xyz(self, d):
        f = d["fits"]["catalytic"]
        dist = np.linalg.norm(np.array([f["x"], f["y"], f["z"]]), axis=0)
        assert np.allclose(dist, f["distance"], atol=2e-3)

    def test_residuos_cataliticos_cerca_del_sustrato(self, d, fam):
        f = d["fits"]["catalytic"]
        for pos in FAMILIAS[fam]["cols"]:
            dist = [x for x, p in zip(f["distance"], f["hmm_position"]) if p == pos]
            assert dist and np.median(dist) < FAMILIAS[fam]["max_dist"]

    def test_estructuras_y_horizonte(self, d, fam):
        assert len(set(d["fits"]["catalytic"]["pdb"])) >= FAMILIAS[fam]["min_salida"]
        assert set(d["fits"]["catalytic"]["pdb"]) == set(d["fits"]["core"]["pdb"])
        h = d["horizon"]
        assert h["grid"] == 26 and len(h["north"]) == len(h["south"]) == 26 * 26
        assert sum(v > 0 for v in h["north"]) > 400
        assert len(d["positions"]) > 150
        assert all({"hmm_position", "n", "rmsf", "distance", "top_residue"} <= set(p) for p in d["positions"])


# ---------- D-017: O1 reconstruido del anómero beta ----------

def _angulo_o1(code, cadena, anillo, vecino):
    from easimap.geometry import beta_o1_direction
    pdb = fetch.pdb_file(code)

    def at(res, atom):
        return atom_from_spec(pdb, dict(chain=cadena, resseq=res, resname="BGC", atom=atom))

    c1 = at(anillo, "C1")
    real = at(vecino, "O4") - c1
    real /= np.linalg.norm(real)
    d = beta_o1_direction(c1, at(anillo, "C2"), at(anillo, "O5"))
    assert np.isclose(np.linalg.norm(d), 1.0)
    return float(np.degrees(np.arccos(np.clip(d @ real, -1, 1))))


def test_o1_reconstruido_acierta_en_anillos_sin_distorsion():
    """Anillos beta de glucosa con vecino enlazado: el O1 reconstruido cae a menos de 5° del O4 real."""
    for code, anillo, vecino in (("3SCU", 4, 3), ("3SCU", 3, 2), ("3SCU", 2, 1), ("3ZMR", 4, 3), ("3ZMR", 3, 2)):
        assert _angulo_o1(code, "C", anillo, vecino) < 5.0


def test_o1_reconstruido_no_valida_en_el_anillo_menos_uno_de_gh1():
    """Resultado de D-017: 14,1° en el anillo -1 del 3SCU, más de 10°. Por eso GH5 conserva el marco de 2020."""
    assert 10.0 < _angulo_o1("3SCU", "C", 5, 4) < 20.0
