"""Réplica del procedimiento de 2020 con Python.

Reproduce lo que hacían los scripts Tcl de VMD, con estas decisiones fijas:

* Solo la cadena A, modelo 1. Se conservan todos los altloc (VMD los cargaba como átomos distintos).
* Superposición sobre N, CA, C y O de los residuos de referencia. 2020 pedía también CB,
  pero `backbone` de VMD no incluye CB; el ajuste real usaba ocho átomos.
* Marco del sustrato en el 3SCU: origen en C1 del anillo -1, y hacia O4 del anillo +1.
* Punto medido: centro de los átomos de la cadena lateral más el CA, sin pesos.
* Una fila por resid >= 0 de la cadena A. Se ordenan como texto, como hacía `lsort` de Tcl.
* Una proteína se omite si sus átomos de ajuste no suman los mismos que los de la referencia.

Diferencia conocida: el HMM de 2020 no se guardó. Aquí se usa PF00232 actual, con el número real de
columna. Los residuos de referencia se fijan por resid en el 3SCU, no por número de columna.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import collections
import numpy as np

from .geometry import kabsch, spherical, substrate_frame, to_frame
from .hmmalign import assign_positions, assign_positions_local, load_hmm
from .structure import AA3TO1, BACKBONE, Residue, frame_atoms, read_pdb_chain


@dataclass
class ChainTable:
    """Una fila por resid de una proteína, en el orden de 2020."""
    pdb: str
    resid: list[int] = field(default_factory=list)
    resname: list[str] = field(default_factory=list)
    hmm_position: list[int] = field(default_factory=list)   # 0 si no tiene o si el resid es ambiguo
    distance: list[float] = field(default_factory=list)
    theta: list[float] = field(default_factory=list)
    phi: list[float] = field(default_factory=list)
    xyz: list[np.ndarray] = field(default_factory=list)


@dataclass
class Skipped:
    pdb: str
    reason: str


def protein_residues(path) -> list[Residue]:
    return [r for r in read_pdb_chain(path, "A") if r.is_protein()]


def number_residues(hmm, residues: list[Residue], local: bool = True) -> list[int]:
    """Posición HMM de cada residuo. Lo que no es de los 20 aminoácidos estándar entra como Z, como en 2020.

    local=True alinea en modo local (como 2020). local=False alinea en modo glocal, que fuerza todo el modelo.
    """
    seq = "".join(AA3TO1.get(r.resname, "Z") for r in residues)
    if local:
        return assign_positions_local(hmm, seq)
    return assign_positions(hmm, seq, all_consensus_cols=True)


def fit_atoms(residues: list[Residue], positions: list[int], cols: list[int], altloc: str | None = None) -> list | None:
    """Átomos N, CA, C, O de los residuos en las columnas dadas, en orden de fichero. None si falta alguno.

    Por defecto entran todos los altloc, como los cargaba VMD (un residuo con dos conformaciones aporta 8 átomos).
    altloc="first" o "last" deja un solo átomo por nombre en cada residuo (el primero o el último del fichero).
    """
    out = []
    for c in cols:
        hit = [r for r, p in zip(residues, positions) if p == c]
        if not hit:
            return None
        for r in hit:
            atoms = [a for a in r.atoms if a.name in BACKBONE]
            if altloc:
                by: dict[str, object] = {}
                for a in atoms:
                    if altloc == "last" or a.name not in by:
                        by[a.name] = a
                atoms = [a for a in atoms if by[a.name] is a]
            out += atoms
    return out


# Átomos que el macro `backbone` de VMD trata como esqueleto: N, CA, C, O, los oxígenos del extremo C,
# el hidrógeno amida, los del extremo N y los HA del carbono alfa (medido: con HA fuera, 5YJ7 y 5NS6 pasan de 0,12 Å a 0,00007 Å).
VMD_BACKBONE = {"N", "CA", "C", "O", "OT1", "OT2", "O1", "O2", "OXT", "H", "HN", "HT1", "HT2", "HT3", "H1", "H2", "H3",
                "HA", "HA2", "HA3"}


def sidechain_point(atoms) -> np.ndarray:
    """Centro de los átomos de la cadena lateral (`sidechain` de VMD) más el CA. Un CA cuenta una vez por altloc."""
    sel = [a.xyz for a in atoms if a.name not in VMD_BACKBONE or a.name == "CA"]
    return np.mean(sel, axis=0)


def substrate_from_reference(path, cfg: dict):
    """Marco del sustrato a partir de los átomos del ligando en la estructura de referencia."""
    return substrate_frame(*frame_atoms(path, cfg["substrate"]))


def reference_columns(hmm, ref_path, resids: list[int], local: bool = True) -> list[int]:
    residues = protein_residues(ref_path)
    pos = number_residues(hmm, residues, local)
    cols = []
    for rid in resids:
        # "186A" = resid 186 con código de inserción A (dos residuos pueden compartir el número)
        num, ic = (int(str(rid)[:-1]), str(rid)[-1]) if isinstance(rid, str) else (rid, None)
        hit = [p for r, p in zip(residues, pos) if r.resseq == num and (ic is None or r.icode == ic) and p > 0]
        if not hit:
            raise ValueError(f"el resid {rid} de la referencia no cae en una columna del HMM")
        cols.append(hit[0])
    return cols


def tabulate(pdb: str, residues, positions, rot, trans, origin, axes) -> ChainTable:
    """Tabla de 2020 para una proteína ya ajustada a la referencia."""
    by: dict[int, list] = collections.OrderedDict()
    pos_by: dict[int, set] = {}
    names: dict[int, str] = {}
    for r, p in sorted(zip(residues, positions), key=lambda rp: str(rp[0].resseq)):
        if r.resseq < 0:
            continue
        by.setdefault(r.resseq, []).extend(r.atoms)
        pos_by.setdefault(r.resseq, set()).add(p)
        names.setdefault(r.resseq, r.resname)
    t = ChainTable(pdb)
    for rid, atoms in by.items():
        pt = to_frame(rot @ sidechain_point(atoms) + trans, origin, axes)[0]
        d, th, ph = spherical(pt)
        ps = pos_by[rid]
        t.resid.append(rid)
        t.resname.append(names[rid])
        t.hmm_position.append(next(iter(ps)) if len(ps) == 1 else 0)
        t.distance.append(d); t.theta.append(th); t.phi.append(ph); t.xyz.append(pt)
    return t


def positions_from_table(residues: list[Residue], rows: list[dict]) -> list[int]:
    """Posición HMM de cada residuo según la propia tabla de 2020 (sin HMM). Un resid ambiguo o fuera de dominio vale 0.

    Las filas van en el orden de `lsort` de los resid (ver `table2020`). Sirve para comprobar el marco, la superposición
    y el punto medido cuando no se tiene el modelo con el que se numeró en 2020.
    """
    by: dict[int, None] = {}
    for r in residues:
        if r.resseq >= 0:
            by.setdefault(r.resseq, None)
    order = sorted(by, key=str)
    n = len(order)
    if len(rows) not in (n, 2 * n):
        raise ValueError("la tabla de 2020 no tiene una fila por resid")
    pos_of = {order[i]: (rows[i]["position"] or 0) for i in range(n)}
    seen: dict[int, int] = {}
    for r in residues:                     # dos residuos con el mismo resid y otro código de inserción: ambiguo
        seen[r.resseq] = seen.get(r.resseq, 0) + 1
    return [0 if seen.get(r.resseq, 0) > 1 else pos_of.get(r.resseq, 0) for r in residues]


def run(family_cfg: dict, pdb_dir: Path, hmm_path: Path, codes: list[str] | None = None, local: bool = True,
        use_resids: bool = False, fit_altloc: str | None = None, table_rows: dict | None = None):
    """Procesa la lista de 2020. Devuelve (tablas por código, omitidas, columnas de referencia).

    use_resids=False usa las columnas de referencia de 2020 tal cual (GH1: 165 y 360). Solo valen con el HMM de 2020.
    use_resids=True busca las columnas de los residuos de referencia (`reference_resids`) en la estructura de
    referencia con el HMM dado. Es lo que hay que hacer con un HMM vigente, que tiene otra numeración.
    Admite cualquier número de residuos de referencia.
    """
    hmm = load_hmm(hmm_path)
    ref = family_cfg["reference_pdb"]
    rcfg = family_cfg["replica2020"]
    ref_path = pdb_dir / f"{ref}.pdb"
    origin, axes = substrate_from_reference(ref_path, family_cfg)
    if table_rows:
        cols = list(rcfg["reference_positions_2020"])
    elif use_resids:
        cols = reference_columns(hmm, ref_path, rcfg["reference_resids"], local)
    else:
        cols = list(rcfg["reference_positions_2020"])
    ref_res = protein_residues(ref_path)
    ref_pos = positions_from_table(ref_res, table_rows[ref]) if table_rows else number_residues(hmm, ref_res, local)
    fit_altloc = fit_altloc or rcfg.get("fit_altloc")
    ref_atoms = fit_atoms(ref_res, ref_pos, cols, fit_altloc)
    if ref_atoms is None:
        raise ValueError(f"la referencia {ref} no tiene un residuo en alguna de las columnas {cols}")
    ref_xyz = np.array([a.xyz for a in ref_atoms])

    tables: dict[str, ChainTable] = {}
    skipped: list[Skipped] = []
    for code in codes or list(rcfg["proteins"]):
        res = protein_residues(pdb_dir / f"{code}.pdb")
        pos = positions_from_table(res, table_rows[code]) if table_rows else number_residues(hmm, res, local)
        if code == ref:
            rot, trans = np.eye(3), np.zeros(3)
        else:
            atoms = fit_atoms(res, pos, cols, fit_altloc)
            if atoms is None:
                skipped.append(Skipped(code, "falta un residuo de referencia en las columnas del HMM"))
                continue
            if len(atoms) != len(ref_atoms):
                skipped.append(Skipped(code, f"átomos de ajuste distintos: {len(atoms)} frente a {len(ref_atoms)}"))
                continue
            rot, trans = kabsch(np.array([a.xyz for a in atoms]), ref_xyz)
        tables[code] = tabulate(code, res, pos, rot, trans, origin, axes)
    return tables, skipped, cols
