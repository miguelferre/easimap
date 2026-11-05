"""Lectura mínima de ficheros PDB: átomos por cadena, agrupados en residuos.

Solo lee lo que EASiMap necesita. No interpreta TER ni CONECT.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

AA3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
    "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
    "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TYR": "Y", "TRP": "W",
}
BACKBONE = ("N", "CA", "C", "O")


@dataclass
class Atom:
    name: str
    altloc: str
    xyz: np.ndarray
    element: str
    hetero: bool


@dataclass
class Residue:
    resname: str
    resseq: int
    icode: str
    atoms: list[Atom] = field(default_factory=list)

    def names(self) -> set[str]:
        return {a.name for a in self.atoms}

    def is_protein(self) -> bool:
        """Como el macro `protein` de VMD: el residuo tiene N, CA, C y O."""
        return set(BACKBONE) <= self.names()


def read_pdb_chain(path, chain: str = "A", model: int = 1) -> list[Residue]:
    """Residuos de una cadena, en el orden del fichero.

    Un residuo nuevo empieza cuando cambia (resseq, icode, resname) entre átomos consecutivos.
    Se leen ATOM y HETATM. Se conservan todos los altloc, como hacía VMD.
    """
    out: list[Residue] = []
    cur: tuple | None = None
    n_model = 1
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            rec = line[:6]
            if rec == "MODEL ":
                n_model = int(line[10:14])
                continue
            if rec == "ENDMDL" and n_model == model:
                break
            if rec not in ("ATOM  ", "HETATM") or n_model != model:
                continue
            if line[21] != chain:
                continue
            key = (int(line[22:26]), line[26].strip(), line[17:20].strip())
            if key != cur:
                out.append(Residue(key[2], key[0], key[1]))
                cur = key
            el = line[76:78].strip() if len(line) >= 78 else ""
            out[-1].atoms.append(Atom(
                line[12:16].strip(), line[16].strip(),
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                el, rec == "HETATM"))
    return out


def read_ligand_atoms(path, resname: str, chain: str = "A", resseq: int | None = None) -> dict[str, np.ndarray]:
    """Coordenadas de los átomos de un ligando, por nombre."""
    res: dict[str, np.ndarray] = {}
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM") or line[21] != chain:
                continue
            if line[17:20].strip() != resname:
                continue
            if resseq is not None and int(line[22:26]) != resseq:
                continue
            res.setdefault(line[12:16].strip(), np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return res


def atom_from_spec(path, spec: dict) -> np.ndarray:
    """Coordenadas de un átomo de ligando. `spec`: chain, resseq, resname, atom y, si hace falta, altloc.

    Sin `altloc`, vale la primera coincidencia. Con `altloc`, solo esa conformación (el marco de GH5 de 2020 usa dos).
    """
    want_alt = spec.get("altloc")
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM") or line[21] != spec["chain"]:
                continue
            if line[17:20].strip() != spec["resname"] or int(line[22:26]) != spec["resseq"]:
                continue
            if line[12:16].strip() != spec["atom"]:
                continue
            if want_alt is not None and line[16].strip() != want_alt:
                continue
            return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    raise LookupError(f"no hay átomo {spec} en {path}")


def frame_atoms(path, substrate: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Los tres átomos del marco: origen, punto sobre el eje y, punto del plano xy (x > 0)."""
    return tuple(atom_from_spec(path, substrate[k]) for k in ("origin", "y_axis", "x_plane"))
