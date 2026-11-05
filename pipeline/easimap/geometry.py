"""Geometría de EASiMap: ajuste rígido, sistema de referencia del sustrato y coordenadas esféricas."""
from __future__ import annotations

import numpy as np


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotación R y traslación t que llevan `mobile` sobre `target` con RMSD mínimo.

    Se aplica como `R @ x + t`. Sin pesos, como `measure fit` de VMD.
    """
    cm, ct = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - cm).T @ (target - ct))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, ct - rot @ cm


def substrate_frame(anomeric_c: np.ndarray, glycosidic_o: np.ndarray, c4: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sistema de referencia del contrato de datos.

    Origen en el carbono anomérico del azúcar -1. El eje y apunta al oxígeno glicosídico.
    El eje x apunta a la proyección del C4 en el plano y = 0. Devuelve (origen, matriz con los ejes por filas).
    """
    y = glycosidic_o - anomeric_c
    y = y / np.linalg.norm(y)
    v = c4 - anomeric_c
    v = v - np.dot(v, y) * y
    x = v / np.linalg.norm(v)
    return anomeric_c, np.array([x, y, np.cross(x, y)])


def to_frame(points: np.ndarray, origin: np.ndarray, axes: np.ndarray) -> np.ndarray:
    return (np.atleast_2d(points) - origin) @ axes.T


def spherical(p: np.ndarray) -> tuple[float, float, float]:
    """(distancia, theta, phi) de un punto en el marco del sustrato.

    theta: ángulo desde +y, de 0 a pi. phi: atan2(z, x), de 0 a 2 pi.
    """
    x, y, z = (float(v) for v in p)
    return float(np.linalg.norm(p)), float(np.arctan2(np.hypot(x, z), y)), float(np.arctan2(z, x) % (2 * np.pi))


TETRAHEDRAL_HALF = np.radians(109.4712 / 2)    # mitad del ángulo entre los dos sustituyentes que faltan en el C1
C_O_BOND = 1.43                                # Å, enlace C1-O1 de un glucósido


def beta_o1_direction(c1, c2, o5) -> np.ndarray:
    """Dirección unitaria del O1 del anómero beta de una D-piranosa, desde C1, C2 y O5 (geometría tetraédrica).

    El C1 tiene cuatro sustituyentes: O5, C2, O1 y H1. Los dos que faltan (O1 y H1) están a ±54,74° de la bisectriz
    opuesta a O5 y C2, en el plano perpendicular al plano O5-C1-C2. Se elige el candidato con la quiralidad del beta:
    el triple producto u(O5) . (u(C2) x u(O1)) es positivo en una D-piranosa beta. El signo no depende de la
    conformación del anillo. Se comprobó en 6 enlaces beta reales (3ZMR y 3SCU): valía de 0,60 a 0,78 en todos.
    """
    u2 = (c2 - c1) / np.linalg.norm(c2 - c1)
    u5 = (o5 - c1) / np.linalg.norm(o5 - c1)
    bis = -(u2 + u5)
    bis /= np.linalg.norm(bis)
    n = np.cross(u2, u5)
    n /= np.linalg.norm(n)
    cands = [bis * np.cos(TETRAHEDRAL_HALF) + s * n * np.sin(TETRAHEDRAL_HALF) for s in (1.0, -1.0)]
    return max(cands, key=lambda d: float(np.dot(u5, np.cross(u2, d))))
