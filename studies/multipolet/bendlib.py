#!/usr/bin/env python
"""bendlib.py -- analysis for the MULTIPOLET bend study.

Adapted from runs/bendtest/validation/bendlib.py. Reads OPALX tracking output
and turns it into the quantities the tests compare against analytic optics:
the 6x6 transfer map, dispersion, the reference-orbit bend angle, the on-axis
field integral, and bunch covariance transport.

Phase-space convention (verified against the tracker output):
  * OPALX (opal-t) dumps per-particle x,y,z in the beam co-moving frame [m] and
    px,py,pz in beta*gamma, rotated into the reference direction.
  * We form (x, x'=px/pz, y, y'=py/pz, z, delta=|p|/p0-1).
  * Different particles in one dump share the same time, not the same s. To
    compare against an s-based transfer matrix we (a) place the read planes in
    the field-free drifts on either side of the bend, (b) project each particle
    onto the reference transverse plane (drift by -z), and (c) strip the known
    drift from each plane to the design face. The bend map is then built by
    centered finite differences from the +/- perturbation pairs.

Sign convention of the analytic map: the map is measured in the co-moving
frame, whose x axis is the entrance x carried along the orbit. Both SBEND and
the curved MULTIPOLET turn the orbit towards -x, so that x points away from the
centre of curvature and the plain sector matrix applies (``x_sign`` = +1). The
conjugation with diag(x_sign, x_sign, 1, 1, 1, 1) is kept so a element with the
other handedness could still be described.

Uses the shared readers in opalxruns/opalx_diagnostics.py.
Run with the conda python (numpy/scipy/h5py/matplotlib). See run_all.sh.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from opalxruns import results
from opalxruns.case import Case as CaseBase
from opalxruns.matrices import (  # noqa: F401  (used as bl.<name>)
    centred_differences, drift_matrix, symplectic_residual, transfer_matrix,
)

HERE = Path(__file__).resolve().parent

COORDS = ["x", "x'", "y", "y'", "z", "delta"]
PAIRS = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12)]  # (plus_idx, minus_idx)

C_LIGHT = 0.299792458  # GeV/c per T m


def load_manifest() -> dict:
    return {m["name"]: m for m in json.loads((HERE / "cases.json").read_text())}


def brho(P0_GeV: float) -> float:
    """Magnetic rigidity B*rho [T m] for a momentum in GeV/c."""
    return P0_GeV / C_LIGHT


# ---------------------------------------------------------------------------
# Analytic matrices
# ---------------------------------------------------------------------------

def _cs(k: float, L: float) -> tuple[float, float]:
    """(C, S) of the focusing block for constant k: cos/sin, cosh/sinh, or drift."""
    if k > 0:
        w = math.sqrt(k)
        return math.cos(w * L), math.sin(w * L) / w
    if k < 0:
        w = math.sqrt(-k)
        return math.cosh(w * L), math.sinh(w * L) / w
    return 1.0, L


def cf_sector_matrix(h: float, k1: float, L: float, gamma: float) -> np.ndarray:
    """6x6 combined-function sector bend, x pointing away from the centre of
    curvature. h = 1/rho, k1 = (q/p) dBy/dx (k1 > 0 focuses in x). kx = h^2 + k1,
    ky = -k1. Longitudinal row signs follow the OPALX co-moving convention
    (R51 = -h S, R52 = -R16, R56 = -h^2 (L - S)/kx + L/gamma^2); the matrix is
    symplectic. With k1 = 0 this is the plain sector matrix."""
    kx, ky = h * h + k1, -k1
    Cx, Sx = _cs(kx, L)
    Cy, Sy = _cs(ky, L)
    M = np.eye(6)
    M[0, 0], M[0, 1] = Cx, Sx
    M[1, 0], M[1, 1] = -kx * Sx, Cx
    M[2, 2], M[2, 3] = Cy, Sy
    M[3, 2], M[3, 3] = -ky * Sy, Cy
    if kx != 0.0:
        D, Dp = h * (1.0 - Cx) / kx, h * Sx
        R56 = -h * h * (L - Sx) / kx
    else:
        D, Dp = h * L * L / 2.0, h * L
        R56 = -h * h * L**3 / 6.0
    M[0, 5], M[1, 5] = D, Dp
    M[4, 0], M[4, 1] = -Dp, -D
    M[4, 5] = R56 + L / gamma**2
    return M


def k1_of(m: dict) -> float:
    """Horizontal focusing constant of the gradient TP[1]: k1 = sign(q) B1 / Brho.

    TP[1] is the physical gradient dBy/dx in the co-moving frame, so this is the
    ordinary rule: a positive gradient defocuses a negative charge in x, the same
    sign convention as MULTIPOLE / QUADRUPOLE K1 (see the multipole study)."""
    if m["element"] == "MULTIPOLET" and m["B1"]:
        return math.copysign(1.0, m["charge"]) * m["B1"] / brho(m["P0_GeV"])
    return 0.0


def analytic_map(m: dict) -> np.ndarray:
    """Ideal map for one manifest entry: sector bend (plus the gradient k1 for
    the combined-function MULTIPOLET cases), Maxwell-consistent (ky = -k1), in
    the co-moving x convention of the element (x_sign)."""
    M = cf_sector_matrix(1.0 / m["rho"], k1_of(m), m["arc"], m["gamma"])
    F = np.diag([m["x_sign"], m["x_sign"], 1.0, 1.0, 1.0, 1.0])
    return F @ M @ F


# ---------------------------------------------------------------------------
# Symplecticity helpers
# ---------------------------------------------------------------------------

def block_dets(M: np.ndarray) -> list[float]:
    return [float(np.linalg.det(M[i:i + 2, i:i + 2])) for i in (0, 2, 4)]


# ---------------------------------------------------------------------------
# One tracked case
# ---------------------------------------------------------------------------

class Case(CaseBase):
    """One case directory: manifest metadata + output files. The shared reading
    (stat table, dumps, read_plane, trajectory) is in opalxruns.case.Case."""

    def __init__(self, name: str, manifest: dict):
        super().__init__(HERE / name, manifest[name], manifest[name]["bg0"])
        self.gamma = self.m["gamma"]

    def field_extent(self, frac: float = 0.01) -> tuple[float, float]:
        """(s_start, s_end) where |By_ref| exceeds ``frac`` of its peak."""
        df = self.stat_df
        s = df["s"].to_numpy()
        by = np.abs(df["By_ref"].to_numpy())
        nz = np.where(by > frac * by.max())[0]
        if len(nz) == 0:
            raise RuntimeError(f"{self.name}: By_ref is zero everywhere (reference "
                               "particle saw no field)")
        return float(s[nz[0]]), float(s[nz[-1]])

    def field_integral(self) -> float:
        """int By ds over the field region [T m]."""
        df = self.stat_df
        s = df["s"].to_numpy()
        by = df["By_ref"].to_numpy()
        s0, s1 = self.field_extent()
        sel = (s >= s0) & (s <= s1)
        return float(np.trapezoid(by[sel], s[sel]))

    def field_profile(self):
        """(s, By_ref / B0) along the line."""
        df = self.stat_df
        return df["s"].to_numpy(), df["By_ref"].to_numpy() / self.m["B0"]

    def planes(self, margin: float):
        steps, sp = self._steps_spos()
        s0, s1 = self.field_extent()
        before = np.where(sp < s0 - margin)[0]
        after = np.where(sp > s1 + margin)[0]
        if len(before) == 0 or len(after) == 0:
            raise RuntimeError(f"{self.name}: no field-free plane (margin too big?)")
        return steps, sp, before[-1], after[0], (s0, s1)

    def transfer_map(self, margin: float = 0.05):
        """6x6 bend map by centered finite differences between a plane in the
        entrance drift and one in the exit drift, drifts stripped to the design
        faces. Returns (M, info)."""
        steps, sp, ie, ix, field = self.planes(margin)
        face_in, face_out = self.m["face_in_s"], self.m["face_out_s"]
        Xin = self.read_plane(steps[ie])
        Xout = self.read_plane(steps[ix])
        Mplanes = transfer_matrix(Xin, Xout, PAIRS)
        Din = drift_matrix(face_in - sp[ie], self.gamma)
        Dout = drift_matrix(sp[ix] - face_out, self.gamma)
        M = np.linalg.inv(Dout) @ Mplanes @ np.linalg.inv(Din)
        info = {"s_in": float(sp[ie]), "s_out": float(sp[ix]), "field": field,
                "eps_in": np.diag(centred_differences(Xin, PAIRS)).copy()}
        return M, info

    def sigma_at_faces(self, margin: float = 0.05):
        """(Sigma_in, Sigma_out, mean_in, mean_out) at the design faces."""
        steps, sp, ie, ix, _ = self.planes(margin)
        face_in, face_out = self.m["face_in_s"], self.m["face_out_s"]
        Xi = drift_matrix(face_in - sp[ie], self.gamma) @ self.read_plane(steps[ie])
        Xo = np.linalg.inv(drift_matrix(sp[ix] - face_out, self.gamma)) @ self.read_plane(steps[ix])
        return np.cov(Xi), np.cov(Xo), Xi.mean(1), Xo.mean(1)

    def momentum_conservation(self) -> dict:
        steps, _ = self._steps_spos()
        p0 = self._pmag(steps[0])
        pN = self._pmag(steps[-1])
        rel = np.abs(pN - p0) / self.bg0
        return {"max_rel": float(rel.max()), "per_particle": rel}

    def survived(self) -> tuple[int, int]:
        """(particles at the first dump, particles at the last dump)."""
        df = self.stat_df
        n = df["numParticles"].to_numpy()
        return int(n[0]), int(n[-1])


def exit_bend_angle(case: Case) -> float:
    """Realized bend angle from the reference momentum at the last stat row."""
    _, _, _, px, pz = case.ref_orbit()
    return math.atan2(abs(px[-1]), pz[-1])


def bend_direction(case: Case) -> float:
    """Sign of the lab-frame x the reference orbit turns towards."""
    _, x, _, _, _ = case.ref_orbit()
    return math.copysign(1.0, x[-1])


# ---------------------------------------------------------------------------
# Covariance transport (bunch tests)
# ---------------------------------------------------------------------------

def emittance(sig: np.ndarray, plane: int) -> float:
    i = 2 * plane
    return math.sqrt(max(sig[i, i] * sig[i + 1, i + 1] - sig[i, i + 1] ** 2, 0.0))


def emittance_6d(sig: np.ndarray) -> float:
    return math.sqrt(max(np.linalg.det(sig), 0.0))


def centroid_orbit_shift(off: Case, design: Case, s_target: float) -> float:
    """Perpendicular separation of the ``off`` reference orbit from the ``design``
    one at path length s_target: the off-energy bunch centroid's x in the design
    co-moving frame (OPALX's reference particle rides with its own bunch)."""
    sd, xd, zd, pxd, pzd = design.ref_orbit()
    sm, xm, zm, _, _ = off.ref_orbit()
    idd = int(np.argmin(np.abs(sd - s_target)))
    idm = int(np.argmin(np.abs(sm - s_target)))
    dx, dz = xm[idm] - xd[idd], zm[idm] - zd[idd]
    p = math.hypot(pxd[idd], pzd[idd])
    tx, tz = pxd[idd] / p, pzd[idd] / p
    return dx * tz + dz * (-tx)


# ---------------------------------------------------------------------------
# Result table / tolerance checking
# ---------------------------------------------------------------------------

class Results(results.Results):
    """The shared result table, with the "analytic" column this study has always printed."""

    def __init__(self):
        super().__init__(name_width=30, expected="analytic", rule=100)
