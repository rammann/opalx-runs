#!/usr/bin/env python
"""fmlib.py -- shared code for the G4beamline field map study.

Three jobs:

  1. Write G4beamline field maps in both formats OPALX can read, from fields
     given in closed form: the cartesian ``grid`` format (G4BL3DMagnetoStatic)
     and the axisymmetric ``cylinder`` format (G4BL2DMagnetoStatic).
  2. Provide the closed-form transfer matrices the tracking is compared with.
  3. Read the OPALX output back and turn it into transfer matrices and
     reference-orbit numbers.

Why the fields are the ones they are: a constant field and a field that is
linear in each coordinate are reproduced by trilinear interpolation with no
discretisation error at all, and the same holds for bilinear interpolation of
a constant field in the cylinder format. So for the dipole, quadrupole and
uniform-Bz cases the map carries the analytic field exactly, and any
disagreement with the closed form is a reader, placement or integration
problem rather than a sampling problem. This is the same argument the unit
test unit_tests/AbsBeamline/TestFieldmapVsAnalytic.cpp makes.

Units in the files are millimetres and Tesla, and the values are absolute:
the G4beamline readers do not normalise, so SCALE on FIELDMAP (and KS on
SOLENOID) is a plain multiplier.

Phase space convention, same as the other studies here: OPALX dumps x, y, z in
the beam co-moving frame in metres and px, py, pz in beta*gamma rotated into
the reference direction. We form (x, x'=px/pz, y, y'=py/pz, z, delta).

Uses the shared readers in opalxruns/opalx_diagnostics.py.
Run with the conda python (numpy/h5py/matplotlib). See run_all.sh.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from opalxruns.case import Case as CaseBase
from opalxruns.g4bl import write_cylinder_map, write_grid_map  # noqa: F401  (used as F.write_*)
from opalxruns.matrices import drift_matrix, symplectic_residual, transfer_matrix  # noqa: F401
from opalxruns.results import Results  # noqa: F401  (used as F.Results)

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Beam and geometry. One place, used by the generator and by the tests.
# ---------------------------------------------------------------------------

EMASS = 0.51099895e-3            # GeV, electron rest mass
EDES = 0.1                       # GeV, kinetic energy
CLIGHT_GEV = 0.299792458         # GeV/c per T m

GAMMA = (EDES + EMASS) / EMASS
BETA = math.sqrt(1.0 - 1.0 / GAMMA**2)
BG0 = GAMMA * BETA               # beta*gamma, what parts.txt carries
P0_GEV = BG0 * EMASS
BRHO = P0_GEV / CLIGHT_GEV       # T m
CHARGE = -1.0                    # electron

Z_MAP = 1.0                      # lab z of the map's own z = 0, for most cases
HALF_L = 0.5                     # map runs from -HALF_L to +HALF_L in its own z
L_FIELD = 2.0 * HALF_L           # length of the field region along z
ZSTOP = 2.0                      # path length at which tracking stops
Z_SHIFT = 0.25                   # how far the placement case moves the map

# Field strengths. The dipole is set from the bend angle we want to see.
BEND_DEG = 10.0
BEND_ANGLE = math.radians(BEND_DEG)
RHO_DIP = L_FIELD / math.sin(BEND_ANGLE)
B_DIP = BRHO / RHO_DIP           # T, uniform By
G_QUAD = 0.25 * BRHO             # T/m, so |k1| = 0.25 m^-2
B_SOL = 0.1                      # T, uniform Bz
E_LONG = 1.0                     # MV/m, uniform Ez for the electric cases

DT_FINE = 1e-13                  # s, the step every case uses ...
DT_COARSE = 2e-13                # ... except the step-size twin

EPS = {"x": 1e-4, "xp": 1e-4, "y": 1e-4, "yp": 1e-4, "z": 1e-4, "delta": 1e-3}

COORDS = ["x", "x'", "y", "y'", "z", "delta"]
PAIRS = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12)]  # (plus, minus) rows
PART_LABELS = ["ref", "x+", "x-", "xp+", "xp-", "y+", "y-", "yp+", "yp-",
               "z+", "z-", "delta+", "delta-"]


def load_manifest() -> dict:
    return {m["name"]: m for m in json.loads((HERE / "cases.json").read_text())}


# ---------------------------------------------------------------------------
# The fields, in closed form. x, y, z in metres (the map's own frame),
# returning (Bx, By, Bz) in Tesla.
# ---------------------------------------------------------------------------

def field_dipole(x, y, z):
    """Uniform By. Constant, so trilinear interpolation is exact."""
    return 0.0 * x, B_DIP + 0.0 * x, 0.0 * x


def field_quad(x, y, z, gradient=None):
    """Quadrupole, By = g x and Bx = g y. Linear in one coordinate each, so
    trilinear interpolation is exact."""
    g = G_QUAD if gradient is None else gradient
    return g * y, g * x, 0.0 * x


def field_solenoid(x, y, z):
    """Uniform Bz. Not divergence free at the two ends, which is deliberate:
    the point is a field the map carries exactly and whose motion has a closed
    form, not a realistic magnet."""
    return 0.0 * x, 0.0 * x, B_SOL + 0.0 * x


def field_efield_long(x, y, z):
    """Uniform Ez, in MV/m as the files carry it. Constant, so trilinear interpolation is
    exact. A particle crossing it gains q E dz of energy no matter what path it takes, which
    is what makes it easy to check."""
    return 0.0 * x, 0.0 * x, E_LONG + 0.0 * x


def ramp_bz_on_axis(z):
    """A solenoid profile that is deliberately not symmetric in z, so that
    reading it back to front is a different field and ZREVERSE has something
    to prove. Small at both map ends (about 7e-4 of the peak)."""
    return B_SOL * np.exp(-((z / 0.18) ** 2)) * (1.0 + 0.6 * np.tanh(z / 0.12))


def ramp_dbz_on_axis(z, h=1e-6):
    return (ramp_bz_on_axis(z + h) - ramp_bz_on_axis(z - h)) / (2.0 * h)


def ramp_cylinder(r, z):
    """(Bz, Br) of the asymmetric profile, to the usual first order in r."""
    return ramp_bz_on_axis(z), -0.5 * r * ramp_dbz_on_axis(z)


# ---------------------------------------------------------------------------
# Writing the map files (the writers themselves are in opalxruns/g4bl.py)
# ---------------------------------------------------------------------------

def _grid_axis(start_mm, step_mm, n):
    return [start_mm + i * step_mm for i in range(n)]


# ---------------------------------------------------------------------------
# Closed-form transfer matrices (drift_matrix and symplectic_residual are in
# opalxruns.matrices)
# ---------------------------------------------------------------------------

def _focus_block(k: float, L: float) -> np.ndarray:
    """2x2 for x'' = -k x: trigonometric for k > 0, hyperbolic for k < 0."""
    if k > 0:
        w = math.sqrt(k)
        C, S = math.cos(w * L), math.sin(w * L) / w
    elif k < 0:
        w = math.sqrt(-k)
        C, S = math.cosh(w * L), math.sinh(w * L) / w
    else:
        C, S = 1.0, L
    return np.array([[C, S], [-k * S, C]])


def k1_quad(gradient: float = None) -> float:
    """Horizontal focusing constant. x'' = -k1 x and y'' = +k1 y follow from
    the Lorentz force with By = g x and Bx = g y, so k1 = sign(q) g / (B rho).
    The electron's negative charge makes a positive gradient defocus in x."""
    g = G_QUAD if gradient is None else gradient
    return math.copysign(1.0, CHARGE) * g / BRHO


def quad_matrix(k1: float, L: float) -> np.ndarray:
    """4x4 transverse matrix of a quadrupole with hard ends."""
    M = np.zeros((4, 4))
    M[0:2, 0:2] = _focus_block(k1, L)
    M[2:4, 2:4] = _focus_block(-k1, L)
    return M


def k_solenoid(bz: float = B_SOL) -> float:
    """Rotation rate of the transverse momentum, rad/m: k = sign(q) Bz/(B rho),
    from x'' = k y' and y'' = -k x'."""
    return math.copysign(1.0, CHARGE) * bz / BRHO


def solenoid_matrix(k: float, L: float) -> np.ndarray:
    """4x4 transverse matrix of a region of uniform Bz with no end field.

    Integrating x'' = k y', y'' = -k x' gives a rotation of (x', y') through
    kL, and x and y follow from integrating that once more. This is NOT the
    usual solenoid matrix: a real solenoid's focusing comes from the radial
    field at its ends, and a map holding a uniform Bz has none, so a particle
    entering parallel to the axis simply carries on straight.
    """
    C, s = math.cos(k * L), math.sin(k * L)
    S = s / k if k != 0.0 else L
    D = (1.0 - C) / k if k != 0.0 else 0.0
    return np.array([
        [1.0, S, 0.0, D],
        [0.0, C, 0.0, s],
        [0.0, -D, 1.0, S],
        [0.0, -s, 0.0, C],
    ])


def dipole_exit(bend_radius: float, length: float) -> tuple[float, float]:
    """(exit angle, transverse offset at the exit face) of a particle entering
    a box of uniform By on the axis: a circular arc of the given radius, left
    when it has travelled ``length`` along z."""
    theta = math.asin(length / bend_radius)
    return theta, bend_radius * (1.0 - math.cos(theta))


def energy_gain(e_mvpm: float, length: float, charge: float = CHARGE) -> float:
    """Energy a particle picks up crossing a uniform longitudinal electric field [GeV].

    The work is q E dz integrated along the path, and for a field with no z dependence
    that is just q E times the z it crossed -- so this holds whether or not the particle
    is also being bent sideways by a magnetic field."""
    return charge * e_mvpm * 1e6 * length * 1e-9


# ---------------------------------------------------------------------------
# Reading the OPALX output
# ---------------------------------------------------------------------------

class Case(CaseBase):
    """One case directory: its manifest entry plus the tracking output. The
    shared reading (stat table, dumps, read_plane, trajectory) is in
    opalxruns.case.Case."""

    def __init__(self, name: str, manifest: dict):
        super().__init__(HERE / name, manifest[name], BG0)

    # -- stat file ---------------------------------------------------------
    def field_window(self, frac: float = 1e-3) -> tuple[float, float]:
        """(first s, last s) at which the reference particle sees any field.

        Only usable where the reference orbit is off the magnetic axis or the
        field is non-zero on it -- a quadrupole's axis field is zero, so the
        tests take the window for those cases from the manifest instead."""
        df = self.stat_df
        s = df["s"].to_numpy()
        b = np.sqrt(df["Bx_ref"].to_numpy()**2 + df["By_ref"].to_numpy()**2
                    + df["Bz_ref"].to_numpy()**2)
        nz = np.where(b > frac * b.max())[0]
        if len(nz) == 0:
            raise RuntimeError(f"{self.name}: the reference particle saw no field")
        return float(s[nz[0]]), float(s[nz[-1]])

    def ref_energy_change(self) -> float:
        """Total energy the reference particle gained over the run [GeV]."""
        df = self.stat_df
        bg = np.sqrt(df["ref_px"].to_numpy()**2 + df["ref_py"].to_numpy()**2
                     + df["ref_pz"].to_numpy()**2)
        W = np.sqrt((bg * EMASS)**2 + EMASS**2)
        return float(W[-1] - W[0])

    def exit_angle(self) -> float:
        """Angle of the reference momentum at the last stat row, rad."""
        _, _, _, px, pz = self.ref_orbit()
        return math.atan2(px[-1], pz[-1])

    def ref_x_at(self, s_target: float) -> float:
        s, x, _, _, _ = self.ref_orbit()
        return float(np.interp(s_target, s, x))

    def ref_x_at_z(self, z_target: float) -> float:
        """Reference orbit x where it reaches a lab z.

        Not the same as reading it at that path length: a bent orbit is longer
        than its extent along z, so on the dipole case the two differ by 5 mm
        at the exit face."""
        _, x, z, _, _ = self.ref_orbit()
        return float(np.interp(z_target, z, x))

    # -- particle dumps ----------------------------------------------------
    def planes(self, margin: float = 0.03):
        """The last dump before the field and the first one after it."""
        steps, sp = self._steps_spos()
        s0, s1 = self.m["field_s0"], self.m["field_s1"]
        before = np.where(sp < s0 - margin)[0]
        after = np.where(sp > s1 + margin)[0]
        if len(before) == 0 or len(after) == 0:
            raise RuntimeError(f"{self.name}: no field-free dump on one side")
        return steps, sp, int(before[-1]), int(after[0])

    def transfer_map(self, margin: float = 0.03):
        """6x6 map of the field region by centred finite differences between a
        dump before it and one after, with the known drifts stripped off."""
        steps, sp, ie, ix = self.planes(margin)
        s0, s1 = self.m["field_s0"], self.m["field_s1"]
        Xin, Xout = self.read_plane(steps[ie]), self.read_plane(steps[ix])
        M_planes = transfer_matrix(Xin, Xout, PAIRS)
        M = (np.linalg.inv(drift_matrix(sp[ix] - s1, GAMMA)) @ M_planes
             @ np.linalg.inv(drift_matrix(s0 - sp[ie], GAMMA)))
        return M, {"s_in": float(sp[ie]), "s_out": float(sp[ix])}

    def transfer_map4(self, margin: float = 0.03) -> np.ndarray:
        M, _ = self.transfer_map(margin)
        return M[0:4, 0:4]

    def final_state(self) -> np.ndarray:
        steps, _ = self._steps_spos()
        return self.read_plane(steps[-1])

    def momentum_error(self) -> float:
        """Largest relative change in |p| between the first and last dump.
        A magnetic field does no work, so this should be at round-off."""
        steps, _ = self._steps_spos()
        return float((np.abs(self._pmag(steps[-1]) - self._pmag(steps[0])) / BG0).max())
