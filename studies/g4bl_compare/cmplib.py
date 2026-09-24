"""Shared pieces for the g4bl_compare study.

Readers for both codes' output, the reduction of each to one set of six
coordinates, the transfer matrix, the tolerance table, and the result table.

Everything here is code-agnostic: it takes a file from G4beamline or from OPALX
and returns the same arrays, so run_tests.py never has to know which code wrote
what.

Units on the way in:
  G4beamline #BLTrackFile   mm, MeV/c, ns
  G4beamline fieldntuple    mm, ns, T
  OPALX monitor .h5         m, beta*gamma, s
  OPALX DUMPEMFIELDS        m, ns, T, MV/m
Units on the way out: m, radians, T, and beta*gamma where a momentum is wanted raw.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from opalxruns.h5 import read_monitor  # noqa: F401  (used as cmplib.read_monitor)
from opalxruns.mue4_beam import (  # noqa: F401  (used as cmplib.<name>)
    BETA0, BG0, BLTRACK_COLS, BRHO, C_MM_PER_NS, E0, MUON_MASS, P0,
    canonical, matched, read_bltrack,
)
from opalxruns.paths import G4BL_FILES

HERE = Path(__file__).resolve().parent
MAPS = G4BL_FILES / "muE4" / "maps"

# The +/- pairs of the 19-particle set, as indices into the id-sorted array.
# Twelve small steps, one pair per coordinate, used for the transfer matrix.
MATRIX_PAIRS = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12)]
# Six large steps in x, y and dp/p, used as the amplitude probe. They are NOT
# part of the matrix: at 20 mm the map is no longer linear, which is the point.
LARGE_PAIRS = [(13, 14), (15, 16), (17, 18)]
STEP_SMALL = 1.0e-3                         # m, rad, and dp/p
STEP_LARGE = 2.0e-2                         # m and dp/p

COORD_NAMES = ["x", "x'", "y", "y'", "zeta", "delta"]

# ---------------------------------------------------------------------------
# Tolerances. Every number has a reason, and every reason was checked against
# the measurement floor first: a comparison can sit inside its tolerance while
# measuring nothing but rounding.
# ---------------------------------------------------------------------------
TOL = {
    # The field tolerance is NOT a fixed number of Tesla, because the floor is
    # not: G4beamline writes 6 significant figures, so a 0.086 T field is
    # resolved to 1e-7 T and a 0.215 T field only to 1e-6 T. A fixed 2e-7 passes
    # the dipoles and fails the quadrupole while both agree perfectly. So the
    # tolerance is three printing steps of that case's own peak field, computed
    # by printing_step(), and these two entries are the multiplier.
    "field_grid_point": 3.0,                # x the printing step at the peak field
    "field_between": 3.0,
    # A map whose rows are all zero must give exactly zero, not nearly zero.
    "field_zero": 1e-12,
    # Nothing has happened yet at the entrance plane. Anything above round-off
    # means the two codes were not handed the same beam.
    "entrance": 1e-9,                       # m
    # The measurement floor is %.6g on a ~1 m coordinate, 1e-6 m, so this is a
    # real bound. The ported asr61_dipole case measures 4.9e-5 m over a 41 degree bend.
    "exit_pos": 1e-4,                       # m
    # x' and y' are momentum ratios, so they carry no frame convention. The floor
    # is one printing step of 28 MeV/c, which is 1e-4 MeV/c and NOT 1e-5 -- six
    # significant figures of 28 is 28.0000. So the floor on a momentum ratio is
    # about 3.6e-6 rad, and 1e-5 is under three floors rather than the twenty-five
    # an earlier comment here claimed. A tighter bound is not measurable from
    # G4beamline's text output.
    "exit_angle": 1e-5,                     # rad
    # A static magnetic field does no work. Checked inside each code as well as
    # between them, so a broken integrator shows up even if both codes break
    # the same way.
    "mom_conserved": 5e-6,                  # relative
    # The transfer matrix is limited by how precisely G4beamline prints its
    # output, not by either code's physics. Six significant figures on the 1.27 m
    # the bent orbit reaches is a 1e-5 m step; divided by the 1e-3 m difference
    # step that is 1e-2 in matrix units. So 2e-2 is two floors, and a matrix
    # cannot be compared better than this without more output precision from
    # G4beamline. printing_floor() computes the floor per case so the table can
    # show it next to the measurement.
    # The primary matrix test is absolute, in units of that floor: the two
    # matrices must agree to within about one printing step. A relative test
    # cannot be the primary one, because dividing a difference that IS the floor
    # by a small matrix element turns it into a large number with no meaning --
    # M[zeta<-x] of the ASR61 group is 3.6% of the largest element, so the floor
    # alone shows up there as 3.4%.
    "matrix_floors": 3.0,                   # multiples of that row's own floor
    # The relative test is kept as a second check, but only on elements at least
    # a hundred times the floor, where the floor can contribute at most 1%.
    "matrix_rel": 2e-2,
    "matrix_big": 100.0,                    # x the floor, to qualify for the relative test
    "gauss_mean": 5e-6,                     # m
    "gauss_rms": 1e-3,                      # relative
    # The 20000-particle stage is tested on the median and the 99th percentile
    # of the per-particle difference, not on its maximum. The difference has a
    # distribution, and the largest of 20000 draws reaches further into its tail
    # than the largest of 19 -- measured on the ASR61 dipole, median 8.0e-6 rad
    # and worst 3.7e-5 rad, where the 19-particle stage's worst is 8.2e-6. So a
    # maximum tested against a tolerance set from 19 particles fails on sample
    # size alone. Halving the time step moves the worst only from 3.7e-5 to
    # 3.1e-5, so it is the tail rather than the step.
    "gauss_median_angle": 1e-5,             # rad, same as the pair stage
    "gauss_p99_angle": 3e-5,                # rad, three times the median tolerance
    "gauss_median_pos": 1e-4,               # m, same as the pair stage
    "gauss_p99_pos": 3e-4,                  # m
}


# ---------------------------------------------------------------------------
# G4beamline readers
# ---------------------------------------------------------------------------
def read_g4bl_field(path) -> dict[str, np.ndarray]:
    """A G4beamline fieldntuple ASCII dump: x y z t Bx By Bz Ex Ey Ez.

    Positions come back in metres; the field is already in Tesla.
    """
    a = np.loadtxt(path, comments="#")
    if a.ndim == 1:
        a = a[None, :]
    return {"x": a[:, 0] * 1e-3, "y": a[:, 1] * 1e-3, "z": a[:, 2] * 1e-3,
            "Bx": a[:, 4], "By": a[:, 5], "Bz": a[:, 6],
            "Ex": a[:, 7], "Ey": a[:, 8], "Ez": a[:, 9]}


# ---------------------------------------------------------------------------
# OPALX readers
# ---------------------------------------------------------------------------
def read_opalx_field(path) -> dict[str, np.ndarray]:
    """An OPALX DUMPEMFIELDS dump: x y z t Bx By Bz Ex Ey Ez, already in m and T.

    The header is a count, then one numbered line per column, then a lone 0.
    Note that OPALX writes this into data/, not the working directory.
    """
    lines = Path(path).read_text().splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.strip() == "0") + 1
    a = np.array([[float(v) for v in ln.split()] for ln in lines[start:] if ln.strip()])
    return {"x": a[:, 0], "y": a[:, 1], "z": a[:, 2],
            "Bx": a[:, 4], "By": a[:, 5], "Bz": a[:, 6],
            "Ex": a[:, 7], "Ey": a[:, 8], "Ez": a[:, 9]}


def read_stat(path):
    """The OPALX .stat file as a DataFrame. ref_x/ref_y/ref_z are the reference
    orbit in the lab frame and Bx_ref/By_ref/Bz_ref the field it sees.

    parse_opal_stat returns (metadata, DataFrame); only the frame is wanted here.
    """
    from opalxruns.opalx_diagnostics import parse_opal_stat
    _meta, df = parse_opal_stat(str(path))
    return df


# ---------------------------------------------------------------------------
# One set of six coordinates from either code
# ---------------------------------------------------------------------------
def time_of_flight(d: dict[str, np.ndarray], code: str) -> float:
    """Arrival time of particle 0, in seconds. Compared between codes on its
    own, because canonical() deliberately subtracts it out."""
    return float(d["t"][0] if code == "g4bl" else d["time"][0])


def duplicates(d: dict) -> list[int]:
    """Ids recorded more than once at one plane.

    A plane inside the field can be crossed twice by the same particle, and the
    monitor records both crossings -- at the ASR61 map centre two of 20000 do it.
    That is the monitor behaving correctly, but it means a plane can hold more
    rows than there are particles, so it is reported rather than absorbed.
    """
    u, n = np.unique(d["id"], return_counts=True)
    return u[n > 1].tolist()


# ---------------------------------------------------------------------------
# Transfer matrix
# ---------------------------------------------------------------------------
def transfer_matrix(Xin: np.ndarray, Xout: np.ndarray,
                    pairs=MATRIX_PAIRS) -> np.ndarray:
    """6 x 6 by centred differences over the +/- pairs, entrance plane to exit.

    Built as (d out) @ inv(d in) rather than assuming the entrance differences
    are exactly the nominal steps, so a plane that sits a little downstream of
    where the particles were made does not bias it.
    """
    din = np.zeros((6, 6))
    dout = np.zeros((6, 6))
    for j, (ip, im) in enumerate(pairs):
        din[:, j] = (Xin[:, ip] - Xin[:, im]) / 2.0
        dout[:, j] = (Xout[:, ip] - Xout[:, im]) / 2.0
    return dout @ np.linalg.inv(din)


def in_field(z_mm: float, boxes, x: float = 0.0, y: float = 0.0) -> bool:
    """Whether a recording plane on the axis sits inside any map's box.

    A plane in the field and a plane in free space do not measure the same
    thing, so they cannot share a tolerance. Position is interpolated onto the
    plane by both codes, but the momentum is reported from the tracking step the
    particle happened to be on, and inside a field the momentum is changing with
    z -- so a fraction of a step of sampling offset turns into a real-looking
    angle difference. Outside the field there is nothing to change and the two
    codes agree to round-off.
    """
    z = z_mm / 1000.0
    for b in boxes:
        if (b[0] <= x <= b[1]) and (b[2] <= y <= b[3]) and (b[4] <= z <= b[5]):
            return True
    return False


def offset_explaining(dangle: float, b_tesla: float) -> float:
    """The longitudinal sampling offset that would produce this angle difference.

    In a field B the angle turns at B / (B*rho) per metre, so an angle difference
    of `dangle` is `dangle * BRHO / B` of z. Turning the number into a length is
    what makes it readable: 2.4e-4 rad in the ASR61 dipole is 0.26 mm, about one
    and a half OPALX steps, rather than a disagreement about the field.
    """
    if b_tesla <= 0:
        return float("nan")
    return dangle * BRHO / b_tesla


def field_at(stat_path, z_m: float) -> float:
    """|B| on the reference orbit nearest this lab z, from the .stat file."""
    df = read_stat(stat_path)
    i = int(np.argmin(np.abs(np.asarray(df["ref_z"]) - z_m)))
    return float(np.hypot(np.hypot(df["Bx_ref"][i], df["By_ref"][i]), df["Bz_ref"][i]))


def printing_step(value: float) -> float:
    """The smallest difference G4beamline's %.6g output can show at this size."""
    v = abs(float(value))
    if v <= 0:
        return 0.0
    return 10.0 ** (np.floor(np.log10(v)) - 5)


def on_face(pos: np.ndarray, boxes, tol: float = 1e-7) -> np.ndarray:
    """Which sample points sit on the bounding face of any map placed in the case.

    Both codes drop a point that is outside a map's box and interpolate one that
    is inside, and on the face itself they disagree about the last bit of the
    float: a sample written as -0.110 and a box edge computed as -110 * 1e-3 are
    not the same double. The difference that produces is the whole field, so a
    handful of face points read as a catastrophic disagreement when nothing is
    wrong. No particle can be affected -- a face has no thickness -- so these
    points are reported separately instead of being tested.

    `boxes` is a list of (xmin, xmax, ymin, ymax, zmin, zmax) in lab metres.
    """
    mask = np.zeros(pos.shape[1], dtype=bool)
    for b in boxes:
        # Inside the box along each axis, allowing the tolerance at both ends.
        inside = [(pos[ax] >= b[2 * ax] - tol) & (pos[ax] <= b[2 * ax + 1] + tol)
                  for ax in range(3)]
        for ax in range(3):
            # On one of this axis's two faces AND within the box in the other two.
            # Testing the face plane alone would flag the whole sample set whenever
            # some map in the case happens to have a boundary at y = 0.
            others = inside[(ax + 1) % 3] & inside[(ax + 2) % 3]
            for edge in (b[2 * ax], b[2 * ax + 1]):
                mask |= (np.abs(pos[ax] - edge) < tol) & others
    return mask


def row_floors(dg: dict, step: float) -> np.ndarray:
    """The smallest matrix difference that means anything, one value per row.

    Derived from G4beamline's PRINTED PRECISION ALONE, from the G4beamline file
    only. It must not be derived from the measured OPALX-minus-G4beamline
    differences: an earlier version of this function did exactly that, taking each
    row's floor to be the worst disagreement in that coordinate, and since the
    matrix difference is built from those same differences the ratio was bounded by
    about one by construction. With a tolerance of 3 the check then could not fail
    -- verified by replacing the OPALX exit state with zeros, with noise, and with
    x doubled, all of which still passed. It measured nothing.

    A row of the matrix is one exit coordinate differenced over the entrance step,
    so its floor is that coordinate's own printing resolution divided by the step.
    The six do not share one: G4beamline prints six significant figures of each
    quantity in its own unit, so a 1.27 m position resolves to 1e-5 m while a
    momentum ratio inherits the 1e-4 MeV/c step on a 28 MeV/c momentum.

    `dg` is the G4beamline exit plane as read_bltrack returns it: m, beta*gamma, s.
    """
    # back to the units G4beamline actually printed, take the step, and return
    dx_m = printing_step(max(abs(dg["x"]).max(), abs(dg["y"]).max()) * 1e3) * 1e-3
    pmax = max(abs(dg[k]).max() for k in ("Px", "Py", "Pz"))          # beta*gamma
    dp_bg = printing_step(pmax * MUON_MASS * 1e3) / (MUON_MASS * 1e3)  # MeV/c -> beta*gamma
    dt_s = printing_step(abs(dg["t"]).max() * 1e9) * 1e-9
    pz = max(abs(dg["Pz"]).min(), 1e-12)
    xp = max(abs(dg["Px"] / dg["Pz"]).max(), abs(dg["Py"] / dg["Pz"]).max())
    f_pos = dx_m / step
    f_ang = dp_bg * (1.0 + xp) / pz / step        # x' = px/pz, both printed
    f_zeta = BETA0 * 299792458.0 * dt_s / step
    f_delta = dp_bg / BG0 / step
    return np.array([f_pos, f_ang, f_pos, f_ang, f_zeta, f_delta])


def significant(M: np.ndarray, floor: float, mult: float = 100.0) -> np.ndarray:
    """Which matrix elements a relative test can say anything about.

    Not `a fraction of the largest element` -- that lets through elements only a
    few times the measurement floor, where the ratio reports the floor rather
    than the codes. An element has to be `mult` times the floor, so the floor
    contributes at most 1/mult of the relative number.
    """
    return np.abs(M) >= mult * floor


def symplectic_residual(M: np.ndarray) -> float:
    """max |M^T J M - J|. A diagnostic, not a test: zeta and delta as defined
    here are conjugate only up to a factor and the differences are taken at a
    finite step, so the residual sits near 1e-3 by construction in both codes.
    Compare the two codes' residuals to each other."""
    J = np.zeros((6, 6))
    for i in range(3):
        J[2 * i, 2 * i + 1] = 1.0
        J[2 * i + 1, 2 * i] = -1.0
    return float(np.abs(M.T @ J @ M - J).max())


def momentum_error(d: dict, code: str) -> np.ndarray:
    """|p|/p0 - 1 per particle. A static magnetic field does no work, so
    comparing this between the entrance and the exit plane tests the integrator
    inside one code, with no reference to the other."""
    if code == "g4bl":
        px, py, pz = d["Px"], d["Py"], d["Pz"]
    else:
        px, py, pz = d["px"], d["py"], d["pz"]
    return np.sqrt(px**2 + py**2 + pz**2) / BG0 - 1.0


# ---------------------------------------------------------------------------
# Field comparison
# ---------------------------------------------------------------------------
def match_field(fg: dict, fo: dict, tol_m: float = 1e-9):
    """Pair up two field dumps by sample position.

    Both dumps are generated from the same specification so the rows ought to
    correspond one to one, but they are matched on position anyway: a silent
    reordering would otherwise read as a field error.
    """
    def key(f):
        return np.round(np.vstack([f["x"], f["y"], f["z"]]).T / tol_m).astype(np.int64)

    kg, ko = key(fg), key(fo)
    index = {tuple(r): i for i, r in enumerate(kg)}
    gi, oi = [], []
    for i, r in enumerate(ko):
        j = index.get(tuple(r))
        if j is not None:
            gi.append(j)
            oi.append(i)
    if not gi:
        raise ValueError("no field sample positions in common between the two dumps")
    gi, oi = np.array(gi), np.array(oi)
    Bg = np.vstack([fg["Bx"][gi], fg["By"][gi], fg["Bz"][gi]])
    Bo = np.vstack([fo["Bx"][oi], fo["By"][oi], fo["Bz"][oi]])
    pos = np.vstack([fo["x"][oi], fo["y"][oi], fo["z"][oi]])
    return pos, Bg, Bo


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def load_manifest(path=None) -> dict:
    path = Path(path or HERE / "cases.json")
    return {c["name"]: c for c in json.loads(path.read_text())}


# ---------------------------------------------------------------------------
# Result table. Lifted from studies/fieldmap_g4bl/fmlib.py so both studies
# report the same way.
# ---------------------------------------------------------------------------
class Results:
    """Rows of (test, quantity, measured, expected, tol) printed with PASS/FAIL.
    A row with tol None is a diagnostic and never fails the suite."""

    def __init__(self):
        self.rows = []

    def check(self, test, name, measured, expected, tol, rel=False, note=""):
        if tol is None:
            ok = None
        elif rel:
            ok = bool(abs(measured - expected) / max(abs(expected), 1e-30) <= tol)
        else:
            ok = bool(abs(measured - expected) <= tol)  # bool(): numpy bools fail `is False`
        self.rows.append(dict(test=test, name=name, measured=measured, expected=expected,
                              tol=tol, rel=rel, ok=ok, note=note))
        return ok

    def at_least(self, test, name, measured, minimum, note=""):
        """A check that something is big enough rather than small enough, for
        the cases where the point is that two results must NOT agree."""
        self.rows.append(dict(test=test, name=name, measured=measured,
                              expected=minimum, tol=None, rel=False,
                              ok=bool(measured >= minimum),
                              note=note or f"must be at least {minimum:g}"))
        return measured >= minimum

    def note(self, test, name, text):
        self.rows.append(dict(test=test, name=name, measured=None, expected=None,
                              tol=None, rel=False, ok=None, note=text))

    @property
    def failed(self) -> int:
        return sum(1 for r in self.rows if r["ok"] is False)

    def print_table(self, title=""):
        if title:
            print(f"\n{'=' * 110}\n{title}\n{'=' * 110}")
        hdr = (f"{'case':16s} {'quantity':34s} {'measured':>13s} {'expected':>13s} "
               f"{'|diff|':>10s} {'tol':>9s}  result")
        print(hdr)
        print("-" * len(hdr))
        for r in self.rows:
            if r["measured"] is None:
                print(f"{r['test']:16s} {r['name']:34s} {r['note']}")
                continue
            diff = abs(r["measured"] - r["expected"])
            res = "diag" if r["ok"] is None else ("PASS" if r["ok"] else "**FAIL**")
            tolstr = "-" if r["tol"] is None else f"{r['tol']:.1e}{'r' if r['rel'] else ''}"
            note = f"  {r['note']}" if r["note"] else ""
            print(f"{r['test']:16s} {r['name']:34s} {r['measured']:>13.6g} "
                  f"{r['expected']:>13.6g} {diff:>10.3g} {tolstr:>9s}  {res}{note}")
