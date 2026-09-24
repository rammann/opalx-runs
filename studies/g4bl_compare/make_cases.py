#!/usr/bin/env python3
"""Write every case in the g4bl_compare study, for both codes, from one description.

One case is one field map, or one group of maps at the relative positions muE4
places them at, isolated on a straight centreline. For each case this writes:

  <case>.g4bl        G4beamline, 19 particles, plus three field sample grids
  <case>.in          OPALX, the same 19 particles, the same three grids
  <case>_gauss.g4bl  G4beamline, 20000 Gaussian muons
  <case>_gauss.in    OPALX, the same 20000 muons
  parts.txt          the 19 particles, OPALX FROMFILE
  beam.txt           the same 19, G4beamline #BLTrackFile
  parts_gauss.txt    the 20000, OPALX
  beam_gauss.txt     the same 20000, G4beamline
  cases.json         one record per case, for run_tests.py and plot_tests.py

The two particle files of a pair are written from one array, so the two codes
cannot be handed different particles.

Both codes are run at a step size where they have converged. That is not a
detail: at G4beamline's default 1 mm and OPALX's 1e-11 s the two disagree by
49 um on the ASR61 dipole, and nearly all of that is each code's own step
error rather than anything to do with the map. See README.md.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from opalxruns.mue4 import (CornerArc, Placement, matmul, parse_rotation,
                            rot_y, to_tait_bryan, walk)
from opalxruns.paths import G4BL_FILES

HERE = Path(__file__).resolve().parent
MAPS = G4BL_FILES / "muE4" / "maps"

MUON_MASS = 0.1056583755
P0 = 0.028
BG0 = P0 / MUON_MASS

# Step sizes at which both codes have converged on the ASR61 dipole, measured
# by halving each in turn until the answer stopped moving. The pair stage uses
# them; the 20000-particle stage uses the coarser pair, where the residual step
# error is 12 um against a 10 mm beam and so cannot affect any distribution number.
DT_FINE, MAXSTEP_FINE = "1e-12", 0.1        # s, mm
DT_COARSE, MAXSTEP_COARSE = "2e-12", 0.5    # s, mm
# A second, finer pair of steps, used only by the cases sensitive enough to need it.
# On the bisector geometry the transverse coordinate stays at 35 mm instead of growing
# to 1.27 m, so G4beamline's six figures resolve 1e-7 m instead of 1e-5 -- and at that
# resolution DT_FINE / MAXSTEP_FINE turn out NOT to be converged after all. Halving
# both halves the disagreement, which is first order in the step and means what is
# left is the integrator, not the map. "Converged" is only ever a statement about the
# resolution you can measure at.
DT_FINER, MAXSTEP_FINER = "5e-13", 0.05    # s, mm

N_GAUSS = 20000
SIG = dict(x=10.0e-3, y=3.0e-3, xp=2.0e-3, yp=2.0e-3)   # m and rad

# ---------------------------------------------------------------------------
# The 19-particle set: the reference, then a small +/- step in each of the six
# coordinates, then a large +/- step in x, y and dp/p.
#
# The small steps give the transfer matrix by centred differences. The large
# ones probe where the map stops being linear, which the small steps cannot see
# and the Gaussian bunch can only see statistically.
#
# The set is exactly symmetric in +/-, so its mean is exactly zero in every
# coordinate. OPALX under FROMFILE takes its reference orbit from the bunch
# mean, and that is what makes the mean coincide with G4beamline's declared
# reference particle. One asymmetric particle puts a silent offset into every
# transverse number.
# ---------------------------------------------------------------------------
STEP_SMALL = 1.0e-3
STEP_LARGE = 2.0e-2


def particle_set() -> np.ndarray:
    """(N, 6) of x, x', y, y', z, dp/p in m, rad and relative."""
    rows = [[0.0] * 6]
    for axis in range(6):                       # x, x', y, y', z, dp/p
        for sign in (+1, -1):
            r = [0.0] * 6
            r[axis] = sign * STEP_SMALL
            rows.append(r)
    for axis in (0, 2, 5):                      # x, y, dp/p
        for sign in (+1, -1):
            r = [0.0] * 6
            r[axis] = sign * STEP_LARGE
            rows.append(r)
    return np.array(rows)


def gauss_set(n: int, seed: int = 20260922) -> np.ndarray:
    """(n, 6). Monoenergetic, and with the sample mean subtracted exactly.

    A raw draw of 20000 is off by sigma/sqrt(n), 0.07 mm in x, which is enough
    to tilt the OPALX reference orbit away from G4beamline's fixed reference
    particle and offset every transverse comparison. Subtracting the sample
    mean, not the population mean, is what removes it.
    """
    rng = np.random.default_rng(seed)
    a = np.zeros((n, 6))
    a[:, 0] = rng.normal(0.0, SIG["x"], n)
    a[:, 1] = rng.normal(0.0, SIG["xp"], n)
    a[:, 2] = rng.normal(0.0, SIG["y"], n)
    a[:, 3] = rng.normal(0.0, SIG["yp"], n)
    a -= a.mean(axis=0)                          # z and dp/p are already all zero
    return a


def write_particles(d: Path, arr: np.ndarray, stem: str) -> None:
    """The same particles as an OPALX FROMFILE list and a G4beamline #BLTrackFile.

    OPALX wants m and beta*gamma, G4beamline mm and MeV/c. x' and y' are turned
    into transverse momenta about the same pz, so |p| carries the dp/p.
    """
    n = len(arr)
    px = np.zeros(n)
    py = np.zeros(n)
    pz = np.zeros(n)
    for i, (x, xp, y, yp, z, dpp) in enumerate(arr):
        p = BG0 * (1.0 + dpp)
        # x' = px/pz, so pz = p / sqrt(1 + x'^2 + y'^2)
        pz[i] = p / math.sqrt(1.0 + xp * xp + yp * yp)
        px[i] = xp * pz[i]
        py[i] = yp * pz[i]
    with open(d / f"{stem}.txt", "w") as f:       # OPALX
        f.write(f"{n}\nx y z px py pz\n")
        for i, (x, xp, y, yp, z, dpp) in enumerate(arr):
            f.write(f"{x:.12e} {y:.12e} {z:.12e} {px[i]:.12e} {py[i]:.12e} {pz[i]:.12e}\n")
    g = "beam" + stem[len("parts"):]
    with open(d / f"{g}.txt", "w") as f:          # G4beamline
        f.write("#BLTrackFile written by make_cases.py\n")
        f.write("#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n")
        f.write("#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n")
        for i, (x, xp, y, yp, z, dpp) in enumerate(arr):
            f.write(f"{x * 1e3:.12g} {y * 1e3:.12g} {z * 1e3:.12g} "
                    f"{px[i] * MUON_MASS * 1e3:.12g} {py[i] * MUON_MASS * 1e3:.12g} "
                    f"{pz[i] * MUON_MASS * 1e3:.12g} 0 -13 {i + 1} 1 0 1\n")


# ---------------------------------------------------------------------------
# Map headers
# ---------------------------------------------------------------------------
def map_header(name: str) -> dict:
    """The grid line of a G4beamline `grid` map: where it starts, how many
    points and how far apart, per axis, in mm."""
    with open(MAPS / name) as f:
        for line in f:
            t = line.split()
            if t and t[0] == "grid":
                kv = dict(p.split("=") for p in t[1:])
                return {k: float(v) if k[0] in "XYZd" and not k.startswith("n") else int(float(v))
                        for k, v in kv.items()}
            if t and t[0] == "cylinder":
                raise ValueError(f"{name}: cylinder maps are covered by wsx_solenoid/")
    raise ValueError(f"{name}: no grid header")


SIGNS = {"": (1, 1, 1), "Y180": (-1, 1, -1), "Z180": (-1, -1, 1), "Y180,Z180": (1, -1, -1)}


def lab_grid(h: dict, rot: str, origin_mm: tuple[float, float, float]):
    """Where the map's own grid points land in the lab, as (start, step, count)
    per axis, in mm.

    The only rotations muE4 uses on these maps are 180 degrees about y or z or
    both, which flip signs without swapping axes, so a map grid point is always
    a lab grid point and `at the map's own grid points` stays meaningful.
    """
    out = []
    for ax, key, nkey, dkey in (("x", "X0", "nX", "dX"), ("y", "Y0", "nY", "dY"),
                                ("z", "Z0", "nZ", "dZ")):
        s, n, d = h[key], h[nkey], h[dkey]
        sign = SIGNS[rot]["xyz".index(ax)]
        o = origin_mm["xyz".index(ax)]
        if sign > 0:
            out.append((s + o, d, n))
        else:
            out.append((-(s + (n - 1) * d) + o, d, n))
    return out


# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------
ASR61_D = "asr61_300d_track.g4blmap"
ASR61_SM = "asr61_300sm_track.g4blmap"
ASR62 = "asr62shim_280_track.g4blmap"
ASR62_SM = "asr62shim_280_sm_track.g4blmap"
QSM = "qsm01a_210_track.g4blmap"

CASES = [
    dict(name="asr61_dipole",
         desc="muE4 ASR61_300d dipole map on its own",
         maps=[dict(file=ASR61_D, rot="", scale=-9.9528e-01, xyz=(0.0, 0.0, 1500.0))],
         monitors=[50, 1000, 1500, 2000, 2500, 2950], zstop=4.0,
         world=dict(h=1000, w=4000, l=4200, x=1000, z=1500), radius_cut=2000),
    dict(name="asr61_300sm",
         desc="muE4 ASR61_300sm map, moved onto the axis so a beam can reach it",
         maps=[dict(file=ASR61_SM, rot="", scale=-9.9528e-01, xyz=(-505.0, 0.0, 1880.0))],
         monitors=[400, 700, 900, 1100, 1330], zstop=1.6,
         world=dict(h=600, w=1200, l=2000, x=0, z=800), radius_cut=500),
    dict(name="asr62_dipole",
         desc="muE4 ASR62 dipole map on its own, at the polarity of the third bend",
         maps=[dict(file=ASR62, rot="", scale=9.6187e-01, xyz=(0.0, 0.0, 1550.0))],
         monitors=[50, 1000, 1550, 2100, 2600, 3050], zstop=4.0,
         world=dict(h=1000, w=4000, l=4400, x=-1000, z=1550), radius_cut=2000),
    dict(name="qsm600_quad",
         desc="muE4 QSM600 quadrupole map, the one file with nine columns",
         maps=[dict(file=QSM, rot="", scale=-9.053e-01, xyz=(0.0, 0.0, 850.0))],
         monitors=[50, 400, 850, 1300, 1650], zstop=2.0,
         world=dict(h=1200, w=1200, l=2200, x=0, z=900), radius_cut=500),
    dict(name="asr61_group",
         desc="ASR61_300d and both ASR61_300sm placements, as muE4 places them",
         maps=[dict(file=ASR61_D, rot="", scale=-9.9528e-01, xyz=(0.0, 0.0, 1500.0)),
               dict(file=ASR61_SM, rot="", scale=-9.9528e-01, xyz=(0.0, 0.0, 1500.0)),
               dict(file=ASR61_SM, rot="Y180,Z180", scale=9.9528e-01, xyz=(0.0, 0.0, 1500.0))],
         monitors=[50, 1000, 1500, 2000, 2500, 2950], zstop=4.0,
         world=dict(h=1000, w=4000, l=4200, x=1000, z=1500), radius_cut=2000),
    dict(name="asr62_d2_group",
         desc="ASR62 and both ASR62_sm placements of the second bend, all turned by Y180",
         maps=[dict(file=ASR62, rot="Y180", scale=-9.6380e-01, xyz=(0.0, 0.0, 1550.0)),
               dict(file=ASR62_SM, rot="Y180", scale=-9.6380e-01, xyz=(0.0, 0.0, 1550.0)),
               dict(file=ASR62_SM, rot="Z180", scale=9.6380e-01, xyz=(0.0, 0.0, 1550.0))],
         monitors=[50, 1000, 1550, 2100, 2600, 3050], zstop=4.0,
         world=dict(h=1000, w=4000, l=4400, x=-1000, z=1550), radius_cut=2000),
    dict(name="asr62_d3_group",
         desc="ASR62 and both ASR62_sm placements of the third bend, none turned",
         maps=[dict(file=ASR62, rot="", scale=9.6187e-01, xyz=(0.0, 0.0, 1550.0)),
               dict(file=ASR62_SM, rot="", scale=9.6187e-01, xyz=(0.0, 0.0, 1550.0)),
               dict(file=ASR62_SM, rot="Y180,Z180", scale=9.6187e-01, xyz=(0.0, 0.0, 1550.0))],
         monitors=[50, 1000, 1550, 2100, 2600, 3050], zstop=4.0,
         world=dict(h=1000, w=4000, l=4400, x=-1000, z=1550), radius_cut=2000),
]


def bisector_layout(case: dict):
    """Resolve a case written the way muE4 writes it -- cornerarc commands and
    placements at centreline z -- into lab positions and rotations.

    This reuses walk() and to_tait_bryan() from the whole-line study rather than
    deriving the frame here, so the case tests the same conversion a real muE4 run
    depends on. A monitor is just a placement with no rotation of its own, which
    means walk() gives it the centreline frame at its z -- and an OPALX MONITOR
    records in its own local frame, so the two codes report in the same frame with
    nothing to transform.

    Returns the Placement objects in deck order, each carrying .lab in mm and
    .pose as (THETA, PHI, PSI).
    """
    arcs, places, line = [], [], 0
    for item in case["sequence"]:
        line += 1
        if item["kind"] == "arc":
            arcs.append(CornerArc(cl_z=item["cl_z"], angle_deg=item["angle_deg"],
                                  radius=item["radius"], line=line))
        else:
            places.append(Placement(
                index=len(places), solid=None, name=item["name"], cl_z=item["cl_z"],
                off_x=item.get("off_x", 0.0), off_y=item.get("off_y", 0.0),
                rotation=item.get("rot", ""), current=item.get("scale"), line=line))
    walk(places, arcs)
    return places


def bisector_boxes(case: dict, places):
    """Each placed map's lab bounding box, from its eight corners put through the
    placement. Axis aligned only because run_tests uses it to spot sample points on
    a face; at 20 degrees the box is not aligned, so this is a bound, not the box."""
    out = []
    by_name = {p.name: p for p in places}
    for item in case["sequence"]:
        if item["kind"] != "map":
            continue
        h = map_header(item["file"])
        p = by_name[item["name"]]
        th = math.radians(p.frame_deg)
        r = matmul(rot_y(th), parse_rotation(item.get("rot", "")))
        xs, ys, zs = [], [], []
        for cx in (h["X0"], h["X0"] + (h["nX"] - 1) * h["dX"]):
            for cy in (h["Y0"], h["Y0"] + (h["nY"] - 1) * h["dY"]):
                for cz in (h["Z0"], h["Z0"] + (h["nZ"] - 1) * h["dZ"]):
                    v = [r[i][0] * cx + r[i][1] * cy + r[i][2] * cz for i in range(3)]
                    xs.append(p.lab[0] + v[0])
                    ys.append(p.lab[1] + v[1])
                    zs.append(p.lab[2] + v[2])
        out.append([min(xs) / 1000, max(xs) / 1000, min(ys) / 1000, max(ys) / 1000,
                    min(zs) / 1000, max(zs) / 1000])
    return out


def bisector_grid(case: dict, places, step_mm: float = 20.0):
    """One lab-frame plane of sample points covering the magnet.

    Deliberately NOT the map's own grid: at 20 degrees the map's points do not land
    on a lab-aligned grid, so nothing here reads the table back -- the straight-frame
    cases do that. What this tests is that both codes put the SAME field at the same
    lab point, which is a direct check of the cornerarc-to-pose conversion with no
    tracking involved. A wrong THETA shows up here immediately.
    """
    b = bisector_boxes(case, places)
    x0 = min(v[0] for v in b) * 1000
    x1 = max(v[1] for v in b) * 1000
    z0 = min(v[4] for v in b) * 1000
    z1 = max(v[5] for v in b) * 1000
    # one step inside every face, so no sample lands on a box edge
    nx = max(2, int((x1 - x0) / step_mm) - 1)
    nz = max(2, int((z1 - z0) / step_mm) - 1)
    return dict(frame=dict(x=(x0 + step_mm, step_mm, nx), y=(0.0, step_mm, 1),
                           z=(z0 + step_mm, step_mm, nz)))


def lab_boxes(case: dict):
    """Every placed map's bounding box in lab metres, for run_tests to spot
    sample points that land on a face."""
    out = []
    for m in case["maps"]:
        h = map_header(m["file"])
        g = lab_grid(h, m["rot"], m["xyz"])
        b = []
        for (s0, d, n) in g:
            b += [s0 / 1000.0, (s0 + (n - 1) * d) / 1000.0]
        out.append(b)
    return out


# muE4's own numbers, from mue4_WsxOn.g4bl lines 143-146. The magnet sits at the
# endpoint of the first arc, which is why it lands on the bisector of the 40 degree
# total bend. Recording planes at centreline z = 1500 and 4600: both clear of the
# field, and both outside the arcs, where a chord frame would not be defined.
ARC_IN = dict(kind="arc", cl_z=2578.0774, angle_deg=20.0, radius=1159.0622)
ARC_OUT = dict(kind="arc", cl_z=2982.6664, angle_deg=20.0, radius=1159.0621)
MON_IN = dict(kind="mon", name="MON_1500", cl_z=1500.0)
MON_OUT = dict(kind="mon", name="MON_4600", cl_z=4600.0)
CL_Z = 2982.6664

BISECTOR_CASES = [
    dict(name="asr61_bisector",
         desc="ASR61_300d on the bisector, in muE4's own cornerarc frame",
         bisector=True, zstop=4.900,
         sequence=[MON_IN, ARC_IN,
                   dict(kind="map", name="ASR61_300d", file=ASR61_D, cl_z=CL_Z,
                        scale=-9.9528e-01, rot=""),
                   ARC_OUT, MON_OUT],
         world=dict(h=1000, w=5000, l=7000, x=1000, z=2500), radius_cut=2500),
    dict(name="asr61_group_bisector",
         desc="ASR61_300d and both ASR61_300sm placements, exactly as muE4 places them",
         bisector=True, zstop=4.900,
         sequence=[MON_IN, ARC_IN,
                   dict(kind="map", name="ASR61_300d", file=ASR61_D, cl_z=CL_Z,
                        scale=-9.9528e-01, rot=""),
                   dict(kind="map", name="ASR61_300sm", file=ASR61_SM, cl_z=CL_Z,
                        scale=-9.9528e-01, rot=""),
                   dict(kind="map", name="ASR61_300sm_m", file=ASR61_SM, cl_z=CL_Z,
                        scale=9.9528e-01, rot="Y180,Z180"),
                   ARC_OUT, MON_OUT],
         world=dict(h=1000, w=5000, l=7000, x=1000, z=2500), radius_cut=2500),
]


def sample_grids(case: dict) -> dict:
    """Three sets of points at which both codes are asked for the field.

    `grid` sits on the map's own grid points, so a reader that loads the table
    correctly returns the tabulated value and nothing is being interpolated.
    `half` sits halfway between them in x and z, where both codes have to
    interpolate, so a difference there with `grid` passing means the two
    interpolation schemes differ rather than the tables.
    `yscan` walks the full height at a few places, which is the only set that
    tests the y direction.

    Sampling stops one point short of the top of every axis: both OPALX readers
    exclude the top face, because interpolation needs a whole cell above the point.
    """
    m = case["maps"][0]
    h = map_header(m["file"])
    (x0, dx, nx), (y0, dy, ny), (z0, dz, nz) = lab_grid(h, m["rot"], m["xyz"])
    g = dict(
        grid=dict(x=(x0, dx, nx - 1), y=(y0 + (ny // 2) * dy, dy, 1), z=(z0, dz, nz - 1)),
        # Offset in ALL THREE axes, y included. With y left on a grid plane the
        # trilinear weight in y is always 0 or 1, so the y direction of the
        # interpolation is never exercised -- and on the midplane of a dipole Bx
        # and Bz are zero by symmetry, so the dense scan would be testing one
        # component of three. Half a step off it, all three are nonzero and all
        # three interpolation weights are in play.
        half=dict(x=(x0 + dx / 2, dx, nx - 1), y=(y0 + (ny // 2) * dy + dy / 2, dy, 1),
                  z=(z0 + dz / 2, dz, nz - 1)),
        yscan=dict(x=three_inside(x0, dx, nx), y=(y0 + dy, dy, max(1, ny - 2)),
                   z=three_inside(z0, dz, nz)),
        # Deliberately ON the top x face. For a single unrotated map both codes
        # return zero there, although the file does hold values -- the ASR61 map has
        # 0.033 T at x = 390. Both exclude the top face: the OPALX readers by
        # `r < end`, never `<=`, because interpolating needs a whole cell above the
        # point, and G4beamline likewise. Where a map is placed with a 180 degree
        # rotation the two codes do NOT agree here, because sin(pi) is 1.22e-16
        # rather than 0 and the sample lands a few 1e-16 m either side of the face.
        # Kept as a diagnostic that reports what each code returned.
        #
        # What the two codes do NOT agree on is a point a few 1e-14 m inside a face:
        # there OPALX interpolates and G4beamline returns zero. That is a knife edge
        # 0.04 picometres wide and no particle can land on it, but a sample grid can,
        # which is why three_inside() keeps the scans a whole step clear of both ends.
        # Sampling x at -190 + 2*290 lands on 389.99999999999996 and reads as a 0.05 T
        # disagreement that is entirely the last bit of a float.
        edge=dict(x=(x0 + (nx - 1) * dx, dx, 1), y=(y0, dy, ny - 1), z=(z0, dz, nz - 1)),
    )
    return g


def three_inside(start, step, n):
    """Three sample positions spread across an axis, never on either end face.

    Index 0 is fine but index n-1 is the top face, which the OPALX readers
    exclude, so a scan that lands there reads as a field error when it is only
    the documented boundary rule.
    """
    if n < 4:
        return (start + step, step, max(1, n - 2))
    return (start + step, step * ((n - 3) // 2), 3)


def g4bl_range(spec):
    s, d, n = spec
    if n == 1:
        return f"{s:.6f}"
    # a hundredth of a step of margin, so the last point is not lost to rounding
    return f"{s:.6f},{s + (n - 1) * d + d / 100.0:.6f},{d:.6f}"


def dump_block(name, spec, fname):
    x, y, z = spec["x"], spec["y"], spec["z"]
    return (f'DUMPEMFIELDS, FILE_NAME = "{fname}", COORDINATE_SYSTEM = CARTESIAN,\n'
            f'    X_START = {x[0] / 1000:.9f}, DX = {x[1] / 1000:.9f}, X_STEPS = {x[2]},\n'
            f'    Y_START = {y[0] / 1000:.9f}, DY = {y[1] / 1000:.9f}, Y_STEPS = {y[2]},\n'
            f'    Z_START = {z[0] / 1000:.9f}, DZ = {z[1] / 1000:.9f}, Z_STEPS = {z[2]},\n'
            f'    T_START = 0.0, DT = 1.0, T_STEPS = 1;\n')


def elements(case: dict):
    """Monitors and maps in one list, sorted by lab z and numbered.

    OPALX sorts a line placed by absolute position alphabetically by element
    name rather than by position, so the number in front is what keeps every
    dump in beam order.
    """
    items = [("MON", z, None) for z in case["monitors"]]
    items += [("MAP", m["xyz"][2], m) for m in case["maps"]]
    items.sort(key=lambda t: t[1])
    return [(f"E{i + 1:02d}", kind, z, m) for i, (kind, z, m) in enumerate(items)]


def write_opalx_bisector(case, d: Path, gauss: bool, finer: bool = False):
    """The OPALX deck for a case written the way muE4 writes it.

    Every position and rotation comes from walk(), the same resolver the whole-line
    study uses, so if the cornerarc-to-pose conversion were wrong this deck would be
    wrong in the same way -- which is the point: the straight-frame cases cannot see
    that conversion at all, because every rotation in them is 0 or 180 degrees.
    """
    stem = case["name"] + ("_gauss" if gauss else "_fine" if finer else "")
    # NOT DT_COARSE for the 20000-particle stage here. The coarse step was chosen
    # because what it leaves, 12 um, cannot reach any number measured against a
    # 1e-5 m floor. On this geometry the floor is 1e-7 m, and the coarse step leaves
    # 2e-4 m -- three orders above it. A step is only "good enough" relative to the
    # resolution of the thing it is being judged by.
    dt = DT_FINER if finer else DT_FINE
    npart = N_GAUSS if gauss else 19
    parts = "parts_gauss.txt" if gauss else "parts.txt"
    places = bisector_layout(case)
    by_name = {p.name: p for p in places}

    lines = [f"/*  {stem}.in -- {case['desc']}.", "",
             f"    Paired with {stem}.g4bl: the same {npart} muons, the same map files, and",
             "    the same geometry muE4 places them in -- the two cornerarc commands that",
             "    turn the centreline 20 degrees before the magnet and 20 degrees after it,",
             "    so the magnet sits on the bisector of the total 40 degree bend and the beam",
             "    enters its map at -20 degrees rather than along the map axis.", "",
             "    G4beamline has cornerarc; OPALX does not, so the frame has to be expressed",
             "    as an absolute position and rotation per element. Those come from",
             "    mue4lib.walk() and to_tait_bryan(), and they reproduce the muE4 lattice's",
             "    own numbers exactly: X = 0.069900004, Z = 2.974500020, THETA = 0.349065850.",
             "    That conversion is what this case tests and the straight-frame cases cannot.",
             "",
             "    Both codes report in the CENTRELINE frame: G4beamline through",
             "    `coordinates=centerline`, OPALX because a MONITOR records in its own local",
             "    frame and walk() gives each monitor the centreline frame at its z. The",
             "    transverse numbers therefore stay in millimetres instead of growing to the",
             "    1.27 m a straight frame would give, and G4beamline's six printed figures",
             "    resolve them about a hundred times more finely.", "",
             "    GENERATED by make_cases.py -- do not edit.", "*/", ""]
    lines += ["OPTION, PSDUMPFREQ    = 100000;",
              "OPTION, STATDUMPFREQ  = 1;",
              "OPTION, BOUNDPDESTROY = 1000000;",
              "OPTION, VERSION       = 10900;", "",
              f'Title, string="g4bl_compare / {stem}";', "",
              "REAL MUON_MASS = 0.1056583755;",
              "REAL P0        = 0.028;", ""]

    names = []
    for i, item in enumerate(case["sequence"]):
        if item["kind"] == "arc":
            lines.append(f"// g4bl: cornerarc z={item['cl_z']} angle={item['angle_deg']} "
                         f"centerRadius={item['radius']}  (no OPALX element; the frame it "
                         f"creates is folded into the poses below)")
            continue
        pl = by_name[item["name"]]
        nm = f"E{len(names) + 1:02d}_{item['name'].upper()}"
        x, y, z = (v / 1000.0 for v in pl.lab)
        th, ph, ps = pl.pose
        pose = (f"X = {x:.9f}, Y = {y:.9f}, Z = {z:.9f},\n"
                f"    THETA = {th:.9f}, PHI = {ph:.9f}, PSI = {ps:.9f},")
        if item["kind"] == "mon":
            lines += [f"// centreline z = {item['cl_z']:.4f} mm, frame {pl.frame_deg:.2f} deg",
                      f"{nm}: MONITOR, {pose}",
                      f'    DELETEONTRANSVERSEEXIT = FALSE, OUTFN = "MON_{item["cl_z"]:.0f}";', ""]
        else:
            lines += [f"// g4bl: place {item['name']} z={item['cl_z']}"
                      + (f" rotation={item['rot']}" if item.get("rot") else "")
                      + f" current={item['scale']:.5E}   (frame {pl.frame_deg:.2f} deg)",
                      f"{nm}: FIELDMAP, {pose}",
                      f'    FMAPFN = "{MAPS / item["file"]}",',
                      f"    SCALE = {item['scale']:.6E};", ""]
        names.append(nm)

    lines += [f"CaseLine : LINE = ({', '.join(names)});", ""]
    if not gauss and not finer:
        g = bisector_grid(case, places)
        lines += ["// A lab-frame plane through the magnet. NOT the map's own grid points --",
                  "// at 20 degrees they do not land on a lab-aligned grid, and reading the",
                  "// table back is what the straight-frame cases do. What this tests is that",
                  "// both codes put the same field at the same lab point, which is the",
                  "// cornerarc-to-pose conversion with no tracking in the way.",
                  dump_block("frame", g["frame"], "opalx_field_frame.dat"), ""]

    lines += ["FS1: FIELDSOLVER, TYPE = NONE, NX = 16, NY = 16, NZ = 16,",
              "     PARFFTX = true, PARFFTY = true, PARFFTZ = true,",
              "     BCFFTX = open, BCFFTY = open, BCFFTZ = open,",
              "     BBOXINCR = 1, GREENSF = INTEGRATED;", "",
              f"REAL n_particles = {npart};",
              f'Dist1: DISTRIBUTION, TYPE = FROMFILE, FNAME = "{parts}", NPARTDIST = n_particles;',
              "ES1: EMISSIONSOURCE, DISTRIBUTION = Dist1;",
              "mySources: EMISSIONSOURCELIST = (ES1);", "",
              "BEAM1: BEAM, PARTICLE = MUON, NALLOC = n_particles,",
              "       BCHARGE = 1e-15, SOURCES = mySources, CHARGE = 1;", "",
              "TRACK, LINE = CaseLine, BEAM = BEAM1,",
              "       MAXSTEPS = {8000000},",
              f"       DT = {{{dt}}},",
              f"       ZSTOP = {{{case['zstop']:.3f}}};",
              'RUN, METHOD = "PARALLEL", FIELDSOLVER = FS1;',
              "ENDTRACK;", "Quit;", ""]
    (d / f"{stem}.in").write_text("\n".join(lines))


def write_g4bl_bisector(case, d: Path, gauss: bool, finer: bool = False):
    """The G4beamline deck: muE4's own cornerarc construction, verbatim."""
    stem = case["name"] + ("_gauss" if gauss else "_fine" if finer else "")
    ms = MAXSTEP_FINER if finer else MAXSTEP_FINE
    beam = "beam_gauss.txt" if gauss else "beam.txt"
    w = case["world"]
    places = bisector_layout(case)
    mons = [it["cl_z"] for it in case["sequence"] if it["kind"] == "mon"]

    lines = [f"# {stem}.g4bl -- {case['desc']}.",
             "#",
             "# The two cornerarc commands are muE4's own, copied with their numbers. They",
             "# turn the CENTRELINE, not the beam: the beam carries straight on until the",
             "# field bends it, so it enters the map at -20 degrees to the map's own axis and",
             "# leaves at +20, which is the trajectory muE4 actually uses. A straight-frame",
             "# case sends the beam in along the map axis instead, and it leaves the box",
             "# sideways.",
             "#",
             "# coordinates=centerline keeps the transverse numbers in millimetres. In the",
             "# global frame six significant figures on a metre-sized coordinate would be a",
             "# 1e-5 m floor, larger than the difference being measured.",
             "#",
             "# GENERATED by make_cases.py -- do not edit.", "",
             "physics QGSP_BERT disable=Decay doStochastics=0",
             "param worldMaterial=Vacuum",
             f"param maxStep={ms} deltaChord=0.001", "",
             f"start x=0 y=0 z=0 initialZ=0 radiusCut={case['radius_cut']}",
             "reference particle=mu+ referenceMomentum=28",
             "trackcuts keep=mu+ killSecondaries=1", "",
             f"beam ascii filename={beam}", "",
             f"box WORLDBOX height={w['h']} width={w['w']} length={w['l']} material=Vacuum",
             f"place WORLDBOX z={w['z']} x={w['x']}", ""]
    seen = {}
    for item in case["sequence"]:
        if item["kind"] == "map" and item["file"] not in seen:
            seen[item["file"]] = f"M{len(seen)}"
            lines.append(f"fieldmap {seen[item['file']]} file={MAPS / item['file']}")
    lines.append("")
    for item in case["sequence"]:
        if item["kind"] == "arc":
            lines.append(f"cornerarc z={item['cl_z']} angle={item['angle_deg']} "
                         f"centerRadius={item['radius']}")
        elif item["kind"] == "map":
            rot = f" rotation={item['rot']}" if item.get("rot") else ""
            lines.append(f"place {seen[item['file']]} z={item['cl_z']}{rot} "
                         f"current={item['scale']:.5E}")
    lines += ["",
              "# Recording planes, in centreline coordinates, both clear of the field and both",
              "# outside the arcs -- a plane inside an arc has no well defined chord frame.",
              f"zntuple z={','.join(f'{z:.0f}' for z in mons)} "
              f"format=ascii coordinates=centerline", ""]
    if not gauss and not finer:
        g = bisector_grid(case, places)["frame"]
        lines += ["# The same lab-frame plane the OPALX deck asks for.",
                  f"fieldntuple FFRAME format=ascii filename=g4bl_field_frame.txt "
                  f"x={g4bl_range(g['x'])} y={g4bl_range(g['y'])} z={g4bl_range(g['z'])} t=0"]
    lines.append("")
    (d / f"{stem}.g4bl").write_text("\n".join(lines))


def write_opalx(case, d: Path, gauss: bool, finer: bool = False):
    stem = case["name"] + ("_gauss" if gauss else "")
    dt = DT_COARSE if gauss else DT_FINE
    npart = N_GAUSS if gauss else 19
    parts = "parts_gauss.txt" if gauss else "parts.txt"
    grids = sample_grids(case)

    lines = [f"/*  {stem}.in -- {case['desc']}.", ""]
    lines += [f"    Paired with {stem}.g4bl: the same {npart} muons and the same map files go",
              "    through both codes. Nothing else can differ -- no space charge, no material,",
              "    no decay, no scraping.", "",
              "    A FIELDMAP element takes neither ELEMEDGE nor L, and OPALX allows one",
              "    placement convention per beamline, so every element here is placed by",
              "    absolute position and rotation. The element's local frame is the map's own",
              "    frame, so Z is the lab position of the map's z = 0.", "",
              *( [f"    DT = {dt} s: the step at which OPALX has stopped moving. At 1e-11 it is",
                  "    still 30 um short of its own limit on this geometry, which would read as a",
                  "    disagreement with G4beamline that has nothing to do with the map.", ""]
                 if not gauss else
                 [f"    DT = {dt} s, one step coarser than the 19-particle deck. What is left",
                  "    of the step error there is 12 um, against a bunch 10 mm across, so it",
                  "    cannot reach any number this stage measures -- and the stage is 1000",
                  "    times more particles.", ""] ),
              "    GENERATED by make_cases.py -- do not edit.", "*/", ""]
    lines += ["OPTION, PSDUMPFREQ    = 100000;",
              "OPTION, STATDUMPFREQ  = 1;        // reference orbit and its field, every step",
              "OPTION, BOUNDPDESTROY = 1000000;  // never trip the statistical loss box",
              "OPTION, VERSION       = 10900;", "",
              f'Title, string="g4bl_compare / {stem}";', "",
              "REAL MUON_MASS = 0.1056583755;              // GeV",
              "REAL P0        = 0.028;                     // GeV/c; FROMFILE forbids PC on BEAM",
              "REAL BG0       = P0 / MUON_MASS;",
              "REAL BRHO      = P0 / 0.299792458;", "",
              "value, {P0, BG0, BRHO};", ""]

    names = []
    for tag, kind, z, m in elements(case):
        if kind == "MON":
            nm = f"{tag}_MON_{z:.0f}"
            lines += [f"{nm}: MONITOR, X = 0.0, Y = 0.0, Z = {z / 1000:.6f},",
                      "    THETA = 0.0, PHI = 0.0, PSI = 0.0,",
                      f'    DELETEONTRANSVERSEEXIT = FALSE, OUTFN = "MON_{z:.0f}";', ""]
        else:
            th, ph, ps = to_tait_bryan(parse_rotation(m["rot"])) if m["rot"] else (0.0, 0.0, 0.0)
            nm = f"{tag}_{Path(m['file']).stem[:12].upper()}"
            x, y, zz = m["xyz"]
            lines += [f"// g4bl: place ... z={zz:.4f}"
                      + (f" rotation={m['rot']}" if m["rot"] else "")
                      + f" current={m['scale']:.5E}",
                      f"{nm}: FIELDMAP, X = {x / 1000:.9f}, Y = {y / 1000:.9f}, "
                      f"Z = {zz / 1000:.9f},",
                      f"    THETA = {th:.9f}, PHI = {ph:.9f}, PSI = {ps:.9f},",
                      f'    FMAPFN = "{MAPS / m["file"]}",',
                      f"    SCALE = {m['scale']:.6E};", ""]
        names.append(nm)

    lines += [f"CaseLine : LINE = ({', '.join(names)});", ""]

    if not gauss:
        lines += ["// Field sampling, in lab metres. OPALX writes these into data/, not here.",
                  "// `grid` lands on the map's own grid points, `half` halfway between them,",
                  "// `yscan` walks the full height. See make_cases.py for why."]
        for k in ("grid", "half", "yscan", "edge"):
            lines.append(dump_block(k, grids[k], f"opalx_field_{k}.dat"))
        lines.append("")

    lines += ["FS1: FIELDSOLVER, TYPE = NONE, NX = 16, NY = 16, NZ = 16,",
              "     PARFFTX = true, PARFFTY = true, PARFFTZ = true,   // required even for NONE",
              "     BCFFTX = open, BCFFTY = open, BCFFTZ = open,",
              "     BBOXINCR = 1, GREENSF = INTEGRATED;", "",
              f"REAL n_particles = {npart};",
              f'Dist1: DISTRIBUTION, TYPE = FROMFILE, FNAME = "{parts}", NPARTDIST = n_particles;',
              "ES1: EMISSIONSOURCE, DISTRIBUTION = Dist1;   // EMISSIONMODEL stays NONE for FROMFILE",
              "mySources: EMISSIONSOURCELIST = (ES1);", "",
              "// CHARGE = 1 is mandatory: ParticleProperties maps MUON to -1 and Beam::execute",
              "// only consults that table when the attribute is absent. Check run.log for",
              '// "CHARGE      +e * 1" -- a charge sign error looks exactly like a map that is',
              "// turned the wrong way round.",
              "BEAM1: BEAM, PARTICLE = MUON, NALLOC = n_particles,",
              "       BCHARGE = 1e-15, SOURCES = mySources, CHARGE = 1;", "",
              "// ZSTOP is accumulated path length, not lab z. The orbit turns, so s runs ahead",
              "// of z and the two stop being the same number once the field starts to bend it.",
              "TRACK, LINE = CaseLine, BEAM = BEAM1,",
              "       MAXSTEPS = {8000000},",
              f"       DT = {{{dt}}},",
              f"       ZSTOP = {{{case['zstop']:.3f}}};",
              'RUN, METHOD = "PARALLEL", FIELDSOLVER = FS1;',
              "ENDTRACK;", "Quit;", ""]
    (d / f"{stem}.in").write_text("\n".join(lines))


def write_g4bl(case, d: Path, gauss: bool, finer: bool = False):
    # A straight-frame case has no fine stage: it cannot resolve below 1e-5 m, so a
    # smaller step would show nothing. finer is accepted only to match the signature.
    stem = case["name"] + ("_gauss" if gauss else "")
    ms = MAXSTEP_COARSE if gauss else MAXSTEP_FINE
    beam = "beam_gauss.txt" if gauss else "beam.txt"
    w = case["world"]
    grids = sample_grids(case)

    lines = [f"# {stem}.g4bl -- {case['desc']}.",
             "#",
             f"# Paired with {stem}.in: the same muons and the same map files go through both",
             "# codes. Vacuum world, no material, no decay, no stochastic processes, so only",
             "# the field can cause a difference.",
             "#",
             "# FRAME: z here [mm] / 1000 == path length s in the OPALX input [m], until the",
             "#   field starts to bend the orbit and s runs ahead of z. The centreline is",
             "#   straight -- no cornerarc -- which is what makes the two frames one to one.",
             "#",
             f"# maxStep={ms} mm: the step at which G4beamline has stopped moving. Its default",
             "#   is 100 mm and this deck's ancestors used 1 mm; at 1 mm the answer is still",
             "#   30 um from its own limit, which is larger than anything being measured.",
             "#   NOTE the capital S -- `maxstep` is a different, unused parameter.",
             "#",
             "# GENERATED by make_cases.py -- do not edit.", "",
             "physics QGSP_BERT disable=Decay doStochastics=0",
             "param worldMaterial=Vacuum",
             f"param maxStep={ms} deltaChord=0.001", "",
             "# `start` must come before every `place` or g4bl aborts with `Invalid start`.",
             f"start x=0 y=0 z=0 initialZ=0 radiusCut={case['radius_cut']}",
             "reference particle=mu+ referenceMomentum=28",
             "trackcuts keep=mu+ killSecondaries=1", "",
             "# renumber=0 (the default) keeps EventID from the file, so EventID == OPALX id + 1.",
             "# Do NOT set beamZ: it forces z=0 on every particle and wipes out any longitudinal",
             "# offset the file carries.",
             f"beam ascii filename={beam}", "",
             "# `fieldmap` places no physical volume, so without this nothing expands the world",
             "# and every track leaves immediately, leaving the ntuples silently empty.",
             f"box WORLDBOX height={w['h']} width={w['w']} length={w['l']} material=Vacuum",
             f"place WORLDBOX z={w['z']} x={w['x']}", ""]

    # One `fieldmap` definition per distinct file, however many times it is placed.
    seen = {}
    for m in case["maps"]:
        if m["file"] not in seen:
            seen[m["file"]] = f"M{len(seen)}"
            lines.append(f"fieldmap {seen[m['file']]} file={MAPS / m['file']}")
    lines.append("")
    lines.append("# current= on the place line is the deck-side multiplier on the tabulated")
    lines.append("# Tesla, the same quantity OPALX spells SCALE. The readers do not normalise.")
    for m in case["maps"]:
        x, y, z = m["xyz"]
        rot = f" rotation={m['rot']}" if m["rot"] else ""
        xs = f" x={x}" if x else ""
        ys = f" y={y}" if y else ""
        lines.append(f"place {seen[m['file']]} z={z}{xs}{ys}{rot} current={m['scale']:.5E}")
    lines.append("")
    lines += ["# Recording planes. zntuple interpolates each track linearly onto the exact",
              "# plane, the same construction an OPALX MONITOR uses. One file per plane,",
              "# named after it: Z50.txt and so on.",
              f"zntuple z={','.join(str(z) for z in case['monitors'])} "
              f"format=ascii coordinates=centerline", ""]
    if not gauss:
        lines += ["# The same three sample sets the OPALX deck asks for, in the same order.",
                  "# No tracking is involved, so a difference here is a field difference and",
                  "# everything below it is meaningless until it passes."]
        for k in ("grid", "half", "yscan", "edge"):
            s = grids[k]
            lines.append(f"fieldntuple F{k.upper()} format=ascii filename=g4bl_field_{k}.txt "
                         f"x={g4bl_range(s['x'])} y={g4bl_range(s['y'])} "
                         f"z={g4bl_range(s['z'])} t=0")
    lines.append("")
    (d / f"{stem}.g4bl").write_text("\n".join(lines))


def main():
    manifest = []
    pset = particle_set()
    gset = gauss_set(N_GAUSS)
    for case in CASES + BISECTOR_CASES:
        d = HERE / case["name"]
        d.mkdir(exist_ok=True)
        write_particles(d, pset, "parts")
        write_particles(d, gset, "parts_gauss")
        bis = case.get("bisector", False)
        for g in (False, True):
            (write_opalx_bisector if bis else write_opalx)(case, d, gauss=g)
            (write_g4bl_bisector if bis else write_g4bl)(case, d, gauss=g)
        if bis:
            # A second pair-stage deck at half the step, so the harness can show
            # what is step error rather than assert it.
            write_opalx_bisector(case, d, gauss=False, finer=True)
            write_g4bl_bisector(case, d, gauss=False, finer=True)
        if bis:
            places = bisector_layout(case)
            maps = [dict(file=i["file"], rot=i.get("rot", ""), scale=i["scale"],
                         cl_z=i["cl_z"]) for i in case["sequence"] if i["kind"] == "map"]
            mons = [int(i["cl_z"]) for i in case["sequence"] if i["kind"] == "mon"]
            grids, boxes = bisector_grid(case, places), bisector_boxes(case, places)
        else:
            maps, mons = case["maps"], case["monitors"]
            grids, boxes = sample_grids(case), lab_boxes(case)
        manifest.append(dict(
            name=case["name"], desc=case["desc"], dir=case["name"], bisector=bis,
            maps=maps, monitors=mons, zstop=case["zstop"],
            n_pair=len(pset), n_gauss=N_GAUSS,
            dt_fine=DT_FINE, dt_coarse=DT_COARSE, dt_finer=DT_FINER,
            maxstep_fine=MAXSTEP_FINE, maxstep_coarse=MAXSTEP_COARSE,
            maxstep_finer=MAXSTEP_FINER, converge=bis,
            grids=grids, map_header=map_header(maps[0]["file"]), boxes=boxes,
            zero_field=all(m["file"] == "asr62shim_280_sm_track.g4blmap" for m in maps),
        ))
        print(f"{case['name']:22s} {len(maps)} map(s), {len(mons)} planes, "
              f"{sum(g['x'][2] * g['y'][2] * g['z'][2] for g in grids.values())} "
              f"field samples" + ("   [muE4 cornerarc frame]" if bis else ""))
    (HERE / "cases.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nwrote cases.json with {len(manifest)} cases")


if __name__ == "__main__":
    main()
