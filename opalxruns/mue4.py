#!/usr/bin/env python
"""mue4.py -- read a G4beamline deck and work out where every element ends up.

G4beamline places elements in *centreline* coordinates: a running coordinate whose z is
path length along the beamline and whose orientation is changed by `cornerarc` commands.
OPALX's FIELDMAP element is placed by an absolute lab pose instead, so reproducing the muE4
channel means replaying that centreline walk and reporting a lab position and orientation
for every element.

Two facts make this exact rather than approximate.

1. `place` resolves against whatever centreline segment is current when the line is read.
   `currentCL` moves only when a `cornerarc` line is executed (BLCoordinates.cc:358), so
   placement is purely deck-order driven and there is no z-based tie-break. The muE4 dipole
   maps sit between two `cornerarc` commands, so they are in the frame after the first one
   and before the second -- exactly the bisector of the total bend.

2. A `cornerarc` is three kinks, not an arc, but its endpoint is exactly the endpoint of the
   true arc. The two straight pieces each have length R*theta/2 and run at theta/2 -+ beta,
   so their sum is 2*R*(theta/2)*sin(theta/2)*cos(beta), and beta is chosen so that
   cos(beta) = sin(theta/2)/(theta/2). That collapses to R*(1-cos theta) across and
   R*sin(theta) along -- the chord of the true arc, for any angle. Only the middle vertex is
   off the arc, and nothing in muE4 is placed there.

So a true-arc transform reproduces G4beamline's own placement, which `test_geometry.py`
checks against the positions g4bl prints in `g4bl.out`.

Coordinates here are millimetres and degrees, matching the G4beamline deck. Conversion to
metres happens only when the OPALX deck is written.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

# The deck's z values are rounded to four decimals, so a placement at the end of a cornerarc
# can miss the computed boundary by tens of nanometres. Anything inside this window counts as
# being on the boundary; anything outside it is genuinely mid-arc and is an error.
SNAP_MM = 1.0e-3


# --------------------------------------------------------------------------------------
# Expression evaluation
# --------------------------------------------------------------------------------------
def evaluate(expr: str, params: dict[str, float]) -> float:
    """Evaluate a G4beamline attribute value.

    The deck uses plain arithmetic with $variables, e.g. `1285-0.5*(1839+29)` and
    `0.5*(190+800)`. Substitute the variables and evaluate; refuse anything containing a name,
    so a typo cannot silently reach through to a builtin.
    """
    text = expr.strip()
    for name, value in sorted(params.items(), key=lambda kv: -len(kv[0])):
        text = text.replace(f"${name}", f"({value!r})")
    # Blank out well-formed numbers, exponents included, before looking for stray names --
    # otherwise the E in `-9.9528E-01` reads as an identifier.
    residue = re.sub(r"\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?", "", text)
    if re.search(r"[A-Za-z_$]", residue):
        raise ValueError(f"unresolved name in expression {expr!r} -> {text!r}")
    return float(eval(text, {"__builtins__": {}}, {}))  # noqa: S307 - digits and operators only


def split_attrs(rest: str) -> list[tuple[str, str]]:
    """Split `k=v k=v ...` into pairs, keeping order and tolerating no spaces around `=`."""
    return re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(\S+)", rest)


# --------------------------------------------------------------------------------------
# Rotations
# --------------------------------------------------------------------------------------
def rot_y(a: float) -> list[list[float]]:
    c, s = math.cos(a), math.sin(a)
    return [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]


def rot_x(a: float) -> list[list[float]]:
    c, s = math.cos(a), math.sin(a)
    return [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]


def rot_z(a: float) -> list[list[float]]:
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def parse_rotation(text: str) -> list[list[float]]:
    """G4beamline `rotation=` string to a matrix.

    Successive rotations of the object about the FIXED axes of the parent, applied left to
    right: "Y180,Z180" means Rz(180) * Ry(180) (BLCommand::stringToRotationMatrix, which
    accumulates `*rot = thisRotation * *rot`).
    """
    result = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    if not text:
        return result
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        axis, angle = token[0].upper(), math.radians(float(token[1:]))
        this = {"X": rot_x, "Y": rot_y, "Z": rot_z}[axis](angle)
        result = matmul(this, result)
    return result


def to_tait_bryan(m: list[list[float]]) -> tuple[float, float, float]:
    """Decompose a rotation into OPALX's (THETA, PHI, PSI), in radians.

    OpalElement::update builds the pose as Ry(THETA) * Rx(PHI) * Rz(PSI), so this inverts
    that product. The muE4 cases are all either a pure y rotation or a y rotation composed
    with a single 180 degree flip, which this covers exactly; the general branch is here so a
    future deck cannot silently get a wrong answer.
    """
    # Ry(t) Rx(p) Rz(s) has m[1][2] = -sin(phi)
    sin_phi = max(-1.0, min(1.0, -m[1][2]))
    phi = math.asin(sin_phi)
    if abs(math.cos(phi)) > 1.0e-9:
        theta = math.atan2(m[0][2], m[2][2])
        psi = math.atan2(m[1][0], m[1][1])
    else:  # gimbal lock: y and z rotations are degenerate, put it all in theta
        theta = math.atan2(-m[2][0], m[0][0])
        psi = 0.0
    # The branch above is one of several equivalent triples. Prove it reproduces the input
    # rather than reasoning about which branch is right: a wrong pose is otherwise invisible
    # until the fields disagree.
    check = matmul(rot_y(theta), matmul(rot_x(phi), rot_z(psi)))
    worst = max(abs(check[i][j] - m[i][j]) for i in range(3) for j in range(3))
    if worst > 1.0e-12:
        raise ValueError(f"rotation decomposition failed, worst element error {worst:.3e}")
    return theta, phi, psi


# --------------------------------------------------------------------------------------
# Deck model
# --------------------------------------------------------------------------------------
@dataclass
class Solid:
    """A `box`, `tubs`, `virtualdetector` or `fieldmap` definition."""

    kind: str
    name: str
    attrs: dict[str, str] = field(default_factory=dict)


@dataclass
class Placement:
    """One `place` line, resolved into the lab frame."""

    index: int
    solid: Solid
    name: str  # after any rename=
    cl_z: float  # centreline z [mm]
    off_x: float  # transverse offset in the centreline frame [mm]
    off_y: float
    rotation: str  # the raw rotation= string
    current: float | None  # deck-side current=, the field scale
    line: int  # source line in the .g4bl file

    # filled in by walk()
    lab: tuple[float, float, float] = (0.0, 0.0, 0.0)  # [mm]
    frame_deg: float = 0.0  # accumulated centreline rotation about y
    pose: tuple[float, float, float] = (0.0, 0.0, 0.0)  # THETA, PHI, PSI [rad]


@dataclass
class CornerArc:
    cl_z: float
    angle_deg: float
    radius: float
    line: int


def parse_deck(path: Path) -> tuple[list[Placement], list[CornerArc], dict[str, Solid]]:
    """Read a G4beamline deck into placements, cornerarcs and solid definitions."""
    params: dict[str, float] = {}
    solids: dict[str, Solid] = {}
    placements: list[Placement] = []
    arcs: list[CornerArc] = []

    definition_kinds = {"box", "tubs", "virtualdetector", "fieldmap", "cylinder"}

    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        text = raw.split("#", 1)[0].strip()
        if not text:
            continue
        head, _, rest = text.partition(" ")

        if head == "param":
            for key, value in split_attrs(rest.replace("-unset", "")):
                try:
                    params[key] = evaluate(value, params)
                except ValueError:
                    pass  # non-numeric params (worldMaterial=Vacuum) are not needed here
            continue

        if head in definition_kinds:
            name, _, attrs = rest.partition(" ")
            solids[name] = Solid(head, name, dict(split_attrs(attrs)))
            continue

        if head == "cornerarc":
            a = dict(split_attrs(rest))
            arcs.append(
                CornerArc(
                    cl_z=evaluate(a["z"], params),
                    angle_deg=evaluate(a["angle"], params),
                    radius=evaluate(a["centerRadius"], params),
                    line=lineno,
                )
            )
            continue

        if head == "place":
            name, _, attrs = rest.partition(" ")
            a = dict(split_attrs(attrs))
            if name not in solids:
                raise KeyError(f"line {lineno}: place of undefined '{name}'")
            placements.append(
                Placement(
                    index=len(placements),
                    solid=solids[name],
                    name=a.get("rename", name),
                    cl_z=evaluate(a["z"], params) if "z" in a else 0.0,
                    off_x=evaluate(a["x"], params) if "x" in a else 0.0,
                    off_y=evaluate(a["y"], params) if "y" in a else 0.0,
                    rotation=a.get("rotation", ""),
                    current=evaluate(a["current"], params) if "current" in a else None,
                    line=lineno,
                )
            )
            continue

    return placements, arcs, solids


def walk(placements: list[Placement], arcs: list[CornerArc]) -> list[Placement]:
    """Resolve every placement into the lab frame.

    Walks the deck in order. Each `cornerarc` advances the running frame by its full angle and
    moves the anchor to the arc's endpoint; each `place` is resolved in whatever frame is
    current at that point in the file, which is what G4beamline itself does.
    """
    # Interleave by source line so deck order is respected.
    events = sorted(
        [("place", p) for p in placements] + [("arc", a) for a in arcs],
        key=lambda e: e[1].line,
    )

    px, pz = 0.0, 0.0  # anchor position [mm]
    angle = 0.0  # accumulated frame rotation about y [rad]
    z_anchor = 0.0  # centreline z the anchor corresponds to [mm]

    for kind, item in events:
        if kind == "arc":
            # Advance along the straight run to the arc's start.
            px += (item.cl_z - z_anchor) * math.sin(angle)
            pz += (item.cl_z - z_anchor) * math.cos(angle)
            # Chord of the arc, in the pre-arc frame, then rotated into the lab.
            t = math.radians(abs(item.angle_deg))
            sign = 1.0 if item.angle_deg > 0 else -1.0
            cx = item.radius * (1.0 - math.cos(t)) * sign
            cz = item.radius * math.sin(t)
            px, pz = (
                px + cx * math.cos(angle) + cz * math.sin(angle),
                pz - cx * math.sin(angle) + cz * math.cos(angle),
            )
            angle += math.radians(item.angle_deg)
            z_anchor = item.cl_z + item.radius * t
            continue

        # A placement. Guard against anything genuinely inside an arc, where a chord
        # transform would be wrong.
        ds = item.cl_z - z_anchor
        if ds < -SNAP_MM:
            raise ValueError(
                f"line {item.line}: '{item.name}' at centreline z={item.cl_z} is inside a "
                f"cornerarc (anchor {z_anchor}); the chord transform does not apply there"
            )
        ds = max(ds, 0.0)

        along = (math.sin(angle), math.cos(angle))  # local +z in the lab
        across = (math.cos(angle), -math.sin(angle))  # local +x in the lab
        item.lab = (
            px + ds * along[0] + item.off_x * across[0],
            item.off_y,
            pz + ds * along[1] + item.off_x * across[1],
        )
        item.frame_deg = math.degrees(angle)
        item.pose = to_tait_bryan(matmul(rot_y(angle), parse_rotation(item.rotation)))

    return placements


def load(deck: Path) -> tuple[list[Placement], dict[str, Solid]]:
    """Parse and walk a deck in one call."""
    placements, arcs, solids = parse_deck(deck)
    return walk(placements, arcs), solids
