#!/usr/bin/env python
"""elements_to_vtk.py -- beamline ELEMENTS and the reference orbit as ParaView
geometry, in absolute LAB coordinates.

Load these together with the lab-frame bunch from ``particles_to_vtk.py``; they
share the same frame.

What it makes
-------------
``<prefix>_elements.vtp``
    One body per element, swept along the element's lab-frame centerline (bent
    dipoles follow their arc). An element WITH an aperture becomes a hollow
    pipe: the inner surface is the real aperture (rectangle or ellipse), the
    outer surface sits ``--wall`` further out, and the ends are closed with
    annular faces -- so the bore is visible and particles fly through the
    opening. An element WITHOUT an aperture (e.g. a hard-edge bend with
    HGAP = 0) has no bore to show and becomes a solid tube of radius
    ``--default-aperture``. Cell arrays: ``element_type`` (colour by type) and
    ``element_id`` (index into the printed name list).
``<prefix>_reforbit.vtp``
    The reference / design orbit as a polyline.

Where the geometry comes from (reuse, don't recompute)
------------------------------------------------------
OPALX already dumps, on every run, into ``<run>/data/``:

* ``<base>_ElementPositions.txt`` -- per-element BEGIN/MID/END (plus ENTRY/EXIT
  EDGE for bends) points, already in lab coordinates. That gives the bent dipole
  centerline for free, so no placement or bend math is redone here.
* ``<base>_DesignPath.dat`` -- the finely sampled reference orbit.

Apertures and element types are read from the input deck. The generated
``<base>_ElementPositions.py --export-vtk`` is deliberately not used: it omits
the dipoles.

Aperture convention (matches OPALX)
-----------------------------------
Same rules as ``OpalElement::getApert()`` / ``OpalSBend::update()``:

    APERTURE = "rectangle(a,b)" / "square(a)"  -> rectangle, half-widths (a/2, b/2)
    APERTURE = "ellipse(a,b)"  / "circle(a)"   -> ellipse,   half-widths (a/2, b/2)
    SBEND/RBEND: APERTURE if given; else HGAP > 0 -> rectangle, half-widths
        (HAPERT if given else unbounded, HGAP); else no aperture (hard edge)
    other elements without APERTURE            -> no aperture

HAPERT/HGAP are already half-values; deck APERTURE arguments are full sizes and
get halved. An unbounded half-width (OPALX default 1e6 m) is clamped to
``--default-aperture`` for display, with a note.

Usage
-----
    python elements_to_vtk.py <run_dir>
    python elements_to_vtk.py <run_dir> --default-aperture 0.05 --wall 0.02
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from vtk import (vtkCellArray, vtkPoints, vtkPolyData, vtkPolyLine, vtkTriangle,
                 vtkXMLPolyDataWriter)
from vtk.util.numpy_support import numpy_to_vtk

sys.path.insert(0, str(Path(__file__).resolve().parent))
from opalx_run import Run, load_designpath, parse_element_positions  # noqa: E402

# element type -> integer code (colour on `element_type` in ParaView)
TYPE_CODE = {"DRIFT": 0, "SOLENOID": 1, "SBEND": 2, "RBEND": 2,
             "QUADRUPOLE": 3, "MONITOR": 4}
TYPE_LEGEND = {0: "DRIFT", 1: "SOLENOID", 2: "DIPOLE", 3: "QUADRUPOLE",
               4: "MONITOR", 5: "OTHER"}


# --------------------------------------------------------------------------- #
# deck parsing: element -> (type, aperture)                                    #
# --------------------------------------------------------------------------- #
def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)   # block comments
    text = re.sub(r"//[^\n]*", "", text)                     # line comments
    return text


def parse_reals(text):
    """Numeric ``REAL name = value;`` assignments (for L_QUAD etc.)."""
    syms = {}
    for m in re.finditer(r"\bREAL\s+(\w+)\s*=\s*([^;]+);", text):
        try:
            syms[m.group(1)] = float(m.group(2))
        except ValueError:
            pass
    return syms


def num(token, syms):
    token = token.strip()
    try:
        return float(token)
    except ValueError:
        return syms.get(token, None)


def parse_aperture_string(s):
    """``'rectangle(0.380,0.210)'`` -> ``('rect', 0.190, 0.105)``."""
    m = re.search(r'"?\s*(square|rectangle|circle|ellipse)\s*\(([^)]*)\)', s, flags=re.IGNORECASE)
    if not m:
        return None
    shape = m.group(1).lower()
    args = [float(a) for a in re.split(r"[,\s]+", m.group(2).strip()) if a]
    if shape == "rectangle":
        return ("rect", args[0] / 2.0, args[1] / 2.0)
    if shape == "square":
        return ("rect", args[0] / 2.0, args[0] / 2.0)
    if shape == "ellipse":
        return ("ell", args[0] / 2.0, args[1] / 2.0)
    if shape == "circle":
        return ("ell", args[0] / 2.0, args[0] / 2.0)
    return None


UNBOUNDED = 1e5  # half-widths at/above this (OPALX default 1e6) count as "no bound"


def parse_deck(deck_path):
    """``{name: {"type", "aper"}}`` for every element the deck declares.

    ``aper`` is ``(shape, hx, hy)`` with OPALX semantics, or ``None`` when the
    element has no transverse aperture at all (rendered as a solid body).
    Unbounded half-widths are kept as-is here and clamped at render time.
    """
    text = strip_comments(Path(deck_path).read_text())
    syms = parse_reals(text)
    elements = {}
    for stmt in text.split(";"):
        m = re.match(r"\s*(\w+)\s*:\s*(\w+)\b(.*)", stmt, flags=re.DOTALL)
        if not m:
            continue
        name, typ, rest = m.group(1), m.group(2).upper(), m.group(3)
        if typ not in TYPE_CODE:
            continue
        am = re.search(r'\bAPERTURE\s*=\s*"?\s*(\w+\s*\([^)]*\))', rest)
        aper = parse_aperture_string(am.group(1)) if am else None
        if typ in ("SBEND", "RBEND") and aper is None:
            # OpalSBend::update(): with HGAP > 0 the poles bound y; HAPERT (if
            # given) bounds x. HGAP = 0 is a hard-edge bend with no aperture.
            hg = re.search(r"\bHGAP\s*=\s*([^\s,]+)", rest)
            hap = re.search(r"\bHAPERT\s*=\s*([^\s,]+)", rest)
            hy = num(hg.group(1), syms) if hg else None
            hx = num(hap.group(1), syms) if hap else None
            if hy is not None and hy > 0.0:
                aper = ("rect", hx if hx is not None else UNBOUNDED, hy)
        # OPALX uppercases object names; _ElementPositions.txt records them uppercase.
        elements[name.upper()] = {"type": typ, "aper": aper}
    return elements


# --------------------------------------------------------------------------- #
# tube meshing                                                                 #
# --------------------------------------------------------------------------- #
def frames_along(P):
    """Unit tangent, right and up at every polyline vertex (elements are planar)."""
    T = np.empty_like(P)
    T[1:-1] = P[2:] - P[:-2]
    T[0] = P[1] - P[0]
    T[-1] = P[-1] - P[-2]
    T /= np.linalg.norm(T, axis=1, keepdims=True)
    world_up = np.array([0.0, 1.0, 0.0])
    right = np.cross(np.broadcast_to(world_up, P.shape), T)
    bad = np.linalg.norm(right, axis=1) < 1e-8            # tangent nearly vertical
    if bad.any():
        right[bad] = np.cross(np.array([1.0, 0.0, 0.0]), T[bad])
    right /= np.linalg.norm(right, axis=1, keepdims=True)
    up = np.cross(T, right)
    up /= np.linalg.norm(up, axis=1, keepdims=True)
    return right, up


def cross_section(aper, n_ell):
    """``(a, b)`` offsets [m] in the (right, up) plane."""
    shape, hx, hy = aper
    if shape == "rect":
        return [(hx, hy), (-hx, hy), (-hx, -hy), (hx, -hy)]
    th = np.linspace(0.0, 2.0 * np.pi, n_ell, endpoint=False)
    return list(zip(hx * np.cos(th), hy * np.sin(th)))


def _sweep_rings(P, sec, verts):
    """Append one ring of the cross-section ``sec`` per path point; return ring starts."""
    right, up = frames_along(P)
    starts = []
    for i in range(len(P)):
        starts.append(len(verts))
        for (a, b) in sec:
            verts.append(P[i] + a * right[i] + b * up[i])
    return starts


def _connect_rings(starts, m, tris, codes, code, flip=False):
    """Triangulate the surface between consecutive rings of ``m`` points each."""
    for r0, r1 in zip(starts[:-1], starts[1:]):
        for j in range(m):
            k = (j + 1) % m
            quads = ((r0 + j, r0 + k, r1 + k), (r0 + j, r1 + k, r1 + j))
            for tri in quads:
                tris.append(tri[::-1] if flip else tri)
                codes.append(code)


def build_solid(P, aper, n_ell, verts, tris, codes, code):
    """One solid body (surface + end caps) — for elements with no aperture."""
    sec = cross_section(aper, n_ell)
    m = len(sec)
    starts = _sweep_rings(P, sec, verts)
    _connect_rings(starts, m, tris, codes, code)
    # end caps (triangle fan from each end-ring centroid)
    for ring_start, P_end in ((starts[0], P[0]), (starts[-1], P[-1])):
        c = len(verts)
        verts.append(P_end)
        for j in range(m):
            tris.append((c, ring_start + j, ring_start + (j + 1) % m))
            codes.append(code)


def build_pipe(P, aper, wall, n_ell, verts, tris, codes, code):
    """One hollow pipe: inner surface = the aperture, outer surface ``wall`` further
    out, ends closed with annular faces. The bore stays open."""
    shape, hx, hy = aper
    inner = cross_section(aper, n_ell)
    outer = cross_section((shape, hx + wall, hy + wall), n_ell)
    m = len(inner)
    si = _sweep_rings(P, inner, verts)
    so = _sweep_rings(P, outer, verts)
    _connect_rings(si, m, tris, codes, code, flip=True)   # inner wall faces the bore
    _connect_rings(so, m, tris, codes, code)
    # end annuli between inner and outer ring
    for ri, ro in ((si[0], so[0]), (si[-1], so[-1])):
        for j in range(m):
            k = (j + 1) % m
            for tri in ((ri + j, ro + j, ro + k), (ri + j, ro + k, ri + k)):
                tris.append(tri)
                codes.append(code)


def write_polydata_triangles(path, verts, tris, codes, ids):
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(verts, dtype=np.float64), deep=1))
    cells = vtkCellArray()
    for (a, b, c) in tris:
        t = vtkTriangle()
        t.GetPointIds().SetId(0, a)
        t.GetPointIds().SetId(1, b)
        t.GetPointIds().SetId(2, c)
        cells.InsertNextCell(t)
    poly = vtkPolyData()
    poly.SetPoints(pts)
    poly.SetPolys(cells)
    for name, values in (("element_type", codes), ("element_id", ids)):
        arr = numpy_to_vtk(np.ascontiguousarray(values, dtype=np.int32), deep=1)
        arr.SetName(name)
        poly.GetCellData().AddArray(arr)
    w = vtkXMLPolyDataWriter()
    w.SetFileName(str(path))
    w.SetInputData(poly)
    w.SetDataModeToBinary()
    if w.Write() != 1:
        raise RuntimeError(f"failed writing {path}")


def write_polyline(path, P):
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(P, dtype=np.float64), deep=1))
    line = vtkPolyLine()
    line.GetPointIds().SetNumberOfIds(len(P))
    for i in range(len(P)):
        line.GetPointIds().SetId(i, i)
    cells = vtkCellArray()
    cells.InsertNextCell(line)
    poly = vtkPolyData()
    poly.SetPoints(pts)
    poly.SetLines(cells)
    w = vtkXMLPolyDataWriter()
    w.SetFileName(str(path))
    w.SetInputData(poly)
    w.SetDataModeToBinary()
    if w.Write() != 1:
        raise RuntimeError(f"failed writing {path}")


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="OPALX run directory (default: cwd)")
    ap.add_argument("--base", default=None,
                    help="run basename, if the directory holds more than one deck")
    ap.add_argument("--deck", default=None, help="input .in deck (default: the run's)")
    ap.add_argument("--output-dir", default=None,
                    help="output directory (default: <run>/paraview)")
    ap.add_argument("--prefix", default=None, help="output prefix (default: run basename)")
    ap.add_argument("--default-aperture", type=float, default=0.05,
                    help="display half-width [m] for solid bodies of elements without an "
                         "aperture, and clamp for unbounded half-widths (default: 0.05)")
    ap.add_argument("--wall", type=float, default=0.01,
                    help="pipe wall thickness [m] for elements with an aperture "
                         "(default: 0.01)")
    ap.add_argument("--ellipse-points", type=int, default=24,
                    help="samples around elliptical apertures (default: 24)")
    ap.add_argument("--orbit-stride", type=int, default=10,
                    help="subsample stride for the reference orbit (default: 10)")
    args = ap.parse_args(argv)

    run = Run(args.run_dir, base=args.base)
    elpos = run.element_positions
    if elpos is None:
        print(f"elements_to_vtk: no *_ElementPositions.txt under {run.dir}",
              file=sys.stderr)
        return 1
    deck = Path(args.deck) if args.deck else run.deck
    if deck is None:
        print(f"elements_to_vtk: no input deck in {run.dir}", file=sys.stderr)
        return 1

    out_dir = Path(args.output_dir) if args.output_dir else run.paraview_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or run.base

    elements = parse_deck(deck)
    polylines = parse_element_positions(elpos)
    print(f"  deck: {deck.name}  ({len(elements)} elements)")
    print(f"  positions: {elpos.name}  ({len(polylines)} placed elements)")

    verts, tris, codes, ids = [], [], [], []
    counts = {}
    names = []
    for name, P in polylines:
        info = elements.get(name.upper(), {"type": "OTHER", "aper": None})
        code = TYPE_CODE.get(info["type"], 5)
        if len(P) < 2:
            continue
        eid = len(names)
        names.append(name)
        n_before = len(tris)
        aper = info["aper"]
        if aper is None:
            # no aperture (hard-edge bend, bare element): solid display body
            solid = ("ell", args.default_aperture, args.default_aperture)
            build_solid(P, solid, args.ellipse_points, verts, tris, codes, code)
            kind = "solid"
        else:
            shape, hx, hy = aper
            clamped = []
            if hx >= UNBOUNDED:
                hx, clamped = args.default_aperture, clamped + ["x"]
            if hy >= UNBOUNDED:
                hy, clamped = args.default_aperture, clamped + ["y"]
            if clamped:
                print(f"  note: {name}: unbounded {'/'.join(clamped)} half-width "
                      f"clamped to {args.default_aperture} m for display")
            build_pipe(P, (shape, hx, hy), args.wall, args.ellipse_points,
                       verts, tris, codes, code)
            kind = f"{shape} {2*hx:.3g} x {2*hy:.3g} m bore"
        ids.extend([eid] * (len(tris) - n_before))
        counts[TYPE_LEGEND[code]] = counts.get(TYPE_LEGEND[code], 0) + 1
        print(f"    [{eid:2d}] {name:12s} {info['type']:10s} {kind}")

    if not tris:
        print("elements_to_vtk: no element geometry to write", file=sys.stderr)
        return 1

    el_out = out_dir / f"{prefix}_elements.vtp"
    write_polydata_triangles(el_out, np.array(verts), tris, codes, ids)
    print(f"  wrote {el_out}  ({len(tris)} triangles)  by type: "
          + ", ".join(f"{k}:{v}" for k, v in sorted(counts.items())))
    print("  colour on `element_type`  "
          + "  ".join(f"{c}={TYPE_LEGEND[c]}" for c in sorted(TYPE_LEGEND))
          + "   (or on `element_id`, indices above)")

    if run.design_path is not None:
        _, R, _ = load_designpath(run.design_path)
        if args.orbit_stride > 1:
            R = R[::args.orbit_stride]
        orb_out = out_dir / f"{prefix}_reforbit.vtp"
        write_polyline(orb_out, R)
        print(f"  wrote {orb_out}  ({len(R)} points)")
    else:
        print("  note: no *_DesignPath.dat -> skipped the reference orbit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
