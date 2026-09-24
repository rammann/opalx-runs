#!/usr/bin/env python
"""particles_to_vtk.py -- OPALX phase-space ``.h5`` dumps to ParaView VTP/PVD.

Writes the bunch in either or both of two frames:

``comoving``
    The raw ``.h5`` contents. OPALX tracks and dumps in the co-moving reference
    frame: the origin sits on the reference particle and local +z points along
    the reference momentum. Nothing is transformed, so the bunch is drawn in its
    own straightened frame -- good for watching the phase space evolve, wrong if
    you want the bunch threaded along a bent orbit.

``lab``
    Absolute lab coordinates, so the bunch lines up with the beamline elements
    written by ``elements_to_vtk.py``::

        R_lab = RefPartR(s) + Q(s) . (x, y, z)_local
        P_lab =               Q(s) . (px, py, pz)_local

    ``RefPartR(s)`` is stored per step. ``Q(s)``, the co-moving -> lab rotation,
    is **not**: the ``TaitBryantAngles`` attribute is written but stubbed to
    ``[0,0,0]`` in ``H5PartWrapperForPT.cpp``. It is rebuilt by parallel
    transport of the reference-momentum direction taken per step from
    ``data/<base>_DesignPath.dat``, which reproduces OPALX's own ``toLabTrafo``
    -- see :func:`opalx_run.build_transport_frames`.

    Without a DesignPath the script falls back to a single shortest-arc rotation
    from +z onto the final reference momentum. That misses the accumulated frame
    roll (~24 deg by the end of muE4) and can misplace wide-beam tail particles
    by ~1 m, so it is only valid for narrow, near-planar beams.

Usage
-----
    python particles_to_vtk.py <run_dir>                    # both frames
    python particles_to_vtk.py <run_dir> --frame lab
    python particles_to_vtk.py <run_dir> --stride 5 --mass 0.000511

Open the ``.pvd`` in ParaView (not the individual ``.vtp``) to get the time
series. Colour by ``Ekin_MeV`` or ``Pmag``; add a Glyph filter on the vector
``P`` for momentum arrows.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from vtk import vtkCellArray, vtkPoints, vtkPolyData, vtkXMLPolyDataWriter
from vtk.util.numpy_support import numpy_to_vtk

from opalxruns.opalx_run import (
    ZHAT,
    Run,
    build_transport_frames,
    h5_steps,
    load_designpath,
    shortest_arc,
)

MUON_MASS_GEV = 0.1056583755

# Datasets making up a vector in the co-moving frame: they get assembled into a
# single VTK vector array, and rotated by Q(s) along with the position in the lab
# frame. Position and momentum are handled separately.
VECTOR_TRIPLES = {
    "Pol": ("polx", "poly", "polz"),   # spin polarization
    "E": ("Ex", "Ey", "Ez"),           # field seen by the particle
    "B": ("Bx", "By", "Bz"),
}
_TRIPLE_MEMBERS = {n for t in VECTOR_TRIPLES.values() for n in t}
_HANDLED = {"x", "y", "z", "px", "py", "pz"} | _TRIPLE_MEMBERS


def _to_1d(group, name):
    return np.asarray(group[name][()]).reshape(-1)


def _extra_arrays(group, n):
    """Everything else in the step, as ``(name, values, is_vector)``.

    Recognised triples (polarization, E, B) come back as ``(n, 3)`` vectors;
    every other 1-D per-particle dataset is passed through as a scalar, so a run
    carrying datasets this script has never heard of still exports them.
    """
    out = []
    for name, cols in VECTOR_TRIPLES.items():
        if all(c in group for c in cols):
            vec = np.column_stack([_to_1d(group, c) for c in cols])
            if len(vec) == n:
                out.append((name, vec, True))
    for name in sorted(group):
        if name in _HANDLED:
            continue
        ds = group[name]
        if getattr(ds, "ndim", 0) != 1:
            continue
        vals = _to_1d(group, name)
        if len(vals) == n:
            out.append((name, vals.astype(np.float64), False))
    return out


def _add_array(polydata, name, values):
    arr = numpy_to_vtk(np.ascontiguousarray(values, dtype=np.float64), deep=1)
    arr.SetName(name)
    polydata.GetPointData().AddArray(arr)


def _write_step_vtp(group, rot, origin, mass_gev, out_file):
    """Write one ``.vtp``. ``rot`` is ``None`` for the untransformed frame.

    Returns ``(n_particles, time)``.
    """
    r_local = np.column_stack([_to_1d(group, n) for n in ("x", "y", "z")])
    p_local = np.column_stack([_to_1d(group, n) for n in ("px", "py", "pz")])
    n = len(r_local)

    if rot is None:
        r_out, p_out = r_local, p_local
    else:
        r_out = origin + r_local @ rot.T
        p_out = p_local @ rot.T

    points = vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(r_out, dtype=np.float64), deep=1))

    verts = vtkCellArray()
    for i in range(n):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)

    poly = vtkPolyData()
    poly.SetPoints(points)
    poly.SetVerts(verts)

    _add_array(poly, "P", p_out)                  # momentum vector, for glyphs
    bg = np.linalg.norm(p_local, axis=1)          # |p| is frame-invariant
    _add_array(poly, "Pmag", bg)                  # in beta*gamma
    _add_array(poly, "Ekin_MeV", mass_gev * 1e3 * (np.sqrt(1.0 + bg * bg) - 1.0))

    # Polarization / E / B are vectors in the same co-moving frame as the
    # momentum, so they rotate with it; scalars are copied as they are.
    for name, values, is_vector in _extra_arrays(group, n):
        if is_vector and rot is not None:
            values = values @ rot.T
        _add_array(poly, name, values)

    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(str(out_file))
    writer.SetInputData(poly)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise RuntimeError(f"failed writing {out_file}")

    t = group.attrs.get("TIME")
    t = float(np.asarray(t).reshape(-1)[0]) if t is not None else 0.0
    return n, t


def _write_pvd(pvd_path, entries):
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
             "  <Collection>"]
    for t, fname in entries:
        lines.append(f'    <DataSet timestep="{t:.16g}" part="0" file="{fname}"/>')
    lines += ["  </Collection>", "</VTKFile>", ""]
    pvd_path.write_text("\n".join(lines), encoding="utf-8")


def _lab_frames(h5_path, designpath):
    """Return ``(s_samples, Q_samples)`` for the lab transform, or ``None``.

    ``None`` means no DesignPath was available and the caller should fall back
    to a per-step shortest-arc rotation.
    """
    if designpath is None or not Path(designpath).exists():
        print("  WARNING: no DesignPath -> momentum-only rotation "
              "(misses frame roll; OK only for narrow near-planar beams)")
        return None
    s_dp, _, p_dp = load_designpath(designpath)
    print(f"  frame: parallel transport from {Path(designpath).name} "
          f"({len(s_dp)} steps)")
    return s_dp, build_transport_frames(p_dp)


def convert(h5_path, out_dir, frame, designpath=None, mass_gev=MUON_MASS_GEV,
            stride=1, prefix=None, quiet=False):
    """Convert one ``.h5`` to a ``.pvd`` + per-step ``.vtp`` series.

    Returns the path of the ``.pvd``.
    """
    h5_path = Path(h5_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = prefix or f"bunch_{frame}"

    transport = _lab_frames(h5_path, designpath) if frame == "lab" else None

    steps = h5_steps(h5_path)
    if not steps:
        raise SystemExit(f"no Step#N groups in {h5_path}")
    steps = steps[::stride]
    print(f"  {h5_path.name}: {len(steps)} steps -> {out_dir}/{prefix}_*.vtp")

    entries = []
    with h5py.File(h5_path, "r") as f:
        for key in steps:
            g = f[key]
            if frame == "lab":
                origin = np.asarray(g.attrs["RefPartR"]).reshape(-1)[:3]
                if transport is not None:
                    s_dp, Qs = transport
                    spos = float(np.asarray(g.attrs["SPOS"]).reshape(-1)[0])
                    rot = Qs[int(np.argmin(np.abs(s_dp - spos)))]
                else:
                    rot = shortest_arc(ZHAT, np.asarray(g.attrs["RefPartP"]).reshape(-1)[:3])
            else:
                origin, rot = None, None

            num = int(key[5:])
            fname = f"{prefix}_step{num:05d}.vtp"
            n, t = _write_step_vtp(g, rot, origin, mass_gev, out_dir / fname)
            entries.append((t, fname))
            if not quiet:
                print(f"    {fname}  ({n} particles, t={t:.4e} s)")

    pvd = out_dir / f"{prefix}.pvd"
    _write_pvd(pvd, entries)
    print(f"  wrote {pvd}  ({len(entries)} steps)")
    return pvd


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="OPALX run directory (default: cwd)")
    ap.add_argument("--frame", choices=("lab", "comoving", "both"), default="both",
                    help="which frame(s) to write (default: both)")
    ap.add_argument("--base", default=None,
                    help="run basename, if the directory holds more than one deck")
    ap.add_argument("--h5", default=None,
                    help="phase-space .h5 to convert (default: the run's <base>.h5)")
    ap.add_argument("--designpath", default=None,
                    help="path to <base>_DesignPath.dat (default: from the run dir)")
    ap.add_argument("--output-dir", default=None,
                    help="output directory (default: <run>/paraview)")
    ap.add_argument("--mass", type=float, default=MUON_MASS_GEV,
                    help=f"species rest mass in GeV for Ekin_MeV "
                         f"(default: muon {MUON_MASS_GEV})")
    ap.add_argument("--stride", type=int, default=1,
                    help="write every Nth dump (default: 1 = all)")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="do not list every step")
    args = ap.parse_args(argv)

    run = Run(args.run_dir, base=args.base)
    h5_path = Path(args.h5) if args.h5 else run.bunch_h5
    if h5_path is None:
        print(f"particles_to_vtk: no phase-space .h5 in {run.dir}", file=sys.stderr)
        return 1

    out_dir = Path(args.output_dir) if args.output_dir else run.paraview_dir
    designpath = args.designpath or run.design_path
    frames = ("lab", "comoving") if args.frame == "both" else (args.frame,)

    for frame in frames:
        convert(h5_path, out_dir, frame, designpath=designpath,
                mass_gev=args.mass, stride=args.stride, quiet=args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
