#!/usr/bin/env python3
"""Compare asr61_dipole.in (OPALX) against asr61_dipole.g4bl, particle by particle.

Run both codes first, from this directory:
    python3 make_parts.py
    g4bl asr61_dipole.g4bl                       > g4bl.log 2>&1
    mpirun -n 1 <opalx> asr61_dipole.in --info 1 > run.log  2>&1
    python3 compare.py

Needs h5py. The system python3 has none; use the conda base env,
/opt/homebrew/Caskroom/miniconda/base/bin/python. Without h5py this falls back
to h5dump, which is enough for the numbers below.

OPALX rows in the monitor file are NOT in id order -- they come out in whatever
order the ranks hold them. Everything is matched on the `id` dataset; comparing
row 0 to row 0 silently pairs different particles and manufactures a difference
of tens of mm.
"""

import re
import subprocess
import sys

MUON_MASS_MEV = 105.6583755


def h5_column(path, dataset):
    """One dataset of Step#0 as a list of floats, via h5py or h5dump."""
    try:
        import h5py

        with h5py.File(path, "r") as f:
            return list(f["Step#0"][dataset][:])
    except ImportError:
        out = subprocess.run(
            ["h5dump", "-d", f"/Step#0/{dataset}", path], capture_output=True, text=True
        ).stdout
        body = re.search(r"DATA \{\n(.*?)\n   \}", out, re.S).group(1)
        return [float(t) for t in re.sub(r"\(\d+\):", "", body).replace(",", " ").split()]


def opalx_plane(path):
    """{id: (x, y, Px, Py, Pz)} in mm and MeV/c. OPALX momenta are beta*gamma."""
    ids = [int(i) for i in h5_column(path, "id")]
    cols = [h5_column(path, c) for c in ("x", "y", "px", "py", "pz")]
    return {
        i: (x * 1e3, y * 1e3, px * MUON_MASS_MEV, py * MUON_MASS_MEV, pz * MUON_MASS_MEV)
        for i, x, y, px, py, pz in zip(ids, *cols)
    }


def g4bl_plane(path):
    """{id: (x, y, Px, Py, Pz)} in mm and MeV/c. EventID == OPALX id + 1."""
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) < 9:
            continue
        out[int(p[8]) - 1] = (float(p[0]), float(p[1]), float(p[3]), float(p[4]), float(p[5]))
    return out


def report(name, opalx, g4bl):
    shared = sorted(set(opalx) & set(g4bl))
    print(f"\n=== {name}  ({len(shared)} particles in both) ===")
    print(f"{'id':>3} {'dx [um]':>10} {'dy [um]':>10} {'dPx [keV]':>11} {'dPz [keV]':>11}")
    worst = 0.0
    for i in shared:
        a, b = opalx[i], g4bl[i]
        row = ((a[0] - b[0]) * 1e3, (a[1] - b[1]) * 1e3, (a[2] - b[2]) * 1e3, (a[4] - b[4]) * 1e3)
        worst = max(worst, abs(row[0]))
        print(f"{i:>3} {row[0]:>10.2f} {row[1]:>10.2f} {row[2]:>11.2f} {row[3]:>11.2f}")
    print(f"worst |dx| = {worst:.2f} um")
    return worst


def field_check():
    """OPALX cannot dump the map on its own, so this only reports g4bl's axis file.

    The field import itself is covered by the unit tests
    (TestG4BL3DMagnetoStatic); this is here so the axis profile can be eyeballed
    next to the orbit.
    """
    try:
        rows = [l.split() for l in open("g4bl_axis.txt") if not l.startswith("#")]
    except FileNotFoundError:
        print("\ng4bl_axis.txt not found -- run g4bl first")
        return
    peak = max(rows, key=lambda r: abs(float(r[5])))
    print(f"\ng4bl on-axis peak By = {float(peak[5]):.6g} T at global z = {float(peak[2]):.0f} mm")


if __name__ == "__main__":
    try:
        report("entrance plane, z = 50 mm", opalx_plane("MON_IN.h5"), g4bl_plane("Z50.txt"))
        report("exit plane, z = 2950 mm", opalx_plane("MON_OUT.h5"), g4bl_plane("Z2950.txt"))
    except FileNotFoundError as error:
        sys.exit(f"missing output: {error}. Run both codes first, see the module docstring.")
    field_check()
