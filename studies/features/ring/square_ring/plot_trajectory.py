#!/usr/bin/env python3
"""Plot the reference-orbit path of the square ring and print the closure offset per turn.

Reads data/square_ring_DesignPath.dat (written by the OrbitThreader) from the run
folder in output/. Columns:
s, Rx, Ry, Rz, Px, Py, Pz, Efx, Efy, Efz, Bfx, Bfy, Bfz, Ekin, t, element names.
The lab-frame x-z path should trace the square ring once per turn (circumference 8 m).
"""

import sys
from pathlib import Path

import numpy as np

from opalxruns.paths import output_dir

OUT = output_dir(Path(__file__).resolve().parent)     # the run, in output/

CIRCUMFERENCE = 8.0


def load_design_path(fname):
    rows = []
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.split()
            if len(tokens) < 15:
                continue
            try:
                rows.append([float(t) for t in tokens[:15]])
            except ValueError:
                continue
    data = np.array(rows)
    data = data[np.argsort(data[:, 0])]  # trackback rows are logged out of order
    return data


def main():
    fname = OUT / "data" / "square_ring_DesignPath.dat"
    if not fname.exists():
        sys.exit(f"not found: {fname} (run it first: python -m opalxruns.run square_ring.in)")

    data = load_design_path(fname)
    s, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    print(f"path length range: {s[0]:.3f} .. {s[-1]:.3f} m "
          f"({s[-1] / CIRCUMFERENCE:.2f} turns)")

    # closure: distance between the position at s = k * circumference and the start
    start = np.array([np.interp(0.0, s, x), np.interp(0.0, s, y), np.interp(0.0, s, z)])
    turn = 1
    while turn * CIRCUMFERENCE <= s[-1]:
        sk = turn * CIRCUMFERENCE
        pos = np.array([np.interp(sk, s, x), np.interp(sk, s, y), np.interp(sk, s, z)])
        offset = np.linalg.norm(pos - start)
        print(f"turn {turn}: position offset to start = {offset * 1e3:.3f} mm")
        turn += 1

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(z, x, lw=0.8)
    ax.plot(z[0], x[0], "go", label="start")
    ax.plot(z[-1], x[-1], "rx", label="end")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("x [m]")
    ax.set_title("square ring, reference-orbit path (lab frame)")
    ax.set_aspect("equal")
    ax.legend()
    out = OUT / "trajectory.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
