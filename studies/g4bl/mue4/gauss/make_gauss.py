#!/usr/bin/env python
"""make_gauss.py -- one Gaussian beam, written for both codes.

The whole point is that both codes get the *same* particles: the distribution is sampled once,
with a fixed seed, and then written twice in the two formats. Nothing is sampled per code.

  parts.txt   OPALX FROMFILE     x y z px py pz, positions m, momenta beta*gamma
  beam.txt    G4beamline         #BLTrackFile, mm and MeV/c

Defaults come from the muE4 deck's own `beam gaussian` line (commented out there in favour of
the phase-space file, but it is the deck author's statement of the beam):

    sigmaX = sigmaY = 0, sigmaXp = sigmaYp = 0.130 rad, sigmaP = 0, P = 28 MeV/c

130 mrad is wide, so this beam does hit the collimators -- which is the point. The hand-picked
set in `cmp/` was chosen to miss everything and test the optics alone; this one tests what the
channel actually does to a beam.

A note on the reference momentum. OPALX takes the reference from the mean of the particle
momenta and cannot be told otherwise under FROMFILE. For a symmetric Gaussian that mean is
close to the design momentum -- the only bias is the cos(theta) foreshortening of the forward
component, about -0.85 % at 130 mrad -- so this beam does not suffer the problem that makes the
real muE4 phase-space file unusable, where the mean forward momentum is 15.5 MeV/c against a
28 MeV/c design.

Usage:  python3 make_gauss.py [-n 2000] [--sigma-xp 0.130] [--seed 12345]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cmp"))
import make_case  # noqa: E402  -- reuse the twin generator

MUON_MASS = 105.6583755  # MeV


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=2000, help="number of particles")
    ap.add_argument("--p0", type=float, default=28.0, help="momentum [MeV/c]")
    ap.add_argument("--sigma-x", type=float, default=0.0, help="[mm]")
    ap.add_argument("--sigma-y", type=float, default=0.0, help="[mm]")
    ap.add_argument("--sigma-xp", type=float, default=0.130, help="[rad]")
    ap.add_argument("--sigma-yp", type=float, default=0.130, help="[rad]")
    ap.add_argument("--sigma-p", type=float, default=0.0, help="relative momentum spread")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n = args.n

    x = rng.normal(0.0, args.sigma_x, n) if args.sigma_x > 0 else np.zeros(n)
    y = rng.normal(0.0, args.sigma_y, n) if args.sigma_y > 0 else np.zeros(n)
    xp = rng.normal(0.0, args.sigma_xp, n)
    yp = rng.normal(0.0, args.sigma_yp, n)
    p = args.p0 * (1.0 + rng.normal(0.0, args.sigma_p, n)) if args.sigma_p > 0 else np.full(n, args.p0)

    # xp, yp are the transverse direction tangents, so the direction vector is
    # (xp, yp, 1) normalised, and |p| is the sampled momentum.
    norm = np.sqrt(1.0 + xp**2 + yp**2)
    px, py, pz = p * xp / norm, p * yp / norm, p / norm

    (HERE / "parts.txt").write_text(
        f"{n}\nx y z px py pz\n"
        + "".join(
            f"{x[i] * 1e-3:.9e} {y[i] * 1e-3:.9e} 0.0 "
            f"{px[i] / MUON_MASS:.9e} {py[i] / MUON_MASS:.9e} {pz[i] / MUON_MASS:.9e}\n"
            for i in range(n)
        )
    )

    (HERE / "beam.txt").write_text(
        "#BLTrackFile written by make_gauss.py\n"
        "#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n"
        "#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n"
        + "".join(
            f"{x[i]:.9g} {y[i]:.9g} 0 {px[i]:.9g} {py[i]:.9g} {pz[i]:.9g} 0 -13 {i + 1} 1 0 1\n"
            for i in range(n)
        )
    )

    # The G4beamline twin: identical geometry to mue4_WsxOn.g4bl, reading beam.txt.
    make_case.HERE = HERE
    make_case.write_twin()

    mean_pz = float(pz.mean())
    print(f"wrote parts.txt and beam.txt: {n} muons, seed {args.seed}")
    print(f"  sigma_xp = sigma_yp = {args.sigma_xp} rad, |p| = {args.p0} MeV/c")
    print(f"  mean forward momentum {mean_pz:.3f} MeV/c "
          f"({100 * (mean_pz / args.p0 - 1):+.2f} % vs design -- this is what OPALX will use "
          f"as its reference)")
    print(f"  rms angle {math.hypot(xp.std(), yp.std()) * 1000:.1f} mrad")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
