#!/usr/bin/env python
"""compare.py -- OPALX against G4beamline for the Gaussian beam.

Two questions are asked separately, because they have different answers.

**Does the tracking agree?** For every particle that reaches a plane in both codes, compare
its transverse position there. This is the field-and-geometry test, and it is unaffected by
anything the collimators do. Same method as `../cmp/compare.py`: centreline coordinates so
G4beamline's 6-significant-figure output resolves ~0.01 um rather than 36 um, and each
G4beamline record drifted onto the plane because a virtualdetector is a 1 mm thick volume
recording on entry.

**Does the transmission agree?** Only partly, and the reason is a known modelling limit rather
than a physics difference. A G4beamline collimator jaw is a finite iron block -- FS61H spans
|x| in [190, 800] mm and |y| < 200 mm -- so a particle at |x| > 800 mm misses it entirely. An
OPALX aperture is a single inside/outside test, so `rectangle(0.380, ...)` kills every particle
with |x| >= 190 whatever its y or how far out it is. With a 130 mrad beam that is a large
population. Measured directly, at the last detector: OPALX keeps 79 with the apertures as
modelled, 152 with them opened, and G4beamline keeps 122 -- i.e. G4beamline sits between the
two, exactly as an over-aggressive aperture model predicts.

A G4beamline virtualdetector is also a finite disc that only records what hits it, while an
OPALX MONITOR is an unbounded plane, so the disc radius is applied to the OPALX records before
any count is compared. At LEMSdet that alone takes OPALX from 199 to 145 against G4beamline's
145.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from opalxruns.g4bl import track_on_plane

HERE = Path(__file__).resolve().parent

# A statistical criterion, not a worst-case one, because this is a 2000-particle beam with a
# 185 mrad angular spread. A handful of large-amplitude particles sample the nonlinear fringes
# of the maps, where a micron of difference amplifies over the remaining metres. Those are
# reported by id and counted rather than hidden inside a wider tolerance.
TOL = {
    "trans_p95_um": 100.0,  # 95th percentile transverse difference
    "trans_p99_um": 200.0,  # 99th percentile
    "outlier_um": 500.0,  # what counts as an outlier
    "outlier_frac": 0.005,  # at most 0.5 % of particles may exceed it
    # |p|, not the momentum vector. The G4beamline record sits 0.5 mm upstream of the plane and
    # only its position is drifted onto it, so inside a field region the vectors are taken at
    # slightly different points and differ by that half-millimetre of kick. |p| is conserved
    # exactly in a magnetostatic field, so it compares cleanly wherever the record was taken --
    # and it is the sharpest test that both codes see the same field.
    "absp_keV": 5.0,
}

# OPALX monitor, G4beamline detector, centreline z of the plane [mm], detector radius [mm].
PLANES = [
    ("E03_LEMSdet", "LEMSdet", 540.0, 112.0),
    ("E04_MiddleOfSolenoidDet", "MiddleOfSolenoidDet", 1285.0, 250.0),
    ("E07_EndOfSolenoidDet", "EndOfSolenoidDet", 1700.0, 250.0),
    ("E42_virtDet", "virtDet", 19250.0, 200.0),
]

MUON_MASS = 105.6583755  # MeV


MUON = 105.6583755


def opalx(name: str, radius: float) -> dict[int, np.ndarray]:
    """Monitor records, cut to the matching G4beamline detector's disc: x, y [mm], p [MeV/c]."""
    with h5py.File(HERE / f"{name}.h5", "r") as h:
        s = h["Step#0"]
        x, y = np.array(s["x"]) * 1000.0, np.array(s["y"]) * 1000.0
        px, py, pz = (np.array(s[k]) * MUON for k in ("px", "py", "pz"))
        ids = np.array(s["id"]).astype(int)
    keep = np.hypot(x, y) < radius
    return {
        int(i): np.array([x[k], y[k], px[k], py[k], pz[k]])
        for k, i in enumerate(ids)
        if keep[k]
    }


def g4bl(name: str, plane_z: float) -> dict[int, np.ndarray]:
    """#BLTrackFile rows in centreline coordinates, drifted onto the plane."""
    return track_on_plane(HERE / f"{name}.txt", plane_z)


def main() -> int:
    failures = []

    print("tracking: transverse difference for particles reaching the plane in both codes")
    print(
        f"{'plane':26s} {'both':>6s} {'p50 um':>9s} {'p95 um':>9s} {'p99 um':>9s} "
        f"{'max um':>9s} {'d|p| keV':>11s}  "
    )
    print("-" * 96)

    for opal_name, g4_name, plane_z, radius in PLANES:
        a, b = opalx(opal_name, radius), g4bl(g4_name, plane_z)
        shared = sorted(set(a) & set(b))
        if not shared:
            failures.append(f"{opal_name}: no particles in common")
            continue

        worst = np.array([np.abs(a[i][:2] - b[i][:2]).max() for i in shared]) * 1000.0  # um
        dmom = np.array(
            [abs(np.linalg.norm(a[i][2:]) - np.linalg.norm(b[i][2:])) for i in shared]
        ) * 1000.0  # keV
        p50, p95, p99, mx = (
            float(np.percentile(worst, 50)),
            float(np.percentile(worst, 95)),
            float(np.percentile(worst, 99)),
            float(worst.max()),
        )
        outliers = [shared[k] for k in np.where(worst > TOL["outlier_um"])[0]]
        frac = len(outliers) / len(shared)

        ok = (
            p95 <= TOL["trans_p95_um"]
            and p99 <= TOL["trans_p99_um"]
            and frac <= TOL["outlier_frac"]
            and dmom.max() <= TOL["absp_keV"]
        )
        if not ok:
            failures.append(
                f"{opal_name}: p95 {p95:.1f}, p99 {p99:.1f} um, "
                f"{len(outliers)} outlier(s), max d|p| {dmom.max():.2f} keV/c"
            )
        print(
            f"{opal_name:26s} {len(shared):6d} {p50:9.1f} {p95:9.1f} {p99:9.1f} {mx:9.1f} "
            f"{dmom.max():11.2f}  {'ok' if ok else 'FAIL'}"
        )
        if outliers:
            print(
                f"  ! {len(outliers)} of {len(shared)} ({100 * frac:.2f} %) exceed "
                f"{TOL['outlier_um']:.0f} um: ids {outliers}"
            )

    print()
    print("transmission: counts, after applying the G4beamline detector disc to OPALX")
    print(f"{'plane':26s} {'g4bl':>6s} {'OPALX':>7s} {'both':>6s} {'g4bl only':>10s}")
    print("-" * 66)
    for opal_name, g4_name, plane_z, radius in PLANES:
        a, b = opalx(opal_name, radius), g4bl(g4_name, plane_z)
        print(
            f"{opal_name:26s} {len(b):6d} {len(a):7d} {len(set(a) & set(b)):6d} "
            f"{len(set(b) - set(a)):10d}"
        )
    print()
    print("  The G4beamline-only column is dominated by the collimator model, not by tracking:")
    print("  an OPALX aperture kills every particle outside it, while a G4beamline jaw is a")
    print("  finite block a wide-angle particle can miss. See the module docstring.")

    print()
    if failures:
        print(f"FAIL ({len(failures)}):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(
        f"PASS: p95 <= {TOL['trans_p95_um']:.0f} um, p99 <= {TOL['trans_p99_um']:.0f} um, "
        f"at most {100 * TOL['outlier_frac']:.1f} % beyond {TOL['outlier_um']:.0f} um, "
        f"and every particle's |p| within {TOL['absp_keV']:.0f} keV/c"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
