#!/usr/bin/env python
"""compare.py -- OPALX against G4beamline, particle by particle.

Both codes report at the same detectors, in the same frame: G4beamline's virtualdetectors are
written with `coordinates=centerline`, and an OPALX MONITOR records in its own local frame,
which for a detector placed on the centreline with no rotation is the same frame. So the
transverse coordinates compare directly, with no lifting.

Two things have to be got right or the answer is meaningless.

**Use centreline coordinates, not global.** G4beamline writes 6 significant figures. In global
coordinates the last detector sits at z = 16549 mm, so the printed resolution there is 0.1 mm
-- which, projected into the tilted detector plane, is 36 um of pure rounding noise. That is
larger than the difference being measured, so a global-frame comparison cannot resolve it at
all. In centreline coordinates the transverse numbers are millimetres and the floor is ~1e-5
mm.

**Drift the G4beamline record onto the plane.** A virtualdetector is a 1 mm thick volume and
records on entry, so its z is half a millimetre upstream of the plane an OPALX monitor
interpolates to. At a 40 mrad angle that is 20 um -- again comparable to the difference. Each
G4beamline record is drifted forward along its own momentum before comparing.

Particles are matched on identity, never on row order: OPALX writes monitor rows in rank
order, and G4beamline's EventID is the OPALX id plus one.

Exit code is non-zero if any tolerance is exceeded.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
POSES = json.loads((HERE.parent / "poses.json").read_text())

TOL = {
    "trans_um": 50.0,  # per-particle transverse position at the detector plane
    "mom_keV": 5.0,  # per-particle momentum
}

# OPALX monitor <-> G4beamline virtualdetector, and the centreline z of the plane [mm].
PLANES = [
    ("E03_LEMSdet", "LEMSdet", 540.0),
    ("E04_MiddleOfSolenoidDet", "MiddleOfSolenoidDet", 1285.0),
    ("E07_EndOfSolenoidDet", "EndOfSolenoidDet", 1700.0),
    ("E42_virtDet", "virtDet", 19250.0),
]

MUON_MASS = 105.6583755  # MeV


def opalx(name: str) -> dict[int, np.ndarray]:
    """Monitor records: x, y [mm] and px, py, pz [MeV/c], in the monitor's own frame."""
    with h5py.File(HERE / f"{name}.h5", "r") as h:
        s = h["Step#0"]
        rows = np.stack(
            [
                np.array(s["x"]) * 1000.0,
                np.array(s["y"]) * 1000.0,
                np.array(s["px"]) * MUON_MASS,
                np.array(s["py"]) * MUON_MASS,
                np.array(s["pz"]) * MUON_MASS,
            ]
        ).T
        ids = np.array(s["id"]).astype(int)
    return {int(i): rows[k] for k, i in enumerate(ids)}


def g4bl(name: str, plane_z: float) -> dict[int, np.ndarray]:
    """#BLTrackFile rows in centreline coordinates, drifted onto the plane."""
    out = {}
    for line in (HERE / f"{name}.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        x, y, z, px, py, pz = (float(v) for v in c[:6])
        dz = plane_z - z  # the detector is 1 mm thick and records on entry
        out[int(c[8]) - 1] = np.array([x + px / pz * dz, y + py / pz * dz, px, py, pz])
    return out


def main() -> int:
    failures = []
    print(
        f"{'plane':26s} {'n':>4s} {'max |dx| um':>12s} {'max |dy| um':>12s} {'max |dp| keV/c':>15s}  "
    )
    print("-" * 82)

    for opal_name, g4_name, plane_z in PLANES:
        try:
            a, b = opalx(opal_name), g4bl(g4_name, plane_z)
        except (FileNotFoundError, OSError) as exc:
            failures.append(f"{opal_name}: {exc}")
            continue

        shared = sorted(set(a) & set(b))
        d = np.array([a[i] - b[i] for i in shared])
        max_dx = float(np.abs(d[:, 0]).max()) * 1000.0
        max_dy = float(np.abs(d[:, 1]).max()) * 1000.0
        max_dp = float(np.abs(np.linalg.norm(d[:, 2:], axis=1)).max()) * 1000.0

        ok = max(max_dx, max_dy) <= TOL["trans_um"] and max_dp <= TOL["mom_keV"]
        if not ok:
            failures.append(
                f"{opal_name}: dx {max_dx:.1f} um, dy {max_dy:.1f} um, dp {max_dp:.2f} keV/c"
            )
        print(
            f"{opal_name:26s} {len(shared):4d} {max_dx:12.1f} {max_dy:12.1f} {max_dp:15.2f}  "
            f"{'ok' if ok else 'FAIL'}"
        )

        # OPALX monitors are known to drop particles. Name them and say what distinguishes
        # them, rather than letting the matched set quietly shrink.
        dropped = sorted(set(b) - set(a))
        if dropped:
            mom = {i: float(np.linalg.norm(b[i][2:])) for i in b}
            lowest = min(mom, key=mom.get)
            note = " -- the lowest-momentum particle" if dropped == [lowest] else ""
            detail = ", ".join(f"id {i} at {mom[i]:.3f} MeV/c" for i in dropped)
            print(f"  ! {opal_name} recorded {len(a)} of {len(b)}: missing {detail}{note}")

    print()
    if failures:
        print(f"FAIL ({len(failures)}):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(
        f"PASS: transverse <= {TOL['trans_um']:.0f} um and momentum <= {TOL['mom_keV']:.0f} keV/c "
        f"per particle, at every plane"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
