"""The muE4 reference muon, and the functions that assume it.

Both codes are given muons at P0 = 28 MeV/c: OPALX as the REAL P0, G4beamline as
referenceMomentum. read_bltrack() turns G4beamline momenta into beta*gamma with the
muon mass, and canonical() and matched() measure zeta and delta against this
reference, so none of them is right for another particle or momentum.

Used by the G4BL comparisons of single elements (cmplib.py) and of the whole muE4 line.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants. These must match the decks: P0 is the value OPALX gets as the REAL
# P0 and G4beamline as referenceMomentum.
# ---------------------------------------------------------------------------
MUON_MASS = 0.1056583755                    # GeV
P0 = 0.028                                  # GeV/c
BG0 = P0 / MUON_MASS                        # beta*gamma = 0.2650050
E0 = np.sqrt(P0**2 + MUON_MASS**2)          # GeV
BETA0 = P0 / E0                             # 0.2561551
BRHO = P0 / 0.299792458                     # 0.093398 T*m
C_MM_PER_NS = 299.792458


BLTRACK_COLS = ["x", "y", "z", "Px", "Py", "Pz", "t",
                "PDGid", "EventID", "TrackID", "ParentID", "Weight"]


def read_bltrack(path) -> dict[str, np.ndarray]:
    """A G4beamline #BLTrackFile, sorted by EventID.

    Returns metres, beta*gamma and seconds, so it lines up with read_monitor.
    G4beamline numbers events from 1 and OPALX ids from 0, so `id` is
    EventID - 1 and is what both sides match on. Matching on row order instead
    pairs different particles: OPALX writes monitor rows in whatever order the
    ranks produce them.
    """
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        rows.append([float(v) for v in parts[:12]])
    if not rows:
        raise ValueError(f"{path}: no particles. A plane outside the world, or past ZSTOP, "
                         f"records nothing and does not complain.")
    a = np.array(rows)
    d = {name: a[:, i] for i, name in enumerate(BLTRACK_COLS)}
    order = np.argsort(d["EventID"])
    out = {k: v[order] for k, v in d.items()}
    out["id"] = out["EventID"].astype(int) - 1
    out["x"] *= 1e-3
    out["y"] *= 1e-3
    out["z"] *= 1e-3
    for k in ("Px", "Py", "Pz"):
        out[k] = out[k] * 1e-3 / MUON_MASS   # MeV/c -> beta*gamma
    out["t"] *= 1e-9                          # ns -> s
    return out


def canonical(d: dict[str, np.ndarray], code: str) -> np.ndarray:
    """(x, x', y, y', zeta, delta) as a 6 x N array, from either code's reader.

    x' and y' are momentum ratios and carry no frame convention, which is what
    makes them the primary comparison. zeta is measured from particle 0 of the
    same file, so a constant time offset between the codes cannot leak into it;
    that offset is reported separately by time_of_flight().
    """
    if code == "g4bl":
        px, py, pz, t = d["Px"], d["Py"], d["Pz"], d["t"]
    else:
        px, py, pz, t = d["px"], d["py"], d["pz"], d["time"]
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    zeta = -BETA0 * 299792458.0 * (t - t[0])
    return np.vstack([d["x"], px / pz, d["y"], py / pz, zeta, pmag / BG0 - 1.0])


def matched(dg: dict, do: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
    """Line the two codes up on particle id.

    Returns (X_g4bl, X_opalx, ids, missing_from_opalx). Dropping to the common
    set quietly would bias every number that follows, so what was dropped is
    handed back and reported.

    Where a particle crossed the plane twice, the FIRST crossing is the one
    compared, in both codes -- searchsorted lands on the first of a run of equal
    ids. duplicates() says when that choice was made.
    """
    common = np.intersect1d(dg["id"], do["id"])
    missing = sorted(set(dg["id"].tolist()) - set(do["id"].tolist()))
    gi = np.searchsorted(dg["id"], common)
    oi = np.searchsorted(do["id"], common)
    Xg = canonical({k: v[gi] for k, v in dg.items()}, "g4bl")
    Xo = canonical({k: v[oi] for k, v in do.items()}, "opalx")
    return Xg, Xo, common, missing
