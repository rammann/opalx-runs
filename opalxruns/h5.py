"""Readers for OPALX .h5 output."""

from __future__ import annotations

import numpy as np


def read_monitor(path) -> dict[str, np.ndarray]:
    """An OPALX MONITOR .h5, sorted by id.

    Particles cross the plane over several tracking steps, so the file holds
    one Step#N group per step that saw a crossing and all of them have to be
    concatenated. Already in m, beta*gamma and seconds.
    """
    import h5py
    keys = ["id", "x", "y", "z", "px", "py", "pz", "time"]
    acc = {k: [] for k in keys}
    with h5py.File(path, "r") as f:
        steps = sorted((k for k in f.keys() if k.startswith("Step#")),
                       key=lambda s: int(s.split("#")[1]))
        if not steps:
            raise ValueError(f"{path}: no Step groups -- the plane recorded nothing.")
        for s in steps:
            g = f[s]
            for k in keys:
                acc[k].append(np.asarray(g[k]))
    out = {k: np.concatenate(v) for k, v in acc.items()}
    order = np.argsort(out["id"])
    out = {k: v[order] for k, v in out.items()}
    out["id"] = out["id"].astype(int)
    return out
