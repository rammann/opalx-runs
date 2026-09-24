"""One case folder of a scripted study and the OPALX output in it."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from opalxruns.opalx_diagnostics import list_h5_steps, parse_opal_stat
from opalxruns.paths import output_dir


class Case:
    """One case folder: its cases.json entry and the OPALX output of its run.

    Holds what the scripted studies share: the .stat table, the reference orbit,
    the bunch dumps in the .h5 file with their path lengths, one dump as the six
    coordinates, and the path of one particle. Each study subclasses it for its
    own checks.

    ``dir`` is the case folder under studies/ (the inputs); the output is read from
    ``out``, by default the case's folder under output/.
    """

    def __init__(self, folder: Path, entry: dict, bg0: float, out: Path | None = None):
        self.dir = Path(folder)
        self.name = self.dir.name
        self.m = entry
        self.out = Path(out) if out is not None else output_dir(self.dir)
        self.h5 = self.out / entry["h5"]
        self.stat = self.out / entry["stat"]
        self.bg0 = bg0                  # reference beta*gamma, for delta
        self._stat = None
        self._steps = None

    # -- stat file ---------------------------------------------------------
    @property
    def stat_df(self):
        if self._stat is None:
            _, self._stat = parse_opal_stat(self.stat)
        return self._stat

    def ref_orbit(self):
        df = self.stat_df
        return (df["s"].to_numpy(), df["ref_x"].to_numpy(), df["ref_z"].to_numpy(),
                df["ref_px"].to_numpy(), df["ref_pz"].to_numpy())

    # -- particle dumps ----------------------------------------------------
    def _steps_spos(self):
        """(dump numbers in order, path length of each dump)."""
        if self._steps is None:
            steps = sorted(list_h5_steps(self.h5))
            with h5py.File(self.h5, "r") as f:
                sp = np.array([float(np.ravel(f[f"Step#{n}"].attrs["SPOS"])[0])
                               for n in steps])
            self._steps = (steps, sp)
        return self._steps

    def read_plane(self, step: int, project: bool = True) -> np.ndarray:
        """6 x Npart phase space (x, x', y, y', z, delta) at one dump, particles
        sorted by id.

        Particles in one dump share a time rather than a path length, so by
        default each is projected back onto the reference transverse plane by
        drifting it through -z. That projection is exact in field-free space."""
        with h5py.File(self.h5, "r") as f:
            g = f[f"Step#{step}"]
            o = np.argsort(np.asarray(g["id"]))
            x = np.asarray(g["x"])[o]
            y = np.asarray(g["y"])[o]
            z = np.asarray(g["z"])[o]
            px = np.asarray(g["px"])[o]
            py = np.asarray(g["py"])[o]
            pz = np.asarray(g["pz"])[o]
        xp, yp = px / pz, py / pz
        delta = np.sqrt(px**2 + py**2 + pz**2) / self.bg0 - 1.0
        if project:
            x = x - xp * z
            y = y - yp * z
        return np.vstack([x, xp, y, yp, z, delta])

    def _pmag(self, step):
        """|p| of every particle at one dump, sorted by id, in beta*gamma."""
        with h5py.File(self.h5, "r") as f:
            g = f[f"Step#{step}"]
            o = np.argsort(np.asarray(g["id"]))
            return np.sqrt(np.asarray(g["px"])[o]**2 + np.asarray(g["py"])[o]**2
                           + np.asarray(g["pz"])[o]**2)

    def trajectory(self, part_index: int = 0):
        """(s, x, y) of one particle over every dump."""
        steps, sp = self._steps_spos()
        xs, ys = [], []
        with h5py.File(self.h5, "r") as f:
            for n in steps:
                g = f[f"Step#{n}"]
                o = np.argsort(np.asarray(g["id"]))
                xs.append(float(np.asarray(g["x"])[o][part_index]))
                ys.append(float(np.asarray(g["y"])[o][part_index]))
        return sp, np.array(xs), np.array(ys)
