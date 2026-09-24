"""Shared layer for the OPALX post-processing scripts in this directory.

Holds two things:

* :class:`Run` -- resolves the files an OPALX run directory contains. Every
  attribute is ``None`` (or an empty list) when the file is absent, so callers
  can skip a step instead of crashing on a run that has no monitors, no
  ``timing.dat``, etc.
* The file readers and frame math the scripts share: ``timing.dat``,
  ``_DesignPath.dat``, ``_ElementPositions.txt``, and the co-moving -> lab
  rotation.

The ``.stat`` and per-particle ``.h5`` readers live in
:mod:`opalxruns.opalx_diagnostics` and are re-exported here for convenience --
study scripts import them from there, so they must keep their home.

Assumed run layout::

    <run>/
        <base>.in  <base>.stat  <base>.h5     deck, SDDS stats, bunch dumps
        MON_*.h5                              optional monitor planes
        timing.dat
        data/<base>_DesignPath.dat            reference particle, per step
            /<base>_ElementPositions.txt      per-element lab points

Run these with a Python that has numpy/h5py/matplotlib/vtk (miniconda base),
not the system ``python3``.
"""

from __future__ import annotations

import re
from functools import cached_property
from pathlib import Path

import numpy as np

from opalxruns.opalx_diagnostics import (  # noqa: F401  (re-exported for the scripts)
    list_h5_datasets,
    list_h5_steps,
    load_h5_dataset,
    parse_opal_stat,
)

ZHAT = np.array([0.0, 0.0, 1.0])


# --------------------------------------------------------------------------- #
# run discovery                                                                #
# --------------------------------------------------------------------------- #
class Run:
    """The files of one OPALX run directory.

    Parameters
    ----------
    directory : Path | str
        The run directory, i.e. the one OPALX was started in.
    base : str, optional
        Run basename to use. Only needed when the directory holds the output of
        more than one deck and the automatic choice picks the wrong one.

    Everything is resolved lazily and cached. Missing files give ``None``
    (single files) or ``[]`` (lists) rather than raising.
    """

    def __init__(self, directory: Path | str, base: str | None = None) -> None:
        self.dir = Path(directory).resolve()
        if not self.dir.is_dir():
            raise NotADirectoryError(self.dir)
        self._base_override = base

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Run({self.dir.name!r}, base={self.base!r})"

    # -- naming ------------------------------------------------------------ #
    @cached_property
    def base(self) -> str | None:
        """Run basename, i.e. the deck filename without ``.in``.

        A run directory can accumulate the output of several decks (an earlier
        deck's ``.stat``/``.h5`` are not cleaned up), so picking the first
        ``*.stat`` alphabetically is not safe. Preference order:

        1. the basename of ``data/*_DesignPath.dat`` / ``*_ElementPositions.txt``
           -- OPALX rewrites those every run, so they name the most recent one;
        2. the directory name, when a matching ``.stat`` or ``.in`` exists;
        3. the newest ``.stat``, else the newest ``.in``.

        Pass ``base=`` to the constructor to override.
        """
        if self._base_override:
            return self._base_override

        data = self.data_dir or self.dir
        for suffix in ("_DesignPath.dat", "_ElementPositions.txt"):
            hits = [p for p in sorted(data.glob(f"*{suffix}")) if _useful(p)]
            if len(hits) == 1:
                return hits[0].name[: -len(suffix)]

        for suffix in (".stat", ".in"):
            if _useful(self.dir / f"{self.dir.name}{suffix}"):
                return self.dir.name

        for suffix in ("*.stat", "*.in"):
            hits = [p for p in self.dir.glob(suffix) if _useful(p)]
            hits.sort(key=lambda p: p.stat().st_mtime)
            if hits:
                return re.sub(r"_c\d+$", "", hits[-1].stem)
        return None

    @cached_property
    def deck(self) -> Path | None:
        """The input deck. Prefers ``<base>.in``, else the only ``*.in``."""
        if self.base:
            cand = self.dir / f"{self.base}.in"
            if cand.is_file():
                return cand
        decks = sorted(self.dir.glob("*.in"))
        return decks[0] if decks else None

    # -- outputs ----------------------------------------------------------- #
    @cached_property
    def stat_files(self) -> list[Path]:
        """This run's ``.stat`` files.

        A multi-container run writes ``<base>_cN.stat``, so there can be
        several. Files belonging to another deck in the same directory are left
        out -- see :attr:`other_stat_files`.
        """
        return [p for p in self._all_stat_files if self._is_ours(p.stem)]

    @cached_property
    def other_stat_files(self) -> list[Path]:
        """``.stat`` files in this directory that belong to another deck."""
        return [p for p in self._all_stat_files if not self._is_ours(p.stem)]

    @cached_property
    def _all_stat_files(self) -> list[Path]:
        return sorted(p for p in self.dir.glob("*.stat") if _useful(p))

    def _is_ours(self, stem: str) -> bool:
        """True if ``stem`` is ``<base>`` or a container variant ``<base>_cN``."""
        if not self.base:
            return True
        return re.fullmatch(rf"{re.escape(self.base)}(_c\d+)?", stem) is not None

    @cached_property
    def bunch_h5(self) -> Path | None:
        """The phase-space dump ``<base>.h5`` (co-moving frame)."""
        if not self.base:
            return None
        cand = self.dir / f"{self.base}.h5"
        return cand if _useful(cand) else None

    @cached_property
    def monitor_h5(self) -> list[Path]:
        """Monitor plane files: every ``*.h5`` that is not a bunch dump.

        OPALX upper-cases element names, so these are usually ``MON_*.h5``, but
        the name is the user's choice. What actually distinguishes them is the
        ``type`` root attribute LossDataSink writes (``b"spatial"``), which is
        checked when h5py can open the file.
        """
        out = []
        for p in sorted(self.dir.glob("*.h5")):
            if self._is_ours(p.stem):
                continue
            if _is_loss_file(p):
                out.append(p)
        return out

    @cached_property
    def timing(self) -> Path | None:
        """``timing.dat``, or ``None`` when it is absent or empty.

        OPALX leaves a zero-byte ``timing.dat`` behind on some runs, which is
        the same as having none.
        """
        cand = self.dir / "timing.dat"
        return cand if _useful(cand) else None

    @cached_property
    def data_dir(self) -> Path | None:
        cand = self.dir / "data"
        return cand if cand.is_dir() else None

    @cached_property
    def design_path(self) -> Path | None:
        """``data/<base>_DesignPath.dat`` -- the sampled reference orbit."""
        return self._in_data("_DesignPath.dat")

    @cached_property
    def element_positions(self) -> Path | None:
        """``data/<base>_ElementPositions.txt`` -- per-element lab points."""
        return self._in_data("_ElementPositions.txt")

    def _in_data(self, suffix: str) -> Path | None:
        for parent in (self.data_dir, self.dir):
            if parent is None:
                continue
            if self.base:
                cand = parent / f"{self.base}{suffix}"
                if _useful(cand):
                    return cand
            hits = [p for p in sorted(parent.glob(f"*{suffix}")) if _useful(p)]
            if hits:
                return hits[0]
        return None

    # -- output directories ------------------------------------------------ #
    @property
    def paraview_dir(self) -> Path:
        """``<run>/paraview`` -- created on access."""
        d = self.dir / "paraview"
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def plots_dir(self) -> Path:
        """``<run>/plots`` -- created on access."""
        d = self.dir / "plots"
        d.mkdir(parents=True, exist_ok=True)
        return d

    # -- reporting --------------------------------------------------------- #
    def inventory(self) -> list[tuple[str, str]]:
        """``(label, description)`` pairs of what was found, for printing."""
        def one(p):
            return p.name if p else "-"

        def many(ps):
            if not ps:
                return "-"
            return f"{len(ps)} file(s): {ps[0].name}" + (" ..." if len(ps) > 1 else "")

        rows = [
            ("base", self.base or "-"),
            ("deck", one(self.deck)),
            ("stat", many(self.stat_files)),
            ("bunch h5", one(self.bunch_h5)),
            ("monitors", many(self.monitor_h5)),
            ("timing", one(self.timing)),
            ("design path", one(self.design_path)),
            ("element pos", one(self.element_positions)),
        ]
        if self.other_stat_files:
            names = ", ".join(p.name for p in self.other_stat_files)
            rows.append(("ignored", f"{names}  (another deck; use --base to select)"))
        return rows


def _useful(path: Path) -> bool:
    """True if ``path`` is a file with something in it.

    OPALX can leave zero-byte outputs behind (an empty ``timing.dat`` shows up
    in several runs here); those are treated as absent so a step skips instead
    of failing on an empty parse.
    """
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _is_loss_file(path: Path) -> bool:
    """True if ``path`` is a LossDataSink dump (monitor / collimator plane)."""
    import h5py

    try:
        with h5py.File(path, "r") as f:
            t = f.attrs.get("type")
    except Exception:  # noqa: BLE001 - unreadable file is not a monitor
        return False
    if t is None:
        return False
    if isinstance(t, bytes):
        t = t.decode()
    return str(np.asarray(t).reshape(-1)[0]).strip("b'\"") == "spatial"


# --------------------------------------------------------------------------- #
# timing.dat                                                                   #
# --------------------------------------------------------------------------- #
_TIMER_ROW = re.compile(r"^(\S.*?)\.*\s{2,}(\d+)\s+(.*)$")


def parse_timing(path: Path | str) -> list[dict]:
    """Parse ``timing.dat`` into one row per timer.

    The file holds two tables: a ``Wall tot`` one and, below it, a
    ``Wall max / min / avg`` one. Both are read and merged by timer name.

    On a multi-rank run the first table often lists only ``mainTimer``, so most
    timers come back with ``wall_tot = None`` and only max/min/avg. That is left
    as-is rather than backfilled -- ``wall_max`` is not a total, and callers
    should decide which figure to use.

    Returns
    -------
    list of dict
        ``{"name", "ranks", "wall_tot", "wall_max", "wall_min", "wall_avg"}``,
        in the order the timers first appear. Missing values are ``None``.
    """
    rows: dict[str, dict] = {}
    order: list[str] = []

    for line in Path(path).read_text().splitlines():
        line = line.rstrip()
        if not line or line.startswith("=") or "ranks" in line:
            continue
        m = _TIMER_ROW.match(line)
        if not m:
            continue
        name = m.group(1).strip().rstrip(".")
        try:
            ranks = int(m.group(2))
            vals = [float(v) for v in m.group(3).split()]
        except ValueError:
            continue

        row = rows.setdefault(
            name,
            {"name": name, "ranks": ranks, "wall_tot": None,
             "wall_max": None, "wall_min": None, "wall_avg": None},
        )
        if name not in order:
            order.append(name)
        row["ranks"] = ranks
        if len(vals) == 1:
            row["wall_tot"] = vals[0]
        elif len(vals) >= 3:
            row["wall_max"], row["wall_min"], row["wall_avg"] = vals[:3]

    return [rows[n] for n in order]


# --------------------------------------------------------------------------- #
# data/<base>_DesignPath.dat                                                   #
# --------------------------------------------------------------------------- #
def load_designpath(path: Path | str):
    """Read the reference orbit: ``(s, R, P)``, sorted by path length.

    Columns of the file are ``s Rx Ry Rz Px Py Pz ...``. ``R`` is the lab
    position [m] and ``P`` the momentum (beta*gamma) of the reference particle,
    both ``(N, 3)``.
    """
    s, R, P = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 7:
                continue
            try:
                s.append(float(f[0]))
                R.append([float(f[1]), float(f[2]), float(f[3])])
                P.append([float(f[4]), float(f[5]), float(f[6])])
            except ValueError:
                continue
    s = np.asarray(s)
    R = np.asarray(R)
    P = np.asarray(P)
    order = np.argsort(s)
    return s[order], R[order], P[order]


# --------------------------------------------------------------------------- #
# data/<base>_ElementPositions.txt                                             #
# --------------------------------------------------------------------------- #
_ELEMPOS_ROW = re.compile(
    r'"([^:]+):\s*(\w+)"\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)'
)


def parse_element_positions(path: Path | str) -> list[tuple[str, np.ndarray]]:
    """Read per-element lab-frame centerlines.

    Each element contributes BEGIN / MID... / END rows (plus ENTRY/EXIT EDGE for
    bends), so a bent dipole comes out as a polyline that already follows its
    arc -- no placement or bend math is redone here.

    File column order is ``Z X Y``, so the points are reordered to ``(x, y, z)``.
    Consecutive duplicate points (ENTRY == BEGIN, END == EXIT) are dropped.

    Returns
    -------
    list of (name, points)
        ``points`` is ``(N, 3)`` in beam order.
    """
    order: list[str] = []
    pts: dict[str, list[tuple[float, float, float]]] = {}
    with open(path) as fh:
        for line in fh:
            m = _ELEMPOS_ROW.match(line)
            if not m:
                continue
            name = m.group(2)
            c1, c2, c3 = float(m.group(3)), float(m.group(4)), float(m.group(5))
            if name not in pts:
                pts[name] = []
                order.append(name)
            pts[name].append((c2, c3, c1))  # Z X Y -> x y z

    out = []
    for name in order:
        P = np.asarray(pts[name], dtype=float)
        keep = [0]
        for i in range(1, len(P)):
            if np.linalg.norm(P[i] - P[keep[-1]]) > 1e-9:
                keep.append(i)
        out.append((name, P[keep]))
    return out


# --------------------------------------------------------------------------- #
# co-moving -> lab frame                                                       #
# --------------------------------------------------------------------------- #
def shortest_arc(a, b) -> np.ndarray:
    """Rotation matrix taking unit vector ``a`` onto ``b`` with no roll."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = np.linalg.norm(v)
    if s < 1e-12:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def build_transport_frames(dhat: np.ndarray) -> np.ndarray:
    """Co-moving -> lab rotation ``Q(k)`` at every reference-orbit sample.

    OPALX realigns its frame every tracking step by the shortest-arc rotation
    putting local +z on the new reference momentum
    (``ParticleContainer::updateRefToLabCSTrafo``). Composing those increments is
    parallel transport of the momentum direction, which reproduces OPALX's own
    ``toLabTrafo`` -- including the accumulated frame roll that a single
    shortest-arc from +z to the final momentum would miss.

    ``Q(k) @ [0,0,1] == dhat[k]`` by construction.
    """
    dhat = dhat / np.linalg.norm(dhat, axis=1, keepdims=True)
    Q = np.empty((len(dhat), 3, 3))
    Q[0] = shortest_arc(ZHAT, dhat[0])
    for k in range(1, len(dhat)):
        Q[k] = shortest_arc(dhat[k - 1], dhat[k]) @ Q[k - 1]
    return Q


def h5_steps(path: Path | str) -> list[str]:
    """Sorted ``Step#N`` group names of an H5Part file."""
    return [f"Step#{n}" for n in list_h5_steps(path)]


__all__ = [
    "Run",
    "ZHAT",
    "build_transport_frames",
    "h5_steps",
    "list_h5_datasets",
    "list_h5_steps",
    "load_designpath",
    "load_h5_dataset",
    "parse_element_positions",
    "parse_opal_stat",
    "parse_timing",
    "shortest_arc",
]
