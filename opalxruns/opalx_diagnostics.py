"""Helpers for interactively browsing OPALX `.stat` and `.h5` outputs.

Filesystem convention: a root folder -- for example a topic's folder under the
repo's ``output/``, such as ``output/elements/aperture`` -- holds one folder per run::

    <root>/
        <run-name>/
            <run-name>[_cN].stat   # per-container SDDS-like statistics
            <run-name>[_cN].h5     # per-container H5Part phase-space dumps

This module provides:

* File-discovery helpers (:func:`discover_runs`, :func:`discover_files`)
* Parsers (:func:`parse_opal_stat`, :func:`list_h5_steps`,
  :func:`list_h5_datasets`, :func:`load_h5_dataset`)
* Cascading-dropdown widget factories (:class:`StatSelector`,
  :class:`H5Selector`) for use inside a notebook with
  ``ipywidgets.interact``.

The widget classes only touch :mod:`ipywidgets` inside their constructors,
so the rest of the module can be imported in a plain Python script
without Jupyter installed.
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# File discovery
# -----------------------------------------------------------------------------


def discover_runs(output_root: Path | str) -> list[Path]:
    """Return sorted run subdirectories of ``output_root``.

    Hidden directories (starting with ``.``) are skipped.
    """
    root = Path(output_root)
    if not root.is_dir():
        return []
    return sorted(
        p for p in root.iterdir() if p.is_dir() and not p.name.startswith(".")
    )


def discover_files(run_dir: Path | str, suffix: str) -> list[Path]:
    """Return sorted files inside ``run_dir`` ending in ``suffix``.

    ``suffix`` should include the leading dot, e.g. ``".stat"`` or ``".h5"``.
    Non-recursive (per-container files live at the top of the run dir).
    """
    d = Path(run_dir)
    if not d.is_dir():
        return []
    return sorted(p for p in d.iterdir() if p.is_file() and p.suffix == suffix)


# -----------------------------------------------------------------------------
# Parsers
# -----------------------------------------------------------------------------


def parse_opal_stat(path: Path | str) -> tuple[dict, pd.DataFrame]:
    """Parse an OPALX SDDS-like ``.stat`` file.

    Returns
    -------
    (metadata, df)
        ``metadata`` is a ``{name: value_string}`` mapping populated from the
        SDDS ``&parameter`` block in the order they appear; ``df`` is a
        :class:`pandas.DataFrame` indexed by row, with columns named after the
        ``&column`` declarations.
    """
    path = Path(path)
    lines = path.read_text().splitlines()

    col_names: list[str] = []
    param_names: list[str] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "&parameter":
            name = None
            i += 1
            while i < len(lines) and lines[i].strip() != "&end":
                m = re.search(r"name\s*=\s*([^,]+)", lines[i])
                if m:
                    name = m.group(1).strip()
                i += 1
            if name:
                param_names.append(name)
        elif line == "&column":
            name = None
            i += 1
            while i < len(lines) and lines[i].strip() != "&end":
                m = re.search(r"name\s*=\s*([^,]+)", lines[i])
                if m:
                    name = m.group(1).strip()
                i += 1
            if name:
                col_names.append(name)
        elif line == "&data":
            i += 1
            while i < len(lines) and lines[i].strip() != "&end":
                i += 1
            i += 1
            break
        i += 1

    body = [ln.strip() for ln in lines[i:] if ln.strip()]
    metadata_values = body[: len(param_names)]
    data_lines = body[len(param_names) :]

    metadata = dict(zip(param_names, metadata_values))

    rows: list[list[float]] = []
    for ln in data_lines:
        toks = ln.split()
        if len(toks) < len(col_names):
            continue
        rows.append(
            [
                np.nan if tok.lower() == "nan" else float(tok)
                for tok in toks[: len(col_names)]
            ]
        )

    return metadata, pd.DataFrame(rows, columns=col_names)


def list_h5_steps(path: Path | str) -> list[int]:
    """Return the sorted list of ``Step#N`` indices in an H5Part file."""
    import h5py  # local import: keeps top-level imports cheap

    with h5py.File(path, "r") as f:
        return sorted(int(name[5:]) for name in f.keys() if re.fullmatch(r"Step#\d+", name))


def list_h5_datasets(path: Path | str, step: int) -> list[str]:
    """Return sorted names of 1-D per-particle datasets in ``/Step#<step>``.

    Returns an empty list if the step group is absent.
    """
    import h5py

    grp_name = f"Step#{int(step)}"
    with h5py.File(path, "r") as f:
        if grp_name not in f:
            return []
        g = f[grp_name]
        return sorted(
            name
            for name in g.keys()
            if isinstance(g[name], h5py.Dataset) and g[name].ndim == 1
        )


def load_h5_dataset(path: Path | str, step: int, name: str) -> np.ndarray:
    """Load ``/Step#<step>/<name>`` from ``path`` into a NumPy array.

    Opens and closes the file each call — not optimised for repeated reads
    on the same step. For interactive notebook use that's irrelevant.
    """
    import h5py

    grp_name = f"Step#{int(step)}"
    with h5py.File(path, "r") as f:
        return np.asarray(f[grp_name][name])


# -----------------------------------------------------------------------------
# Widget factories
# -----------------------------------------------------------------------------


def _dropdown_options(paths: Iterable[Path], root: Path) -> list[tuple[str, str]]:
    """Make ``(label, value)`` pairs for an ipywidgets ``Dropdown``.

    ``label`` is the path relative to ``root`` when possible (shorter), and
    ``value`` is the absolute path string (what the callback consumes).
    """
    out: list[tuple[str, str]] = []
    for p in paths:
        try:
            label = str(p.relative_to(root))
        except ValueError:
            label = p.name
        out.append((label, str(p)))
    return out


class StatSelector:
    """Cascading ``run → .stat file → column`` dropdown set.

    Parameters
    ----------
    output_root : Path | str
        Directory containing run subdirectories. Resolved on construction.
    preferred_columns : Sequence[str], optional
        Column names tried in order for the initial column selection.
        Defaults favour the polarization columns.

    Attributes
    ----------
    widget : ipywidgets.VBox
        Composite layout suitable for :func:`IPython.display.display`.
    run, file, column : ipywidgets.Dropdown
        Individual dropdowns, exposed so the notebook can call
        ``widgets.interact(plot_fn, **selector.interact_kwargs)``.
    """

    def __init__(
        self,
        output_root: Path | str,
        preferred_columns: Sequence[str] = (
            "mean_polz",
            "mean_pol_mag",
            "mean_polx",
            "mean_poly",
            "numParticles",
            "energy",
        ),
    ) -> None:
        import ipywidgets as widgets

        self._root = Path(output_root).resolve()
        self._preferred_columns = list(preferred_columns)

        runs = discover_runs(self._root)
        self.run = widgets.Dropdown(
            options=_dropdown_options(runs, self._root),
            description="run:",
            layout=widgets.Layout(width="80%"),
        )
        self.file = widgets.Dropdown(
            options=[], description=".stat:", layout=widgets.Layout(width="80%")
        )
        self.column = widgets.Dropdown(
            options=[], description="column:", layout=widgets.Layout(width="60%")
        )

        self.run.observe(self._on_run_change, names="value")
        self.file.observe(self._on_file_change, names="value")

        # Initialise the dependent dropdowns.
        self._on_run_change({"new": self.run.value})

        self.widget = widgets.VBox([self.run, self.file, self.column])

    @property
    def interact_kwargs(self) -> dict:
        """Kwargs dict for :func:`ipywidgets.interact` plot callbacks."""
        return {"stat_path": self.file, "column": self.column}

    def attach_plot(self, plot_fn):
        """Bind ``plot_fn(stat_path, column)`` to this selector and display both.

        Uses :func:`ipywidgets.interactive_output` so the dropdowns are shown
        once (via :attr:`widget`) rather than duplicated by ``interact``.
        Returns the output widget so the caller can re-display or clear it.
        """
        import ipywidgets as widgets
        from IPython.display import display

        out = widgets.interactive_output(plot_fn, self.interact_kwargs)
        display(self.widget, out)
        return out

    # ------------------------------------------------------------------
    # Observers
    # ------------------------------------------------------------------
    def _on_run_change(self, change: dict) -> None:
        new_run = change.get("new")
        if not new_run:
            self.file.options = []
            return
        files = discover_files(Path(new_run), ".stat")
        self.file.options = _dropdown_options(files, self._root)
        # Triggers _on_file_change via the new value assignment.

    def _on_file_change(self, _change=None) -> None:
        path = self.file.value
        if not path:
            self.column.options = []
            return
        try:
            _, df = parse_opal_stat(path)
        except Exception:  # noqa: BLE001  (notebook context: surface but keep going)
            self.column.options = []
            return
        cols = [c for c in df.columns if c not in ("t", "s")]
        # Preserve user selection across file changes when the same column
        # exists in the new file; only rewrite options if the list differs.
        current = self.column.value
        if list(self.column.options) != cols:
            self.column.options = cols
        if current in cols:
            self.column.value = current
            return
        if cols:
            for name in self._preferred_columns:
                if name in cols:
                    self.column.value = name
                    return
            self.column.value = cols[0]


class H5Selector:
    """Cascading ``run → .h5 file → step → dataset`` selector.

    Parameters
    ----------
    output_root : Path | str
        Directory containing run subdirectories.
    preferred_datasets : Sequence[str], optional
        Dataset names tried in order for the initial dataset selection.
        Defaults favour polarization, then position.

    Attributes
    ----------
    widget : ipywidgets.VBox
        Composite layout for ``display(...)``.
    run, file, step, dataset
        Individual widgets, exposed for ``widgets.interact``.
    """

    def __init__(
        self,
        output_root: Path | str,
        preferred_datasets: Sequence[str] = (
            "polx",
            "poly",
            "polz",
            "x",
            "px",
        ),
    ) -> None:
        import ipywidgets as widgets

        self._root = Path(output_root).resolve()
        self._preferred_datasets = list(preferred_datasets)

        runs = discover_runs(self._root)
        self.run = widgets.Dropdown(
            options=_dropdown_options(runs, self._root),
            description="run:",
            layout=widgets.Layout(width="80%"),
        )
        self.file = widgets.Dropdown(
            options=[], description=".h5:", layout=widgets.Layout(width="80%")
        )
        self.step = widgets.IntSlider(
            description="step:", value=0, min=0, max=0, continuous_update=False
        )
        self.dataset = widgets.Dropdown(
            options=[], description="dataset:", layout=widgets.Layout(width="60%")
        )

        self.run.observe(self._on_run_change, names="value")
        self.file.observe(self._on_file_change, names="value")
        self.step.observe(self._on_step_change, names="value")

        self._on_run_change({"new": self.run.value})

        self.widget = widgets.VBox([self.run, self.file, self.step, self.dataset])

    @property
    def interact_kwargs(self) -> dict:
        return {"h5_path": self.file, "step": self.step, "dataset": self.dataset}

    def attach_plot(self, plot_fn):
        """Bind ``plot_fn(h5_path, step, dataset)`` and display the controls + plot."""
        import ipywidgets as widgets
        from IPython.display import display

        out = widgets.interactive_output(plot_fn, self.interact_kwargs)
        display(self.widget, out)
        return out

    # ------------------------------------------------------------------
    # Observers
    # ------------------------------------------------------------------
    def _on_run_change(self, change: dict) -> None:
        new_run = change.get("new")
        if not new_run:
            self.file.options = []
            return
        files = discover_files(Path(new_run), ".h5")
        self.file.options = _dropdown_options(files, self._root)

    def _on_file_change(self, _change=None) -> None:
        path = self.file.value
        if not path:
            self.step.min = self.step.max = self.step.value = 0
            self.dataset.options = []
            return
        try:
            steps = list_h5_steps(path)
        except Exception:  # noqa: BLE001
            steps = []
        if not steps:
            self.step.min = self.step.max = self.step.value = 0
            self.dataset.options = []
            return
        self.step.min = steps[0]
        self.step.max = steps[-1]
        self.step.step = max(1, steps[1] - steps[0]) if len(steps) > 1 else 1
        self.step.value = steps[-1]
        # The step change above fires _on_step_change which refreshes the
        # dataset dropdown. Force one in case the value was already equal.
        self._on_step_change(None)

    def _on_step_change(self, _change=None) -> None:
        path = self.file.value
        if not path:
            self.dataset.options = []
            return
        names = list_h5_datasets(path, self.step.value)
        # In H5Part the dataset set is identical across steps in any sane
        # file. Short-circuit when the new list matches the current options
        # so the dropdown's selected value is not reset on every step move.
        if list(self.dataset.options) == names:
            return
        current = self.dataset.value
        self.dataset.options = names
        if current in names:
            self.dataset.value = current
            return
        if names:
            for n in self._preferred_datasets:
                if n in names:
                    self.dataset.value = n
                    return
            self.dataset.value = names[0]


__all__ = [
    "discover_runs",
    "discover_files",
    "parse_opal_stat",
    "list_h5_steps",
    "list_h5_datasets",
    "load_h5_dataset",
    "StatSelector",
    "H5Selector",
]
