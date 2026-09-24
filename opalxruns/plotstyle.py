"""Colours and matplotlib settings shared by the study plots."""

from __future__ import annotations

from pathlib import Path

# Two series, so two categorical colours, assigned to the code and never cycled.
# Checked with the palette validator: adjacent-pair separation dE 25.7 under
# protanopia and 32.5 in normal vision, both well above the floor.
OPALX = "#1f5fd1"       # OPALX data
G4BL = "#c1432d"        # G4beamline data
INK = "#222222"         # text
MUTED = "#777777"       # axes and ticks
GRID = "#dddddd"        # grid lines

STYLE = {
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": MUTED, "text.color": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.axisbelow": True,
    "legend.frameon": False, "figure.facecolor": "white",
}


def use(changes: dict | None = None) -> None:
    """Apply STYLE to matplotlib, with ``changes`` on top."""
    import matplotlib.pyplot as plt

    plt.rcParams.update({**STYLE, **(changes or {})})


def save_with_footer(fig, path, footer: str, dpi: int = 140) -> None:
    """Lay the figure out with room for one line of small monospace text at the
    bottom, write it to ``path`` and close it."""
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 1.0))
    fig.text(0.01, 0.012, footer, fontsize=6, family="monospace", color="#555555")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
