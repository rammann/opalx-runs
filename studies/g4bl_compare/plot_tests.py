#!/usr/bin/env python3
"""Figures for the g4bl_compare study.

The subject is the difference between the two codes, so most of these plot a
difference rather than a value. Only the first is context.

Called from run_tests.py, inside a try/except, so that a broken figure can
never hide a failing test. Can also be run on its own:
    python plot_tests.py [case ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.colors import LogNorm    # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmplib                            # noqa: E402
from cmplib import TOL                   # noqa: E402

PLOTS = HERE / "plots"

# Two series, so two categorical colours, assigned to the code and never cycled.
# Checked with the palette validator: adjacent-pair separation dE 25.7 under
# protanopia and 32.5 in normal vision, both well above the floor.
OPALX = "#1f5fd1"
G4BL = "#c1432d"
INK, MUTED, GRIDC = "#222222", "#777777", "#dddddd"
# Magnitude of a difference is a one-directional quantity, so one hue, light to
# dark. Signed differences use a two-hue scale with a neutral middle.
SEQ, DIV = "Blues", "RdBu_r"
# One fixed hue per case, assigned in order and never cycled. Checked with the
# palette validator: worst adjacent separation dE 9.1 under protanopia and 19.6
# in normal vision, both above the floor. Three of the seven sit below 3:1
# contrast on white, so every series also carries its own marker shape and the
# same numbers are printed in results.txt.
CASE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
               "#008300", "#4a3aa7"]
CASE_MARKERS = ["o", "s", "^", "D", "v", "<", ">"]

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "legend.frameon": False, "figure.facecolor": "white",
})


def _save(fig, name, footer):
    PLOTS.mkdir(exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 1.0))
    fig.text(0.01, 0.012, footer, fontsize=6, family="monospace", color="#555555")
    fig.savefig(PLOTS / f"{name}.png", dpi=140)
    plt.close(fig)
    print(f"  plots/{name}.png")


def _grid_axes(n, ncol=3, w=3.6, h=2.6, squeeze=False):
    nrow = int(np.ceil(n / ncol))
    fig, ax = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    return fig, ax.ravel(), nrow * ncol


def _field(d: Path, kind: str, boxes=None, drop_faces=True):
    """Both codes' field at the same points, by default with the box-face points
    removed -- the same points run_tests.py drops. On a face the two codes round
    the last bit of the coordinate differently, which shows up as a difference
    the size of the whole field and would otherwise dominate every figure."""
    fg = d / "pair" / f"g4bl_field_{kind}.txt"
    fo = d / "pair" / f"opalx_field_{kind}.dat"
    if not (fg.exists() and fo.exists()):
        return None
    pos, Bg, Bo = cmplib.match_field(cmplib.read_g4bl_field(fg),
                                     cmplib.read_opalx_field(fo))
    if drop_faces and boxes:
        keep = ~cmplib.on_face(pos, boxes)
        if keep.any():
            pos, Bg, Bo = pos[:, keep], Bg[:, keep], Bo[:, keep]
    return pos, Bg, Bo


def _field_tol(Bg) -> float:
    """This case's field tolerance: three printing steps of its own peak field."""
    return TOL["field_between"] * cmplib.printing_step(float(np.abs(Bg).max()))


def _raw_g4bl(d: Path, z):
    """The G4beamline plane as read, for the floor calculation, which needs the
    printed quantities rather than the six coordinates."""
    return cmplib.read_bltrack(d / "pair" / f"Z{z}.txt")


def _planes(d: Path, c: dict, sub: str):
    """Every plane of one stage, as {z: (Xg, Xo, ids)}."""
    out = {}
    for z in c["monitors"]:
        try:
            dg = cmplib.read_bltrack(d / sub / f"Z{z}.txt")
            do = cmplib.read_monitor(d / sub / f"MON_{z}.h5")
            Xg, Xo, ids, _ = cmplib.matched(dg, do)
            out[z] = (Xg, Xo, ids)
        except Exception:
            continue
    return out


# ---------------------------------------------------------------------------
def plot_field_on_axis(man, names):
    """Context, not a result: what field each case actually has.

    Value on top, difference below, sharing the axis. Two panels rather than two
    y scales on one -- a second scale on the same frame invites reading one
    curve against the other's axis.
    """
    fig, axes, slots = _grid_axes(len(names) * 2, ncol=2, w=4.6, h=2.2)
    for k, n in enumerate(names):
        c, d = man[n], HERE / man[n]["dir"]
        m = _field(d, "grid", c["boxes"])
        a, b = axes[2 * k], axes[2 * k + 1]
        if m is None:
            a.set_title(f"{n} -- not run")
            a.axis("off"); b.axis("off")
            continue
        pos, Bg, Bo = m
        # the sampled line closest to the map's own axis
        xs = np.unique(pos[0])
        sel = np.isclose(pos[0], xs[np.argmin(np.abs(xs))])
        o = np.argsort(pos[2][sel])
        z = pos[2][sel][o] * 1000.0
        comp = int(np.argmax(np.abs(Bg[:, sel]).max(axis=1)))
        lbl = "Bx By Bz".split()[comp]
        a.plot(z, Bg[comp, sel][o], color=G4BL, lw=1.6, label="G4beamline")
        a.plot(z, Bo[comp, sel][o], color=OPALX, lw=1.0, ls="--", label="OPALX")
        a.set_title(f"{n}   {lbl} at x = {xs[np.argmin(np.abs(xs))] * 1000:.0f} mm, y = 0")
        a.set_ylabel(f"{lbl} [T]")
        a.legend(loc="best")
        b.plot(z, (Bo[comp, sel] - Bg[comp, sel])[o], color=INK, lw=0.9)
        b.axhline(0, color=MUTED, lw=0.5)
        b.set_ylabel("OPALX - G4BL [T]")
        b.set_xlabel("lab z [mm]")
    for j in range(len(names) * 2, slots):
        axes[j].axis("off")
    _save(fig, "01_field_on_axis",
          "the line of the y=0 sample plane nearest each map's axis; "
          "dominant component only")


def plot_field_diff_plane(man, names):
    """Where in the map the two codes disagree, if anywhere.

    One hue, light to dark: this is a magnitude, so it has no sign to show. The
    colour scale is shared down each column so cases can be read against each
    other, and floored at the level G4beamline's six printed figures can produce.
    """
    fig, axes, slots = _grid_axes(len(names) * 2, ncol=2, w=4.6, h=2.4)
    for k, n in enumerate(names):
        d = HERE / man[n]["dir"]
        for j, kind in enumerate(("grid", "half")):
            ax = axes[2 * k + j]
            m = _field(d, kind, man[n]["boxes"])
            if m is None:
                ax.axis("off")
                continue
            pos, Bg, Bo = m
            diff = np.abs(Bo - Bg).max(axis=0)
            xs, zs = np.unique(pos[0]), np.unique(pos[2])
            img = np.full((len(xs), len(zs)), np.nan)
            ix = np.searchsorted(xs, pos[0])
            iz = np.searchsorted(zs, pos[2])
            img[ix, iz] = diff
            floor = 1e-8
            im = ax.imshow(np.maximum(img, floor), origin="lower", aspect="auto",
                           cmap=SEQ, norm=LogNorm(vmin=floor, vmax=max(diff.max(), 1e-7)),
                           extent=[zs[0] * 1000, zs[-1] * 1000, xs[0] * 1000, xs[-1] * 1000])
            ax.set_title(f"{n}, {'at the map grid points' if kind == 'grid' else 'halfway between'}"
                         f"\nworst {diff.max():.1e} T")
            ax.set_xlabel("lab z [mm]")
            ax.set_ylabel("lab x [mm]")
            cb = fig.colorbar(im, ax=ax)
            cb.set_label("|OPALX - G4BL| [T]", fontsize=7)
    for j in range(len(names) * 2, slots):
        axes[j].axis("off")
    _save(fig, "02_field_diff_plane",
          "max over Bx,By,Bz of |OPALX - G4BL| on the y=0 sample plane, log scale")


def plot_field_diff_hist(man, names):
    """Whether a difference is in the table or in the interpolation.

    Sampling at the map's own grid points reads the table back; sampling halfway
    between makes both codes interpolate. If the two histograms sit on top of
    each other, the interpolation schemes agree as well as the tables do.
    """
    fig, axes, slots = _grid_axes(len(names), ncol=3)
    bins = np.logspace(-12, -4, 40)
    for k, n in enumerate(names):
        ax, d = axes[k], HERE / man[n]["dir"]
        drew, tol = False, None
        for kind, col, ls in (("grid", G4BL, "-"), ("half", OPALX, "--")):
            m = _field(d, kind, man[n]["boxes"])
            if m is None:
                continue
            _, Bg, Bo = m
            tol = _field_tol(Bg)
            v = np.abs(Bo - Bg).max(axis=0)
            ax.hist(np.clip(v, bins[0], bins[-1]), bins=bins, histtype="step",
                    color=col, lw=1.4, ls=ls,
                    label=f"{'at grid points' if kind == 'grid' else 'halfway between'}")
            drew = True
        if not drew:
            ax.axis("off")
            continue
        if tol:
            ax.axvline(tol, color=INK, lw=1.0, ls=":")
            ax.text(tol, ax.get_ylim()[1] * 0.6, f" tolerance {tol:.0e} T",
                    fontsize=6, color=INK, rotation=90, va="top")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(n)
        ax.set_xlabel("|OPALX - G4BL| [T]")
        ax.set_ylabel("sample points")
        ax.legend(loc="upper left", fontsize=6)
    for j in range(len(names), slots):
        axes[j].axis("off")
    _save(fig, "03_field_diff_hist",
          "sample points of the y=0 plane clear of a box face; the tolerance is "
          "three printing steps of each case's own peak field")


def plot_exit_diff(man, names):
    """Per-particle difference at the last plane, for the 19-particle set.

    The last six bars of each panel are the large-amplitude particles. If they
    stand out against the first twelve, the map stops agreeing away from the axis.
    """
    fig, axes, slots = _grid_axes(len(names), ncol=3, h=2.8)
    for k, n in enumerate(names):
        ax, c = axes[k], man[n]
        p = _planes(HERE / c["dir"], c, "pair")
        if not p:
            ax.axis("off")
            continue
        z = max(p)
        Xg, Xo, ids = p[z]
        dpos = np.abs(Xo[[0, 2]] - Xg[[0, 2]]).max(axis=0)
        dang = np.abs(Xo[[1, 3]] - Xg[[1, 3]]).max(axis=0)
        w = 0.4
        ax.bar(ids - w / 2, dpos, w, color=OPALX, label="|dx|,|dy| [m]")
        ax.bar(ids + w / 2, dang, w, color=G4BL, label="|dx'|,|dy'| [rad]")
        ax.axhline(TOL["exit_pos"], color=OPALX, lw=0.9, ls=":")
        ax.axhline(TOL["exit_angle"], color=G4BL, lw=0.9, ls=":")
        if len(ids) > 13:
            ax.axvline(12.5, color=MUTED, lw=0.8)
            ax.text(12.6, ax.get_ylim()[1], " large steps", fontsize=6,
                    color=MUTED, va="top")
        ax.set_yscale("log")
        ax.set_title(f"{n}, plane z = {z} mm")
        ax.set_xlabel("particle id")
        ax.legend(loc="lower right", fontsize=6)
    for j in range(len(names), slots):
        axes[j].axis("off")
    _save(fig, "04_exit_diff_bars",
          "19-particle stage; dotted lines are the tolerances; ids 13-18 are the "
          "20 mm and 2% steps")


def plot_matrix_diff(man, names):
    """The transfer matrix difference, absolute and relative.

    Signed on the left, so a two-hue scale with a neutral middle. Relative on the
    right, but only where the element is large enough for a ratio to mean
    anything -- the rest are zero by symmetry and are left blank rather than
    shown as enormous.
    """
    fig, axes, slots = _grid_axes(len(names) * 2, ncol=4, w=3.0, h=2.6)
    lab = cmplib.COORD_NAMES
    for k, n in enumerate(names):
        c = man[n]
        p = _planes(HERE / c["dir"], c, "pair")
        a, b = axes[2 * k], axes[2 * k + 1]
        if len(p) < 2 or len(p[min(p)][2]) != c["n_pair"]:
            a.axis("off"); b.axis("off")
            continue
        Xg_in, Xo_in, _ = p[min(p)]
        Xg_out, Xo_out, _ = p[max(p)]
        Mg = cmplib.transfer_matrix(Xg_in, Xg_out)
        Mo = cmplib.transfer_matrix(Xo_in, Xo_out)
        diff = Mo - Mg
        v = np.abs(diff).max() or 1e-12
        im = a.imshow(diff, cmap=DIV, vmin=-v, vmax=v)
        a.set_title(f"{n}\nOPALX - G4BL, worst {v:.1e}")
        fig.colorbar(im, ax=a)
        floor = float(cmplib.row_floors(_raw_g4bl(HERE / c['dir'], max(p)),
                                        cmplib.STEP_SMALL).max())
        sig = cmplib.significant(Mg, floor, TOL["matrix_big"])
        rel = np.where(sig, np.abs(diff) / np.where(sig, np.abs(Mg), 1), np.nan)
        im = b.imshow(rel, cmap=SEQ, vmin=0, vmax=max(np.nanmax(rel), 1e-12))
        b.set_title(f"relative, elements above 100x the floor\nworst "
                    f"{np.nanmax(rel) if sig.any() else float('nan'):.1e}, "
                    f"tol {TOL['matrix_rel']:.0e}")
        fig.colorbar(im, ax=b)
        for ax in (a, b):
            ax.set_xticks(range(6), lab, fontsize=6)
            ax.set_yticks(range(6), lab, fontsize=6)
            ax.grid(False)
    for j in range(len(names) * 2, slots):
        axes[j].axis("off")
    _save(fig, "05_matrix_diff",
          "6x6 by centred differences over the 1e-3 steps; blank cells are zero "
          "by symmetry")


def plot_orbit_vs_s(man, names):
    """Where along the case the difference appears.

    One point per recording plane, so this is coarse, but it is the figure that
    says whether a disagreement builds up through the map or arrives all at once
    at one place in it.
    """
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for k, n in enumerate(names):
        c = man[n]
        p = _planes(HERE / c["dir"], c, "pair")
        if not p:
            continue
        zs = sorted(p)
        worst = [np.abs(p[z][1][[0, 2]] - p[z][0][[0, 2]]).max() for z in zs]
        ax.plot(zs, np.maximum(worst, 1e-18), marker=CASE_MARKERS[k % 7], ms=5,
                lw=1.4, color=CASE_COLORS[k % 7], label=n)
    ax.axhline(TOL["exit_pos"], color=INK, lw=1.0, ls=":")
    ax.text(ax.get_xlim()[1], TOL["exit_pos"], " tolerance", fontsize=7,
            color=INK, va="bottom", ha="right")
    ax.set_yscale("log")
    ax.set_xlabel("lab z of the recording plane [mm]")
    ax.set_ylabel("worst |dx|,|dy| over the 19 particles [m]")
    ax.set_title("Where the two codes start to differ")
    ax.legend(loc="lower right", ncol=2, fontsize=7)
    _save(fig, "06_orbit_diff_vs_s",
          "19-particle stage, every recording plane of every case")


def plot_gauss_per_particle(man, names):
    """The statistical figure: does the difference depend on amplitude?

    20000 particles, matched by id, so this is still a per-particle difference
    and not a comparison of distribution widths. A map that is read correctly
    near the axis and wrongly further out shows up here as a rising trend and
    nowhere else.
    """
    fig, axes, slots = _grid_axes(len(names), ncol=3, h=2.8)
    for k, n in enumerate(names):
        ax, c = axes[k], man[n]
        p = _planes(HERE / c["dir"], c, "gauss")
        if not p:
            ax.set_title(f"{n} -- gauss stage not run")
            ax.axis("off")
            continue
        z0, z1 = min(p), max(p)
        Xg0, _, ids0 = p[z0]
        Xg1, Xo1, ids1 = p[z1]
        common = np.intersect1d(ids0, ids1)
        i0 = np.searchsorted(ids0, common)
        i1 = np.searchsorted(ids1, common)
        r = np.hypot(Xg0[0][i0], Xg0[2][i0]) * 1000.0          # start amplitude
        dv = np.abs(Xo1[[0, 2]][:, i1] - Xg1[[0, 2]][:, i1]).max(axis=0)
        ax.scatter(r, np.maximum(dv, 1e-12), s=3, alpha=0.18, color=OPALX,
                   edgecolors="none")
        # a median in amplitude bins, so the trend is readable through the cloud
        b = np.linspace(0, r.max(), 21)
        idx = np.clip(np.digitize(r, b) - 1, 0, len(b) - 2)
        med = np.array([np.median(dv[idx == i]) if (idx == i).any() else np.nan
                        for i in range(len(b) - 1)])
        ax.plot(0.5 * (b[1:] + b[:-1]), med, color=G4BL, lw=1.8, label="median")
        ax.axhline(TOL["exit_pos"], color=INK, lw=0.9, ls=":")
        ax.set_yscale("log")
        ax.set_title(f"{n}\nz = {z1} mm, {len(common)} particles")
        ax.set_xlabel("starting radius [mm]")
        ax.set_ylabel("|dx|,|dy| [m]")
        # bottom right: the tolerance line and the dense core are both at the top
        ax.legend(loc="lower right", fontsize=6)
    for j in range(len(names), slots):
        axes[j].axis("off")
    _save(fig, "07_gauss_per_particle",
          "20000-particle stage, last plane; dotted line is the tolerance. Flat "
          "with radius means the map is read as well off axis as on it; the "
          "quadrupole rises because its field does.")


def plot_gauss_moments(man, names):
    """Whether the bunch as a whole comes out the same size in the same place."""
    fig, axes, slots = _grid_axes(len(names), ncol=3, h=2.8)
    for k, n in enumerate(names):
        ax, c = axes[k], man[n]
        p = _planes(HERE / c["dir"], c, "gauss")
        if not p:
            ax.set_title(f"{n} -- gauss stage not run")
            ax.axis("off")
            continue
        zs = sorted(p)
        for j, lbl, ls in ((0, "x", "-"), (2, "y", "--")):
            ax.plot(zs, [p[z][0][j].std() * 1000 for z in zs], color=G4BL, ls=ls,
                    lw=1.6, label=f"rms {lbl}, G4beamline")
            ax.plot(zs, [p[z][1][j].std() * 1000 for z in zs], color=OPALX, ls=ls,
                    lw=1.0, marker="o", ms=3, label=f"rms {lbl}, OPALX")
        ax.set_title(n)
        ax.set_xlabel("lab z of the plane [mm]")
        ax.set_ylabel("rms [mm]")
        ax.legend(fontsize=5.5, ncol=2)
    for j in range(len(names), slots):
        axes[j].axis("off")
    _save(fig, "08_gauss_moments", "20000-particle stage, every recording plane")


def plot_summary(man, names):
    """One bar per case per quantity, as a fraction of that quantity's tolerance.

    Plotted as a ratio, not as a raw value, for two reasons. The three
    quantities are a field, a length and an angle, and stacking Tesla against
    metres against radians on one axis compares nothing. And the field tolerance
    is not the same number for every case -- it is three printing steps of that
    case's own peak field, which is 3e-7 T for the dipoles and 3e-6 T for the
    quadrupole. As a fraction of tolerance there is one line, at 1.
    """
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    rows = {"field": [], "exit position": [], "exit angle": []}
    for n in names:
        c, d = man[n], HERE / man[n]["dir"]
        m = _field(d, "half", c["boxes"])
        if m is not None:
            _, Bg, Bo = m
            tol = _field_tol(Bg)
            rows["field"].append(float(np.abs(Bo - Bg).max()) / tol if tol else np.nan)
        else:
            rows["field"].append(np.nan)
        p = _planes(d, c, "pair")
        if p:
            Xg, Xo, _ = p[max(p)]
            rows["exit position"].append(
                float(np.abs(Xo[[0, 2]] - Xg[[0, 2]]).max()) / TOL["exit_pos"])
            rows["exit angle"].append(
                float(np.abs(Xo[[1, 3]] - Xg[[1, 3]]).max()) / TOL["exit_angle"])
        else:
            rows["exit position"].append(np.nan)
            rows["exit angle"].append(np.nan)

    i = np.arange(len(names))
    w = 0.27
    for k, (lbl, vals) in enumerate(rows.items()):
        v = np.array(vals, dtype=float)
        ax.bar(i + (k - 1) * w, np.nan_to_num(v, nan=1e-6), w,
               color=CASE_COLORS[k], label=lbl)
        for xi, vi in zip(i + (k - 1) * w, v):
            if np.isfinite(vi):
                # Three of the palette slots are below 3:1 on white, so every bar
                # carries its number rather than relying on the colour alone.
                ax.text(xi, vi * 1.15, f"{vi:.2g}", ha="center", fontsize=5.5,
                        rotation=90, color=INK)
    ax.axhline(1.0, color=INK, lw=1.2, ls=":")
    ax.text(len(names) - 0.5, 1.0, " tolerance", fontsize=7, color=INK, va="bottom",
            ha="right")
    ax.set_yscale("log")
    ax.set_ylim(top=3.0)
    ax.set_xticks(i, names, rotation=20, ha="right", fontsize=7)
    ax.set_ylabel("worst difference / tolerance")
    ax.set_title("Every case against its own tolerance")
    ax.legend(ncol=3, fontsize=7, loc="lower left")
    _save(fig, "09_summary",
          "field from the halfway-between samples, off the box faces; tracking "
          "from the last plane of the 19-particle stage")


FIGURES = [plot_field_on_axis, plot_field_diff_plane, plot_field_diff_hist,
           plot_exit_diff, plot_matrix_diff, plot_orbit_vs_s,
           plot_gauss_per_particle, plot_gauss_moments, plot_summary]


def make_all(man, names=None):
    names = [n for n in (names or list(man)) if n in man]
    print("\nplots:")
    for fn in FIGURES:
        try:
            fn(man, names)
        except Exception as e:
            print(f"  {fn.__name__} FAILED: {type(e).__name__}: {e}")


if __name__ == "__main__":
    m = cmplib.load_manifest()
    make_all(m, sys.argv[1:] or None)
