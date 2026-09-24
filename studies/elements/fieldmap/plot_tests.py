#!/usr/bin/env python
"""plot_tests.py -- figures for the G4beamline field map study.

Called at the end of run_tests.py, or on its own:
    python plot_tests.py

Writes plots/<name>.png. Every figure carries a footer saying which case it
came from, so a stray PNG can always be traced back.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import fmlib as F  # noqa: E402
from opalxruns import plotstyle  # noqa: E402

HERE = Path(__file__).resolve().parent
PLOTS = HERE / "plots"

COL = {"grid": "#1f77b4", "cylinder": "#d62728", "analytic": "#333333"}


def _save(fig, name, footer):
    plotstyle.save_with_footer(fig, PLOTS / f"{name}.png", footer)


def _shade(ax, case):
    """Grey the stretch of the line the map covers."""
    ax.axvspan(case.m["field_s0"], case.m["field_s1"], color="0.9", zorder=0,
               label="map extent")


def plot_field_profiles(cases):
    """On-axis field the reference particle sees, for one case per field."""
    picks = [("grid_dipole", "By_ref", "3D grid, uniform By"),
             ("grid_sol", "Bz_ref", "3D grid, uniform Bz"),
             ("cyl_sol", "Bz_ref", "2D cylinder, uniform Bz"),
             ("cyl_ramp", "Bz_ref", "2D cylinder, asymmetric Bz")]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    for ax, (name, col, title) in zip(axes.ravel(), picks):
        c = cases[name]
        df = c.stat_df
        _shade(ax, c)
        ax.plot(df["s"], df[col], color=COL[c.m["fmt"]], lw=1.2)
        ax.set_title(f"{title}  ({name})", fontsize=9)
        ax.set_ylabel(f"{col} [T]", fontsize=8)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("s [m]", fontsize=8)
    fig.suptitle("Field on the reference orbit: the map sits where it was placed",
                 fontsize=11)
    _save(fig, "01_field_profiles", "from the .stat files of grid_dipole, grid_sol, "
          "cyl_sol, cyl_ramp")


def plot_dipole(cases):
    """The bent reference orbit against the circular arc."""
    c = cases["grid_dipole"]
    s, x, _, _, _ = c.ref_orbit()
    s0, s1 = c.m["field_s0"], c.m["field_s1"]
    theta, offset = F.dipole_exit(F.RHO_DIP, F.L_FIELD)

    # The arc, against path length rather than z: at path length sigma from the
    # entrance the orbit has turned through sigma/rho. The arc is 1.0051 m long
    # while the field is 1 m deep in z, so the exit face sits at s0 + rho*theta,
    # not at s1.
    s_exit = s0 + F.RHO_DIP * theta
    ss = np.linspace(s0, s_exit, 400)
    arc = F.RHO_DIP * (1.0 - np.cos((ss - s0) / F.RHO_DIP))
    after = np.linspace(s_exit, s[-1], 100)
    # Against path length the straight section rises as sin(theta), not
    # tan(theta): the path is the hypotenuse, not the z step.
    line = offset + (after - s_exit) * math.sin(theta)

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    # Shaded against path length, so the exit face sits at s0 + rho*theta.
    for a in (ax, ax2):
        a.axvspan(s0, s_exit, color="0.9", zorder=0)
    ax.axvspan(np.nan, np.nan, color="0.9", label="map extent")
    ax.plot(s, np.abs(x), color=COL["grid"], lw=1.5, label="OPALX reference orbit")
    ax.plot(ss, arc, "--", color=COL["analytic"], lw=1.2,
            label="arc of radius rho")
    ax.axvline(s_exit, color="0.5", lw=0.8, ls=":")
    ax.plot(after, line, ":", color=COL["analytic"], lw=1.2, label="straight after it")
    ax.set_ylabel("|x| [m]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title(f"Uniform By = {F.B_DIP:.6f} T from a 3D grid map: "
                 f"{F.BEND_DEG:g} deg arc, rho = {F.RHO_DIP:.4f} m", fontsize=10)

    model_s = np.concatenate([[0.0], ss, after])
    model_x = np.concatenate([[0.0], arc, line])
    resid = np.abs(x) - np.interp(s, model_s, model_x)
    ax2.plot(s, resid * 1e6, color=COL["grid"], lw=1.0)
    ax2.set_ylabel("difference [um]")
    ax2.set_xlabel("s [m]")
    ax2.grid(alpha=0.3)
    _save(fig, "02_dipole_orbit", "grid_dipole: reference orbit vs the arc")


def _matrix_figure(name, M, A, title, footer):
    labels = ["x", "x'", "y", "y'"]
    idx = [(i, j) for i in range(4) for j in range(4)]
    keep = [(i, j) for i, j in idx if abs(A[i, j]) > 1e-12 or abs(M[i, j]) > 1e-12]
    names = [f"R{i+1}{j+1}" for i, j in keep]
    mv = [M[i, j] for i, j in keep]
    av = [A[i, j] for i, j in keep]
    pos = np.arange(len(keep))

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9, 5.5),
                                  gridspec_kw={"height_ratios": [2, 1]})
    ax.bar(pos - 0.2, av, 0.4, color=COL["analytic"], label="closed form")
    ax.bar(pos + 0.2, mv, 0.4, color=COL["grid"], label="from the map")
    ax.set_xticks(pos)
    ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("matrix element")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    ax.set_title(title, fontsize=10)

    rel = [(m - a) / a if abs(a) > 1e-12 else 0.0 for m, a in zip(mv, av)]
    ax2.bar(pos, np.abs(rel), 0.5, color="#ff7f0e")
    ax2.set_yscale("log")
    ax2.set_xticks(pos)
    ax2.set_xticklabels(names, fontsize=8)
    ax2.set_ylabel("|relative\ndifference|", fontsize=8)
    ax2.grid(alpha=0.3, axis="y")
    ax2.set_xlabel(f"transverse coordinates {labels}", fontsize=8)
    _save(fig, name, footer)


def plot_quad(cases):
    c = cases["grid_quad"]
    A = F.quad_matrix(F.k1_quad(), F.L_FIELD)
    _matrix_figure("03_quad_matrix", c.transfer_map4(), A,
                   f"Quadrupole from a 3D grid map, k1 = {F.k1_quad():.4f} m^-2",
                   "grid_quad: 4x4 transverse matrix vs the thick-quadrupole matrix")


def plot_solenoid(cases):
    c = cases["grid_sol"]
    A = F.solenoid_matrix(F.k_solenoid(), F.L_FIELD)
    _matrix_figure("04_solenoid_matrix", c.transfer_map4(), A,
                   f"Uniform Bz = {F.B_SOL:g} T from a 3D grid map, "
                   f"k L = {F.k_solenoid() * F.L_FIELD:.4f} rad",
                   "grid_sol: 4x4 transverse matrix vs the uniform-Bz matrix")


def plot_formats_agree(cases):
    """The two formats holding the same field, and what is left between them."""
    g, c = cases["grid_sol"], cases["cyl_sol"]
    dg = np.abs(g.transfer_map4() - c.transfer_map4())

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    for case, style in ((g, "-"), (c, "--")):
        df = case.stat_df
        ax.plot(df["s"], df["Bz_ref"], style, lw=1.4, color=COL[case.m["fmt"]],
                label=f"{case.m['fmt']} ({case.name})")
    _shade(ax, g)
    ax.set_xlabel("s [m]")
    ax.set_ylabel("Bz on the reference orbit [T]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("The same uniform Bz, both formats", fontsize=10)

    lab = ["x", "x'", "y", "y'"]
    if dg.max() == 0.0:
        # Nothing to colour: every element came out bit for bit the same.
        ax2.axis("off")
        ax2.text(0.5, 0.58, "every element of the 4x4\nidentical to the last digit",
                 ha="center", va="center", fontsize=12, color="#1a7f37")
        ax2.text(0.5, 0.30, "the two readers put the same field\nin the same place",
                 ha="center", va="center", fontsize=9, color="#555555")
        ax2.set_title("Difference between the two 4x4 matrices", fontsize=10)
    else:
        im = ax2.imshow(np.log10(np.maximum(dg, 1e-18)), cmap="viridis")
        ax2.set_xticks(range(4)); ax2.set_yticks(range(4))
        ax2.set_xticklabels(lab); ax2.set_yticklabels(lab)
        for i in range(4):
            for j in range(4):
                ax2.text(j, i, f"{dg[i, j]:.0e}", ha="center", va="center",
                         fontsize=7, color="w")
        ax2.set_title("log10 |difference| of the two 4x4 matrices", fontsize=10)
        fig.colorbar(im, ax=ax2, shrink=0.8)
    _save(fig, "05_formats_agree", "grid_sol vs cyl_sol: the same field in the two "
          "G4beamline formats")


def plot_zreverse(cases):
    """The asymmetric profile, reversed by the reader and reversed in the file."""
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    styles = {"cyl_ramp": ("-", "#1f77b4", "as written"),
              "cyl_ramp_zrev": ("--", "#d62728", "ZREVERSE = TRUE"),
              "cyl_ramp_mirror": (":", "#2ca02c", "reversed in the file")}
    for name, (ls, col, lab) in styles.items():
        df = cases[name].stat_df
        ax.plot(df["s"], df["Bz_ref"], ls, color=col, lw=1.6, label=lab)
    _shade(ax, cases["cyl_ramp"])
    ax.set_ylabel("Bz on the reference orbit [T]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("Turning an asymmetric solenoid map round", fontsize=10)

    a, b = cases["cyl_ramp_zrev"], cases["cyl_ramp_mirror"]
    sa = a.stat_df["s"].to_numpy()
    d = a.stat_df["Bz_ref"].to_numpy() - np.interp(
        sa, b.stat_df["s"].to_numpy(), b.stat_df["Bz_ref"].to_numpy())
    ax2.plot(sa, np.abs(d), color="#ff7f0e", lw=1.0)
    if np.abs(d).max() > 0:
        ax2.set_yscale("log")
    else:
        ax2.set_ylim(-1e-18, 1e-18)
        ax2.text(0.5, 0.5, "identical to the last digit", ha="center", va="center",
                 transform=ax2.transAxes, fontsize=9, color="#555555")
    ax2.set_ylabel("|difference| [T]", fontsize=8)
    ax2.set_xlabel("s [m]")
    ax2.grid(alpha=0.3)
    _save(fig, "06_zreverse", "cyl_ramp / cyl_ramp_zrev / cyl_ramp_mirror")


def make_all(manifest=None):
    man = manifest or F.load_manifest()
    cases = {n: F.Case(n, man) for n in man}
    for fn in (plot_field_profiles, plot_dipole, plot_quad, plot_solenoid,
               plot_formats_agree, plot_zreverse):
        fn(cases)
    print(f"wrote {len(list(PLOTS.glob('*.png')))} figures to {PLOTS}")


if __name__ == "__main__":
    make_all()
