#!/usr/bin/env python
"""plot_tests.py -- one figure per test, per bend angle, for the MULTIPOLET study.

Adapted from runs/bendtest/validation/plot_tests.py. For 30 and 60 deg a set of
figures test_01 .. test_11 is written to plots/<angle>deg/, each showing the
measured quantity against the analytic prediction for the MULTIPOLET case and,
where there is one, the SBEND reference. Called by run_tests.py, or run alone:

    python plot_tests.py

Run with the conda python (numpy/scipy/h5py/matplotlib).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bendlib as bl  # noqa: E402

PLOTS = HERE / "plots"
ANGLES = (30, 60)
TAGS = {"mt_fringe": "MULTIPOLET (tanh 0.02 m)", "sbend_enge": "SBEND (Enge)"}
COL = {"mt_fringe": "#d62728", "mt_sharp": "#ff7f0e", "sbend_enge": "#1f77b4"}


def _analytic_D(s, m, R16, R26):
    """Design dispersion D(s): 0 in the lead-in, sector build-up in the body,
    then a drift, in the element's own x convention."""
    s = np.atleast_1d(s).astype(float)
    D = np.zeros_like(s)
    body = (s > m["face_in_s"]) & (s <= m["face_out_s"])
    D[body] = m["x_sign"] * m["rho"] * (1 - np.cos((s[body] - m["face_in_s"]) / m["rho"]))
    aft = s > m["face_out_s"]
    D[aft] = R16 + R26 * (s[aft] - m["face_out_s"])
    return D


def _shade(ax, case):
    s0, s1 = case.field_extent()
    ax.axvspan(s0, s1, color="0.86", zorder=0, label="field region (+ fringe)")
    ax.axvspan(case.m["face_in_s"], case.m["face_out_s"], facecolor="#fdd0a2",
               edgecolor="#e6550d", lw=0.8, alpha=0.5, zorder=0, label="geometric element")


def _geom_info(man, adeg):
    c = bl.Case(f"mt_fringe_{adeg}", man)
    return (f"run:  {adeg} deg sector bend, arc L = {c.m['arc']:g} m, rho = {c.m['rho']:.4f} m, "
            f"B0 = {abs(c.m['B0']):.4f} T  ·  e-, E_kin = 100 MeV, P0 = {c.m['P0_GeV']:.4f} GeV/c, "
            f"Brho = {bl.brho(c.m['P0_GeV']):.4f} T m\n"
            f"MULTIPOLET: TP = {{{c.m['B0']:+.4f}}} T, LFRINGE = RFRINGE = {c.m['lambda']:g} m "
            f"(sharp: {man[f'mt_sharp_{adeg}']['lambda']:g} m)      "
            f"SBEND: HGAP = {man[f'sbend_enge_{adeg}']['hgap']:g} m")


def test01_energy(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for tag in TAGS:
        c = bl.Case(f"{tag}_{adeg}", man)
        E = c.stat_df["energy"].to_numpy()
        ax.plot(c.stat_df["s"], E / E[0] - 1.0, "-", color=COL[tag],
                label=f"{TAGS[tag]}  (max |d|p||/p0 = {c.momentum_conservation()['max_rel']:.1e})")
    _shade(ax, bl.Case(f"mt_fringe_{adeg}", man))
    ax.axhline(0, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("path length s [m]"); ax.set_ylabel("mean KE:  E/E(0) - 1")
    ax.set_title(f"Test 1: energy conservation ({adeg} deg)")
    ax.legend(); ax.grid(alpha=0.3)
    _save(fig, "test_01_energy", plt, man, adeg)


def test02_bend_angle(man, plt, adeg):
    """Top: accumulated bend angle from int B ds vs the orbit angle. Bottom: the
    on-axis field profile By/B0 near the two faces (fringe shapes)."""
    from scipy.integrate import cumulative_trapezoid
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for col, tag in enumerate(TAGS):
        c = bl.Case(f"{tag}_{adeg}", man)
        df = c.stat_df; s = df["s"].to_numpy()
        th = np.abs(cumulative_trapezoid(df["By_ref"].to_numpy(), s, initial=0)) / bl.brho(c.m["P0_GeV"])
        orb = np.abs(np.arctan2(df["ref_px"].to_numpy(), df["ref_pz"].to_numpy()))
        ax = axes[0, col]
        ax.plot(s, np.degrees(th), "-", color=COL[tag], lw=2, label="int B ds / Brho")
        ax.plot(s, np.degrees(orb), "--", color="k", lw=1.2, label="orbit angle atan2(px,pz)")
        ax.axhline(adeg, color="0.5", lw=0.8, ls=":", label="design ANGLE")
        _shade(ax, c)
        ax.set_title(f"{TAGS[tag]}: exit {math.degrees(orb[-1]):.4f} deg")
        ax.set_xlabel("s [m]"); ax.set_ylabel("bend angle [deg]"); ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    for tag in ("mt_fringe", "mt_sharp", "sbend_enge"):
        c = bl.Case(f"{tag}_{adeg}", man)
        s, prof = c.field_profile()
        for ax, face in zip(axes[1], (c.m["face_in_s"], c.m["face_out_s"])):
            sel = np.abs(s - face) < 0.15
            ax.plot(s[sel] - face, np.abs(prof[sel]), "-", color=COL[tag], label=tag)
    for ax, lbl in zip(axes[1], ("entrance", "exit")):
        ax.axvline(0, color="0.4", lw=0.8, ls=":")
        ax.set_xlabel(f"s - s_{lbl} [m]"); ax.set_ylabel("|By_ref| / B0"); ax.grid(alpha=0.3)
        ax.set_title(f"on-axis field at the {lbl} face"); ax.legend(fontsize=8)
    fig.suptitle(f"Test 2: bend angle from the field integral, and the fringe profiles ({adeg} deg)")
    _save(fig, "test_02_bend_angle", plt, man, adeg)


def _matrix_panels(fig, axes, c, title):
    M, _ = c.transfer_map()
    A = bl.analytic_map(c.m)
    for ax, mat, t in zip(axes, (M, A, M - A), ("measured M", "analytic M", "M - analytic")):
        vmax = np.abs(M - A).max() if t.startswith("M -") else np.abs(M).max()
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(f"{title}: {t}")
        ax.set_xticks(range(6)); ax.set_yticks(range(6))
        ax.set_xticklabels(bl.COORDS, rotation=45); ax.set_yticklabels(bl.COORDS)
        for i in range(6):
            for j in range(6):
                if abs(mat[i, j]) > 1e-4:
                    ax.text(j, i, f"{mat[i, j]:.3f}", ha="center", va="center", fontsize=6)
        fig.colorbar(im, ax=ax, shrink=0.75)


def test03_matrix(man, plt, adeg):
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for row, tag in enumerate(TAGS):
        _matrix_panels(fig, axes[row], bl.Case(f"{tag}_{adeg}", man), tag)
    fig.suptitle(f"Test 3: transfer matrix vs analytic sector bend ({adeg} deg)")
    _save(fig, "test_03_matrix", plt, man, adeg)


def test04_dispersion(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(8.5, 5))
    for tag in TAGS:
        c = bl.Case(f"{tag}_{adeg}", man)
        i_p = c.m["part_labels"].index("delta+"); i_r = c.m["part_labels"].index("ref")
        sp, xp, _ = c.trajectory(i_p)
        _, xr, _ = c.trajectory(i_r)
        Dtraj = (xp - xr) / c.m["eps"]["delta"]
        A = bl.analytic_map(c.m)
        sg = np.linspace(sp.min(), sp.max(), 400)
        step = max(1, len(sp) // 35)
        ax.plot(sp[::step], Dtraj[::step], "o", ms=5, markerfacecolor="white",
                markeredgecolor=COL[tag], markeredgewidth=1.3, zorder=3,
                label=f"{TAGS[tag]} tracked x_delta/delta")
        ax.plot(sg, _analytic_D(sg, c.m, A[0, 5], A[1, 5]), "-", color=COL[tag], lw=1.8,
                zorder=5, label=f"{TAGS[tag]} analytic D(s)")
    _shade(ax, bl.Case(f"mt_fringe_{adeg}", man))
    ax.set_xlabel("path length s [m]"); ax.set_ylabel("dispersion D [m]")
    ax.set_title(f"Test 4: dispersion orbit vs analytic ({adeg} deg)\n"
                 "MULTIPOLET's local x points towards the centre of curvature, SBEND's away")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    _save(fig, "test_04_dispersion", plt, man, adeg)


def test05_rmsz(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    c = bl.Case(f"mt_bunch_pure_{adeg}", man); df = c.stat_df
    M, _ = bl.Case(f"mt_fringe_{adeg}", man).transfer_map()
    Sin, _, _, _ = c.sigma_at_faces(); sdel = math.sqrt(Sin[5, 5])
    ax.plot(df["s"], df["rms_s"] * 1e6, "-", color=COL["mt_fringe"], label="rms_z tracked")
    ax.axhline(abs(M[4, 5]) * sdel * 1e6, color="k", ls=":", lw=1.2,
               label=f"|R56| sigma_delta = {abs(M[4,5])*sdel*1e6:.2f} um")
    _shade(ax, c)
    ax.set_xlabel("path length s [m]"); ax.set_ylabel("rms_z [um]")
    ax.set_title(f"Test 5: RMS z from R56, energy spread only ({adeg} deg)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    _save(fig, "test_05_rmsz", plt, man, adeg)


def test06_vertical(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for tag in TAGS:
        c = bl.Case(f"{tag}_{adeg}", man)
        i_y = c.m["part_labels"].index("y+")
        s, _, y = c.trajectory(i_y)
        M, _ = c.transfer_map()
        ax.plot(s, y * 1e3, "-", color=COL[tag], label=f"{TAGS[tag]}  (R43 = {M[3,2]:+.4f})")
    _shade(ax, bl.Case(f"mt_fringe_{adeg}", man))
    ax.set_xlabel("path length s [m]"); ax.set_ylabel("y [mm]")
    ax.set_title(f"Test 6: vertical plane of a sector bend stays flat ({adeg} deg)")
    ax.legend(); ax.grid(alpha=0.3)
    _save(fig, "test_06_vertical", plt, man, adeg)


def test07_symplectic(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    names, vals, colors = [], [], []
    for tag in ("mt_fringe", "mt_sharp", "sbend_enge"):
        M, _ = bl.Case(f"{tag}_{adeg}", man).transfer_map()
        names.append(tag); vals.append(bl.symplectic_residual(M)); colors.append(COL[tag])
    ax.bar(range(len(names)), vals, color=colors)
    ax.set_yscale("log"); ax.axhline(1e-4, color="k", ls="--", lw=0.8, label="tol 1e-4")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names)
    ax.set_ylabel("max|M^T J M - J|")
    ax.set_title(f"Test 7: symplecticity ({adeg} deg)")
    ax.legend()
    _save(fig, "test_07_symplectic", plt, man, adeg)


def test08_emittance(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    c = bl.Case(f"mt_bunch_emit_{adeg}", man); df = c.stat_df
    ax.plot(df["s"], df["emit_x"] / df["emit_x"].iloc[0], "-", color=COL["mt_fringe"],
            label="emit_x / emit_x(0)")
    ax.plot(df["s"], df["emit_y"] / df["emit_y"].iloc[0], "--", color=COL["mt_fringe"],
            label="emit_y / emit_y(0)")
    _shade(ax, c)
    ax.axhline(1.0, color="0.4", lw=0.8, ls=":")
    ax.set_xlabel("path length s [m]"); ax.set_ylabel("emit / emit(entrance)")
    ax.set_title(f"Test 8: emittance returns to its entrance value ({adeg} deg)\n"
                 "(the in-bend bump is a co-moving-frame transient that closes)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    _save(fig, "test_08_emittance", plt, man, adeg)


def test09_rmsx(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    c = bl.Case(f"mt_bunch_disp_{adeg}", man); df = c.stat_df
    M, _ = bl.Case(f"mt_fringe_{adeg}", man).transfer_map()
    Sin, Sout, _, _ = c.sigma_at_faces()
    pred = math.sqrt((M @ Sin @ M.T)[0, 0])
    ax.plot(df["s"], df["rms_x"] * 1e3, "-", color=COL["mt_fringe"], label="rms_x tracked")
    ax.plot(c.m["face_out_s"], pred * 1e3, "*", color="k", ms=13,
            label=f"pred sqrt((M Sig M^T)11) = {pred*1e3:.3f} mm")
    _shade(ax, c)
    ax.set_xlabel("path length s [m]"); ax.set_ylabel("rms_x [mm]")
    ax.set_title(f"Test 9: RMS x grows from dispersion ({adeg} deg)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    _save(fig, "test_09_rmsx", plt, man, adeg)


def test10_centroid(man, plt, adeg):
    fig, ax = plt.subplots(figsize=(6, 4.8))
    M, _ = bl.Case(f"mt_fringe_{adeg}", man).transfer_map(); D = M[0, 5]
    mean = bl.Case(f"mt_bunch_mean_{adeg}", man)
    meand = mean.m["bunch"]["mean_delta"]
    shift = bl.centroid_orbit_shift(mean, bl.Case(f"mt_bunch_emit_{adeg}", man),
                                    mean.m["face_out_s"])
    ax.bar([0], [D * meand * 1e3], 0.5, color=COL["mt_fringe"], alpha=0.5, label="D <delta> (analytic)")
    ax.bar([1], [shift * 1e3], 0.5, color=COL["mt_fringe"], label="measured orbit shift")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["D <delta>", "measured"])
    ax.set_ylabel("centroid <x> at exit [mm]")
    ax.set_title(f"Test 10: centroid shift from mean energy offset ({adeg} deg)")
    ax.legend(); ax.grid(alpha=0.3, axis="y")
    _save(fig, "test_10_centroid", plt, man, adeg)


def test11_combined_function(man, plt, adeg):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3))
    c = bl.Case(f"mt_cf_{adeg}", man)
    _matrix_panels(fig, axes, c, "mt_cf")
    fig.suptitle(f"Test 11: combined-function MULTIPOLET (TP[1] = {c.m['B1']:g} T/m, "
                 f"k1 = {bl.k1_of(c.m):+.3f} 1/m^2) vs analytic ({adeg} deg)")
    _save(fig, "test_11_combined_function", plt, man, adeg)


def _save(fig, name, plt, man, adeg):
    fig.text(0.5, 0.01, _geom_info(man, adeg), ha="center", va="bottom", fontsize=7,
             family="monospace",
             bbox=dict(boxstyle="round,pad=0.4", fc="#f6f6f4", ec="0.75", lw=0.7))
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    outdir = PLOTS / f"{adeg}deg"
    fig.savefig(outdir / f"{name}.png", dpi=140)
    plt.close(fig)
    print(f"  {adeg}deg/{name}.png")


_TESTS = (test01_energy, test02_bend_angle, test03_matrix, test04_dispersion, test05_rmsz,
          test06_vertical, test07_symplectic, test08_emittance, test09_rmsx, test10_centroid,
          test11_combined_function)


def make_all(man):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for adeg in ANGLES:
        (PLOTS / f"{adeg}deg").mkdir(parents=True, exist_ok=True)
        print(f"per-test plots ({adeg} deg):")
        for fn in _TESTS:
            fn(man, plt, adeg)


if __name__ == "__main__":
    make_all(bl.load_manifest())
