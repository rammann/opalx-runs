#!/usr/bin/env python
"""run_tests.py -- the MULTIPOLET bend-study tests.

Adapted from runs/bendtest/validation/run_tests.py. Reads the tracked output of
the generated cases (see make_decks.py), prints every measured quantity next to
its analytic value with the difference, and PASSES/FAILS against a tolerance.
Writes diagnostic plots to output/elements/multipolet/plots/. Exit code is non-zero if any test fails.

Cases
  mt_fringe_{30,60,90}   MULTIPOLET, tanh fringe 0.02 m      -- primary
  mt_sharp_{30,60,90}    MULTIPOLET, tanh fringe 0.005 m     -- sharper edge
  sbend_enge_{30,60,90}  SBEND, Enge fringe HGAP = 0.01 m    -- reference bend
  mt_cf_{30,60}          MULTIPOLET with a gradient TP[1]    -- combined function
  mt_bunch_*_{30,60}     Gaussian bunches for tests 5, 8-10

Tests
   1 energy conservation        |p| per particle constant first -> last dump
   2 bend angle                 theta from int By ds vs the orbit exit angle,
                                and the orbit exit angle vs the design ANGLE
   3 transfer matrix            bend-plane block, dispersion, vertical block
   4 dispersion                 R16 = +/- rho (1 - cos), R26 = +/- sin
   5 RMS z from R56             energy-spread-only bunch
   6 vertical plane             sector bend: R43 ~ 0, R34 ~ L, MULTIPOLET vs SBEND
   7 symplecticity              M^T J M = J, block determinants
   8 emittance invariance       eps_x, eps_y, 6D emittance in vs out
   9 RMS x from energy spread   sigma_x^2 = betatron + (D sigma_delta)^2
  10 centroid from mean energy  orbit shift = D <delta>
  11 combined function          full map vs the combined-function sector matrix

Run with the conda python (numpy/scipy/h5py/matplotlib). See run_all.sh.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bendlib as bl

TOL = {
    "energy": 1e-6,          # |dp|/p0
    "field_angle": 2e-3,     # rad, |theta_field - exit_angle|
    "design_angle": 1e-3,    # rad, |exit_angle - ANGLE|
    "matrix": 3e-2,          # relative, matrix elements (fringe-aware)
    "matrix_abs": 2.0e-2,    # absolute floor for near-zero matrix elements (= R43_flat)
    "disp": 2e-2,            # relative, dispersion
    "symplectic": 1e-4,      # max|M^T J M - J|
    # |R43| of a sector bend, as a fraction of the rectangular-bend edge focusing
    # tan(theta/2)/rho. A thin-edge sector bend has no vertical focusing at all, but a
    # real extended fringe leaves a small residual that grows with the bend angle, so the
    # bound has to be relative to that scale rather than a fixed number. Measured: about
    # 1 % for the tanh fringe, 2.5 % for the Enge fringe, at every angle.
    "R43_flat": 5e-2,
    "emit": 5e-3,            # relative emittance change
    "emit6d": 1e-2,          # relative 6D emittance change
    "cov": 2e-2,             # relative, RMS growth predictions
    "centroid": 1e-2,        # relative, centroid shift
}

FULL_ANGLES = (30, 60)
XCHECK_ANGLES = (90,)          # map tests only
ALL_ANGLES = FULL_ANGLES + XCHECK_ANGLES
MAP_TAGS = ("mt_fringe", "mt_sharp", "sbend_enge")


def _mtol(a: float) -> float:
    return max(TOL["matrix"] * abs(a), TOL["matrix_abs"])


def _r43_tol(case) -> float:
    """How much vertical focusing a fringe may leave in a sector bend: a fraction of
    the rectangular-bend edge focusing tan(theta/2)/rho of the same geometry."""
    return TOL["R43_flat"] * math.tan(0.5 * case.m["angle"]) / case.m["rho"]


# ---------------------------------------------------------------------------
# Single-particle tests
# ---------------------------------------------------------------------------

def test1_energy(R, case):
    mc = case.momentum_conservation()
    R.check("1", f"|p| conserv [{case.name}]", mc["max_rel"], 0.0, TOL["energy"],
            note="magnetic field does no work")


def test2_bend_angle(R, case):
    theta_field = abs(case.field_integral()) / bl.brho(case.m["P0_GeV"])
    exit_ang = bl.exit_bend_angle(case)
    R.check("2", f"theta from int B ds [{case.name}]", theta_field, exit_ang,
            TOL["field_angle"], note="vs orbit exit angle")
    R.check("2", f"exit angle vs ANGLE [{case.name}]", exit_ang, case.m["angle"],
            TOL["design_angle"],
            note=f"orbit turns towards {'+' if bl.bend_direction(case) > 0 else '-'}x")


def test3_matrix(R, case):
    M, _ = case.transfer_map()
    A = bl.analytic_map(case.m)
    for nm, i, j in (("R11", 0, 0), ("R12", 0, 1), ("R21", 1, 0), ("R22", 1, 1),
                     ("R16", 0, 5), ("R26", 1, 5), ("R33", 2, 2), ("R34", 2, 3),
                     ("R44", 3, 3)):
        R.check("3", f"{nm} [{case.name}]", M[i, j], A[i, j], _mtol(A[i, j]))
    # R43 is zero in the thin-edge sector matrix, so it gets the fringe-scaled bound.
    R.check("3", f"R43 [{case.name}]", M[3, 2], A[3, 2], _r43_tol(case),
            note="fringe residual, bound relative to tan(theta/2)/rho")


def test4_dispersion(R, case):
    M, _ = case.transfer_map()
    A = bl.analytic_map(case.m)
    R.check("4", f"x  = D delta  [{case.name}]", M[0, 5], A[0, 5], TOL["disp"], rel=True,
            note=f"D = {case.m['x_sign']:+d} rho(1-cos)")
    R.check("4", f"x' = D' delta [{case.name}]", M[1, 5], A[1, 5], TOL["disp"], rel=True)


def test6_vertical(R, mt, sb):
    Mm, _ = mt.transfer_map()
    Ms, _ = sb.transfer_map()
    L = mt.m["arc"]
    R.check("6", f"R43 [{mt.name}]", Mm[3, 2], 0.0, _r43_tol(mt),
            note="sector: no vertical focusing beyond the fringe residual")
    R.check("6", f"R43 [{sb.name}]", Ms[3, 2], 0.0, _r43_tol(sb))
    R.check("6", f"R34 = L [{mt.name}]", Mm[2, 3], L, TOL["matrix"], rel=True)
    R.check("6", f"R34 = L [{sb.name}]", Ms[2, 3], L, TOL["matrix"], rel=True)
    scale = math.tan(0.5 * mt.m["angle"]) / mt.m["rho"]
    R.note("6", "MULTIPOLET vs SBEND R43",
           f"{Mm[3, 2]:+.4f} vs {Ms[3, 2]:+.4f}  (tan(theta/2)/rho = {scale:.4f})")


def test7_symplectic(R, case, diagnostic=False):
    M, _ = case.transfer_map()
    res = bl.symplectic_residual(M)
    tol = None if diagnostic else TOL["symplectic"]
    R.check("7", f"symplectic res [{case.name}]", res, 0.0, tol)
    for lbl, d in zip(("x", "y", "z"), bl.block_dets(M)):
        R.check("7", f"det {lbl}-block [{case.name}]", d, 1.0, tol)


def test11_combined_function(R, case):
    M, _ = case.transfer_map()
    A = bl.analytic_map(case.m)
    k1 = bl.k1_of(case.m)
    R.note("11", f"[{case.name}]", f"k1 = {k1:+.4f} 1/m^2, kx = 1/rho^2 + k1 = "
           f"{1/case.m['rho']**2 + k1:+.4f}, ky = {-k1:+.4f}")
    for nm, i, j in (("R11", 0, 0), ("R12", 0, 1), ("R21", 1, 0), ("R22", 1, 1),
                     ("R33", 2, 2), ("R34", 2, 3), ("R43", 3, 2), ("R44", 3, 3),
                     ("R16", 0, 5), ("R26", 1, 5), ("R56", 4, 5)):
        R.check("11", f"{nm} [{case.name}]", M[i, j], A[i, j], _mtol(A[i, j]))
    R.check("11", f"symplectic res [{case.name}]", bl.symplectic_residual(M), 0.0,
            TOL["symplectic"])
    # A source-free gradient focuses in one plane and defocuses in the other, so the
    # vertical block has to follow ky = -k1. Reported next to the wrong-sign value a
    # mirrored transverse field would give, since that was the defect here.
    Cy, _ = bl._cs(k1, case.m["arc"])
    R.note("11", f"R33 [{case.name}]",
           f"{M[2, 2]:+.4f}; ky = -k1 gives {bl._cs(-k1, case.m['arc'])[0]:+.4f}, "
           f"ky = +k1 would give {Cy:+.4f}")


# ---------------------------------------------------------------------------
# Bunch statistics tests
# ---------------------------------------------------------------------------

def test5_rmsz(R, bunch_pure, mapcase):
    M, _ = mapcase.transfer_map()
    R56, D = M[4, 5], M[0, 5]
    Sin, Sout, _, _ = bunch_pure.sigma_at_faces()
    sdel = math.sqrt(Sin[5, 5])
    pred_z = math.sqrt(Sin[4, 4] + (R56 * sdel) ** 2)
    R.check("5", f"sigma_z = |R56| sd [{bunch_pure.name}]", math.sqrt(Sout[4, 4]), pred_z,
            TOL["cov"], rel=True, note=f"R56={R56:.4f}")
    R.check("5", "sigma_x = |D| sd (pure)", math.sqrt(Sout[0, 0]), abs(D) * sdel,
            TOL["cov"], rel=True)


def test8_emittance(R, bunch_emit):
    Sin, Sout, _, _ = bunch_emit.sigma_at_faces()
    for p, nm in ((0, "x"), (1, "y")):
        R.check("8", f"emit_{nm} [{bunch_emit.name}]", bl.emittance(Sout, p),
                bl.emittance(Sin, p), TOL["emit"], rel=True)
    R.check("8", f"emit_6D [{bunch_emit.name}]", bl.emittance_6d(Sout), bl.emittance_6d(Sin),
            TOL["emit6d"], rel=True, note="invariant of any symplectic map")


def test9_rmsx(R, bunch_disp, mapcase):
    M, _ = mapcase.transfer_map()
    D = M[0, 5]
    Sin, Sout, _, _ = bunch_disp.sigma_at_faces()
    sdel = math.sqrt(Sin[5, 5])
    Spred = M @ Sin @ M.T
    R.check("9", f"sigma_x (M Sig M^T) [{bunch_disp.name}]", math.sqrt(Sout[0, 0]),
            math.sqrt(Spred[0, 0]), TOL["cov"], rel=True)
    Sin_bet = Sin.copy(); Sin_bet[5, :] = 0; Sin_bet[:, 5] = 0
    pred = math.sqrt((M @ Sin_bet @ M.T)[0, 0] + (D * sdel) ** 2)
    R.check("9", "sigma_x = sqrt(bet+(D sd)^2)", math.sqrt(Sout[0, 0]), pred, TOL["cov"],
            rel=True, note=f"D={D:.4f} sigma_delta={sdel:.2e}")


def test10_centroid(R, bunch_mean, bunch_design, mapcase):
    M, _ = mapcase.transfer_map()
    D = M[0, 5]
    meand = bunch_mean.m["bunch"]["mean_delta"]
    shift = bl.centroid_orbit_shift(bunch_mean, bunch_design, bunch_mean.m["face_out_s"])
    R.check("10", f"<x> = D <delta> [{bunch_mean.name}]", shift, D * meand, TOL["centroid"],
            rel=True, note=f"<delta>={meand:.1e} D={D:.4f}")


# ---------------------------------------------------------------------------

def main():
    man = bl.load_manifest()
    R = bl.Results()
    C = lambda n: bl.Case(n, man)  # noqa: E731

    maps = {t: {a: C(f"{t}_{a}") for a in ALL_ANGLES} for t in MAP_TAGS}
    cf = {a: C(f"mt_cf_{a}") for a in FULL_ANGLES}
    map_cases = [maps[t][a] for a in ALL_ANGLES for t in MAP_TAGS]

    for c in map_cases + list(cf.values()):
        test1_energy(R, c)
    for c in map_cases:
        test2_bend_angle(R, c)
    for c in map_cases:
        test3_matrix(R, c)
    for c in map_cases:
        test4_dispersion(R, c)
    for a in FULL_ANGLES:
        test5_rmsz(R, C(f"mt_bunch_pure_{a}"), maps["mt_fringe"][a])
    for a in ALL_ANGLES:
        test6_vertical(R, maps["mt_fringe"][a], maps["sbend_enge"][a])
    for c in map_cases:
        test7_symplectic(R, c)
    for a in FULL_ANGLES:
        test8_emittance(R, C(f"mt_bunch_emit_{a}"))
        test9_rmsx(R, C(f"mt_bunch_disp_{a}"), maps["mt_fringe"][a])
        test10_centroid(R, C(f"mt_bunch_mean_{a}"), C(f"mt_bunch_emit_{a}"),
                        maps["mt_fringe"][a])
    for a in FULL_ANGLES:
        test11_combined_function(R, cf[a])

    import contextlib, io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        R.print_table("MULTIPOLET (curved) bend study -- measured vs analytic")
        ntot = sum(1 for r in R.rows if r["ok"] is not None)
        print(f"\n{'='*100}")
        print(f"SUMMARY: {ntot - R.failed}/{ntot} checks passed, {R.failed} failed.")
        print(f"{'='*100}")
    print(buf.getvalue())
    (Path(__file__).resolve().parent / "results.txt").write_text(buf.getvalue())

    try:
        import plot_tests
        plot_tests.make_all(man)
    except Exception as e:  # plotting must never mask a test failure
        import traceback
        traceback.print_exc()
        print(f"(plotting skipped: {e})")

    return 1 if R.failed else 0


if __name__ == "__main__":
    sys.exit(main())
