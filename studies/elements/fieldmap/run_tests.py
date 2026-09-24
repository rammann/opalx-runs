#!/usr/bin/env python
"""run_tests.py -- the checks for the G4beamline field map study.

Tracking every case first:  ./run_all.sh
Re-running just these:      ./run_all.sh --test-only

The checks fall into two kinds, and they are worth telling apart.

Checks against a closed form (tests 3, 5, 6, 7) ask whether OPALX moves a
particle the way the field in the map says it should. They cannot be exact:
the map has hard ends, and the time step resolves where those ends are only to
within one step, so a small residual is expected and its size is set by the
step rather than by the reader. Test 12 is what says so, by running the same
case with twice the step and watching the residual grow with it.

Checks between two runs (tests 4, 8, 9, 10, 11) ask whether two ways of saying
the same thing give the same answer -- the two formats, the two placement
conventions, the three ways of scaling a map, the tidy file against the awkward
one, ZREVERSE against a map written out already reversed. Both sides have the
same hard ends and the same step, so those cancel and the agreement should be
at round-off. These are the sharp tests.

Exit code is non-zero if any check fails.
"""

from __future__ import annotations

import io
import sys
import traceback
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

import fmlib as F

HERE = Path(__file__).resolve().parent

TOL = {
    # A magnetic field does no work.
    "momentum": 1e-12,
    # Where the reference particle first and last sees field, against the
    # placement. The stat file is written every 50 steps, which is 1.5 mm, so
    # this cannot be tighter than that.
    "window": 5e-3,          # m
    # Reference orbit of the uniform dipole against the arc.
    "angle": 1e-4,           # rad
    "offset": 2e-4,          # m
    # Transfer matrix against the closed form. Measured residual is about
    # 5e-6 relative, set by the step resolving the two hard map ends.
    "matrix": 5e-5,          # relative, on elements above the floor
    "matrix_abs": 1e-5,      # absolute, for elements that should be zero
    # A quadrupole is divergence free and curl free, so its matrix in these
    # coordinates is symplectic. The uniform-Bz cases are not checked for this:
    # a box of uniform Bz has no radial field at its ends, is not a physical
    # field, and its slope coordinates are not canonical.
    "symplectic": 1e-7,
    # Two runs that should be doing exactly the same arithmetic.
    "identical": 1e-11,
    # Two runs of the same field through different reader code paths, or with
    # the scaling applied at a different point, so the rounding differs.
    "equivalent": 1e-8,
    # Reversing an asymmetric map has to actually change the tracking, or the
    # comparison against the pre-reversed file would prove nothing.
    "must_differ": 1e-2,
    # Energy picked up crossing the electric field, against q E dz. Same story as the
    # matrix: limited by the step resolving where the map's two ends are.
    "energy": 2e-4,          # relative
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def matrix_rel(measured: np.ndarray, analytic: np.ndarray) -> float:
    """Largest relative difference over the elements that are not essentially
    zero. The zero elements are covered by matrix_abs instead."""
    big = np.abs(analytic) > 0.01
    if not big.any():
        return 0.0
    return float((np.abs(measured - analytic)[big] / np.abs(analytic)[big]).max())


def matrix_abs(measured: np.ndarray, analytic: np.ndarray) -> float:
    """Largest absolute difference over the elements that should be zero."""
    small = np.abs(analytic) <= 0.01
    if not small.any():
        return 0.0
    return float(np.abs(measured - analytic)[small].max())


def state_diff(a: F.Case, b: F.Case) -> float:
    """Largest difference between the two runs' final phase space. Only
    meaningful where both stopped at the same place."""
    return float(np.abs(a.final_state() - b.final_state()).max())


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test1_momentum(R, cases):
    """A magnetic field does no work, so |p| must not move. The cases with a live
    electric field are left out: there |p| is supposed to change, and test 14 is what
    checks by how much."""
    for c in cases.values():
        if not c.m["magnetic_only"]:
            continue
        R.check("1", f"|p| conserved [{c.name}]", c.momentum_error(), 0.0,
                TOL["momentum"])


def test2_window(R, cases):
    """The field must start and stop where the placement says it does."""
    for c in cases.values():
        if not c.m["sharp_edges"]:
            continue  # nothing to find: see the note on sharp_edges in make_inputs.py
        s0, s1 = c.field_window()
        R.check("2", f"field starts [{c.name}]", s0, c.m["field_s0"], TOL["window"])
        R.check("2", f"field ends [{c.name}]", s1, c.m["field_s1"], TOL["window"])


def test3_dipole(R, cases):
    """Uniform By in the 3D grid format against the circular arc."""
    c = cases["grid_dipole"]
    theta, offset = F.dipole_exit(F.RHO_DIP, F.L_FIELD)
    R.check("3", "bend angle [grid_dipole]", abs(c.exit_angle()), theta, TOL["angle"],
            note=f"{F.BEND_DEG:g} deg arc of radius {F.RHO_DIP:.4f} m")
    # Measured where the orbit reaches the exit face in lab z. Reading it at
    # that path length instead would be wrong by 5 mm, because the arc is
    # 1.0051 m long while the field is 1 m deep in z.
    R.check("3", "offset at field exit", abs(c.ref_x_at_z(c.m["field_s1"])), offset,
            TOL["offset"])


def test4_placement(R, cases):
    """Posing the same map further downstream moves the field and nothing else."""
    a, b = cases["grid_dipole"], cases["grid_dipole_shift"]
    R.check("4", "field start moves by Z_SHIFT",
            b.field_window()[0] - a.field_window()[0], F.Z_SHIFT, TOL["window"])
    R.check("4", "field end moves by Z_SHIFT",
            b.field_window()[1] - a.field_window()[1], F.Z_SHIFT, TOL["window"])
    # Not exact: the shift is not a whole number of steps, so the two runs
    # resolve the map ends at a slightly different point within a step.
    R.check("4", "bend angle unchanged by the shift",
            abs(b.exit_angle()), abs(a.exit_angle()), TOL["angle"])


def test5_quad(R, cases):
    """Quadrupole in the 3D grid format against the thick-quadrupole matrix."""
    for name, gradient in (("grid_quad", F.G_QUAD), ("grid_quad_2x", 2 * F.G_QUAD)):
        c = cases[name]
        A = F.quad_matrix(F.k1_quad(gradient), F.L_FIELD)
        M = c.transfer_map4()
        R.check("5", f"4x4 vs analytic [{name}]", matrix_rel(M, A), 0.0, TOL["matrix"],
                note=f"k1 = {F.k1_quad(gradient):.4f} m^-2")
        R.check("5", f"4x4 zero elements [{name}]", matrix_abs(M, A), 0.0,
                TOL["matrix_abs"])


def test6_symplectic(R, cases):
    """A quadrupole field is divergence free and curl free, so its matrix in
    these coordinates is symplectic."""
    for name in ("grid_quad", "grid_quad_2x"):
        R.check("6", f"symplecticity [{name}]",
                F.symplectic_residual(cases[name].transfer_map4()), 0.0,
                TOL["symplectic"])


def test7_solenoid(R, cases):
    """Uniform Bz, in both formats, against the closed form for that field."""
    A = F.solenoid_matrix(F.k_solenoid(), F.L_FIELD)
    for name in ("grid_sol", "cyl_sol", "cyl_sol_solenoid"):
        c = cases[name]
        M = c.transfer_map4()
        R.check("7", f"4x4 vs analytic [{name}]", matrix_rel(M, A), 0.0, TOL["matrix"],
                note=f"k L = {F.k_solenoid() * F.L_FIELD:.6f} rad")
        R.check("7", f"4x4 zero elements [{name}]", matrix_abs(M, A), 0.0,
                TOL["matrix_abs"])


def test8_formats_agree(R, cases):
    """The same uniform Bz written in both formats has to track the same. This
    is the sharp comparison of the two readers: both sides have the same hard
    ends and the same step, so the only thing left is the reading."""
    g, c = cases["grid_sol"], cases["cyl_sol"]
    R.check("8", "grid vs cylinder, 4x4",
            float(np.abs(g.transfer_map4() - c.transfer_map4()).max()), 0.0,
            TOL["equivalent"], note="same field, the two formats")
    R.check("8", "grid vs cylinder, final state", state_diff(g, c), 0.0,
            TOL["equivalent"])


def test9_placement_conventions(R, cases):
    """The same cylinder map on a posed FIELDMAP and on a SOLENOID placed by
    ELEMEDGE. Both put the map's own z = 0 at the same lab position, and KS on
    the solenoid is the same plain multiplier as SCALE."""
    a, b = cases["cyl_sol"], cases["cyl_sol_solenoid"]
    R.check("9", "FIELDMAP vs SOLENOID, 4x4",
            float(np.abs(a.transfer_map4() - b.transfer_map4()).max()), 0.0,
            TOL["identical"])
    R.check("9", "FIELDMAP vs SOLENOID, final state", state_diff(a, b), 0.0,
            TOL["identical"])


def test10_scaling(R, cases):
    """Three ways of asking for twice the field: doubling every value in the
    file, SCALE = 2 on the element, and normB/current = 2 on the param line."""
    ref = cases["grid_quad_2x"]
    for name, how in (("grid_quad_scale2", "SCALE = 2 on the element"),
                      ("grid_quad_param", "normB = 4, current = 2 in the file")):
        c = cases[name]
        R.check("10", f"{name} vs grid_quad_2x, 4x4",
                float(np.abs(c.transfer_map4() - ref.transfer_map4()).max()), 0.0,
                TOL["equivalent"], note=how)
        R.check("10", f"{name} vs grid_quad_2x, state", state_diff(c, ref), 0.0,
                TOL["equivalent"])


def test11_awkward_file(R, cases):
    """A comment block in front, the rows out of order, and nine columns with a
    zero electric field: the same field, so the same tracking."""
    a, b = cases["grid_quad"], cases["grid_quad_messy"]
    R.check("11", "tidy vs awkward file, 4x4",
            float(np.abs(a.transfer_map4() - b.transfer_map4()).max()), 0.0,
            TOL["identical"], note="comments, shuffled rows, nine columns")
    R.check("11", "tidy vs awkward file, state", state_diff(a, b), 0.0,
            TOL["identical"])


def test12_step_size(R, cases):
    """Where the residual against the closed form comes from. Doubling the step
    should roughly double it, which says it is the step resolving the two hard
    map ends and not something the reader did."""
    A = F.solenoid_matrix(F.k_solenoid(), F.L_FIELD)
    fine = matrix_rel(cases["grid_sol"].transfer_map4(), A)
    coarse = matrix_rel(cases["grid_sol_dt2"].transfer_map4(), A)
    R.check("12", "residual grows with the step", coarse, fine, None,
            note=f"{F.DT_FINE:g} s -> {fine:.3g},  {F.DT_COARSE:g} s -> {coarse:.3g}")
    R.at_least("12", "coarse residual is the larger", coarse, fine,
               note="the residual follows the step, so it is not the reading")
    # Part of the residual does not move with the step: it is where the two
    # hard ends fall inside a step, which also shifts the drifts stripped off
    # when the matrix is built. Letting the ends move by less than one step
    # takes grid_sol from 2.5e-05 down to 4.6e-06. Reported, not checked.
    R.check("12", "ratio of the two residuals", coarse / fine if fine else 0.0, 0.0,
            None, note="not a clean factor of 2: see the comment in test12")


def test13_zreverse(R, cases):
    """ZREVERSE against the same field written out already reversed, and a
    check that reversing this map changes anything at all."""
    a, b = cases["cyl_ramp_zrev"], cases["cyl_ramp_mirror"]
    R.check("13", "ZREVERSE vs pre-reversed file, 4x4",
            float(np.abs(a.transfer_map4() - b.transfer_map4()).max()), 0.0,
            TOL["equivalent"], note="mirror in z, negate Bz, leave Br")
    R.check("13", "ZREVERSE vs pre-reversed file, state", state_diff(a, b), 0.0,
            TOL["equivalent"])
    plain = cases["cyl_ramp"]
    R.at_least("13", "reversing changes the tracking",
               float(np.abs(a.transfer_map4() - plain.transfer_map4()).max()),
               TOL["must_differ"],
               note="otherwise the profile was symmetric and this proves nothing")


def test14_efield(R, cases):
    """The electric half of the nine-column form: is it read, converted from MV/m to V/m,
    and applied? The energy a particle gains crossing a uniform Ez is q E dz, and because
    the field has no z dependence that only depends on the z it crossed -- so the same
    number applies to the bent case, where a magnetic field is acting at the same time."""
    expected = F.energy_gain(F.E_LONG, F.L_FIELD)
    for name in ("grid_efield", "grid_eb"):
        c = cases[name]
        R.check("14", f"energy change [{name}]", c.ref_energy_change(), expected,
                TOL["energy"], rel=True,
                note=f"q E dz = {expected * 1e3:.4f} MeV")

    # ESCALE is a plain multiplier, so twice the scale is twice the energy.
    twice = cases["grid_efield_escale2"].ref_energy_change()
    R.check("14", "ESCALE = 2 doubles the energy change", twice,
            2.0 * cases["grid_efield"].ref_energy_change(), TOL["energy"], rel=True)

    # The sharp one: the same file with the electric half switched off has to track like
    # the magnet-only map. That is what says the two scales really are independent.
    a, b = cases["grid_eb_bonly"], cases["grid_dipole"]
    R.check("14", "ESCALE = 0 reproduces grid_dipole",
            float(np.abs(a.transfer_map4() - b.transfer_map4()).max()), 0.0,
            TOL["equivalent"], note="same map, electric half switched off")
    R.check("14", "ESCALE = 0, final state", state_diff(a, b), 0.0, TOL["equivalent"])
    R.check("14", "ESCALE = 0 does no work", a.momentum_error(), 0.0, TOL["momentum"])


TESTS = [test1_momentum, test2_window, test3_dipole, test4_placement, test5_quad,
         test6_symplectic, test7_solenoid, test8_formats_agree,
         test9_placement_conventions, test10_scaling, test11_awkward_file,
         test12_step_size, test13_zreverse, test14_efield]


def main() -> int:
    man = F.load_manifest()
    missing = [n for n in man if not (HERE / n / man[n]["h5"]).exists()]
    if missing:
        print("ERROR: no tracking output for: " + ", ".join(missing), file=sys.stderr)
        print("Run ./run_all.sh first.", file=sys.stderr)
        return 2

    cases = {n: F.Case(n, man) for n in man}
    R = F.Results()
    for t in TESTS:
        try:
            t(R, cases)
        except Exception:
            traceback.print_exc()
            R.note(t.__name__[:2], t.__name__, "RAISED -- see the traceback above")
            R.check(t.__name__[:2], f"{t.__name__} completed", 1.0, 0.0, 0.5)

    buf = io.StringIO()
    with redirect_stdout(buf):
        R.print_table("G4beamline field maps, 2D cylinder and 3D grid -- "
                      "measured vs expected")
        ntot = sum(1 for r in R.rows if r["ok"] is not None)
        print(f"\nSUMMARY: {ntot - R.failed}/{ntot} checks passed, {R.failed} failed.")
    print(buf.getvalue())
    (HERE / "results.txt").write_text(buf.getvalue())

    try:
        import plot_tests
        plot_tests.make_all(man)
    except Exception:  # plotting must never mask a test failure
        traceback.print_exc()
    return 1 if R.failed else 0


if __name__ == "__main__":
    sys.exit(main())
