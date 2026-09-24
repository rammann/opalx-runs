#!/usr/bin/env python3
"""Compare every case's two runs and print one table with a pass/fail count.

Three stages per case, in the order in which a failure makes the next one
meaningless:

  field    the two codes are asked for the field at the same points, with no
           tracking at all. If this fails nothing below it means anything.
  pair     19 particles, compared one by one, plus a transfer matrix from each
           code by centred differences.
  gauss    20000 Gaussian muons, matched by id, so the comparison is still
           per particle and not only a comparison of distribution widths.

Run through run_all.sh, or on its own once both codes have been run:
    python run_tests.py [case ...]
"""

from __future__ import annotations

import io
import sys
import traceback
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmplib                                          # noqa: E402
from cmplib import TOL, Results                        # noqa: E402


def _plane_files(d: Path, sub: str, z: int):
    return d / sub / f"Z{z}.txt", d / sub / f"MON_{z}.h5"


def read_plane(d: Path, sub: str, z: int):
    """Both codes' record of one plane, lined up on particle id."""
    g4, op = _plane_files(d, sub, z)
    dg = cmplib.read_bltrack(g4)
    do = cmplib.read_monitor(op)
    return dg, do


# ---------------------------------------------------------------------------
# Stage A: the field, with no tracking
# ---------------------------------------------------------------------------
def test_field(R: Results, c: dict, d: Path) -> None:
    name = c["name"]
    boxes = c["boxes"]
    # The sample sets differ by case kind: a straight-frame case has grid/half/
    # yscan/edge, a bisector case has one lab-frame plane. Iterate what the case
    # actually declares rather than a fixed list, or a new kind is silently untested.
    tolmap = {"grid": TOL["field_grid_point"], "half": TOL["field_between"],
              "yscan": TOL["field_between"], "edge": None,
              "frame": TOL["field_between"]}
    for k in c["grids"]:
        tolmul = tolmap.get(k, TOL["field_between"])
        fg_p, fo_p = d / "pair" / f"g4bl_field_{k}.txt", d / "pair" / f"opalx_field_{k}.dat"
        if not (fg_p.exists() and fo_p.exists()):
            R.note(name, f"field {k}", "not run")
            continue
        fg = cmplib.read_g4bl_field(fg_p)
        fo = cmplib.read_opalx_field(fo_p)
        pos, Bg, Bo = cmplib.match_field(fg, fo)
        diff = np.abs(Bo - Bg).max(axis=0)
        face = cmplib.on_face(pos, boxes)
        peak = float(np.abs(Bg).max())
        # Three printing steps of this case's own peak field. See TOL in cmplib
        # for why the tolerance cannot be a fixed number of Tesla.
        tol = tolmul * cmplib.printing_step(peak) if tolmul else None

        if k == "edge":
            # This set sits on a box face on purpose. Report what each code
            # actually returned rather than asserting they agree: for a single
            # unrotated map both drop the face and it comes out zero, but where a
            # map is placed with a 180 degree rotation the face lands a few 1e-16 m
            # either side of itself and the two codes fall on opposite sides of it.
            nz_g = int((np.abs(Bg) > 0).any(axis=0).sum())
            nz_o = int((np.abs(Bo) > 0).any(axis=0).sum())
            R.check(name, "field on the top x face", float(diff.max()), 0.0, None,
                    note=f"{Bg.shape[1]} points, g4bl nonzero at {nz_g}, OPALX at "
                         f"{nz_o}; largest value returned {peak:.4g} T")
            continue

        off = ~face
        if off.any():
            R.check(name, f"field, {k} [T]", float(diff[off].max()), 0.0, tol,
                    note=f"peak {peak:.4g} T, {int(off.sum())} points clear of a face")
        if face.any():
            R.check(name, f"field, {k}, on a box face [T]", float(diff[face].max()),
                    0.0, None,
                    note=f"{int(face.sum())} points; the two codes round the edge "
                         f"differently, no particle is affected")

    if c.get("zero_field"):
        fo = cmplib.read_opalx_field(d / "pair" / "opalx_field_grid.dat")
        B = np.abs(np.vstack([fo["Bx"], fo["By"], fo["Bz"]])).max()
        R.check(name, "map of zeros gives zero", float(B), 0.0, TOL["field_zero"])


# ---------------------------------------------------------------------------
# Stage B: 19 particles, one by one
# ---------------------------------------------------------------------------
def test_pair(R: Results, c: dict, d: Path) -> None:
    name = c["name"]
    z_in, z_out = c["monitors"][0], c["monitors"][-1]
    try:
        dg_in, do_in = read_plane(d, "pair", z_in)
        dg_out, do_out = read_plane(d, "pair", z_out)
    except Exception as e:
        R.note(name, "pair stage", f"not run ({e})")
        return

    Xg_in, Xo_in, ids_in, miss_in = cmplib.matched(dg_in, do_in)
    Xg_out, Xo_out, ids_out, miss_out = cmplib.matched(dg_out, do_out)

    # Losing particles quietly would bias every number below, so the count is a
    # test and anything missing is named.
    for z, ids, miss, do in ((z_in, ids_in, miss_in, do_in),
                             (z_out, ids_out, miss_out, do_out)):
        dup = cmplib.duplicates(do)
        R.check(name, f"particles at z={z} mm", len(ids), c["n_pair"], 0,
                note=(f"missing from OPALX: {miss}" if miss else "")
                     + (f"  {len(dup)} crossed the plane twice" if dup else ""))
    if len(ids_in) < 2 or len(ids_out) < 2:
        return

    # Nothing has happened yet at the entrance plane. This is the check that the
    # two codes were handed the same beam, and it has to pass before any exit
    # number is worth reading.
    R.check(name, "entrance, worst |dx|,|dy| [m]",
            float(np.abs(Xo_in[[0, 2]] - Xg_in[[0, 2]]).max()), 0.0, TOL["entrance"])

    R.check(name, "exit, worst |dx|,|dy| [m]",
            float(np.abs(Xo_out[[0, 2]] - Xg_out[[0, 2]]).max()), 0.0, TOL["exit_pos"],
            note=f"orbit reaches |x| = {np.abs(Xg_out[0]).max():.4g} m")
    R.check(name, "exit, worst |dx'|,|dy'| [rad]",
            float(np.abs(Xo_out[[1, 3]] - Xg_out[[1, 3]]).max()), 0.0, TOL["exit_angle"])
    R.check(name, "exit, worst |d delta|",
            float(np.abs(Xo_out[5] - Xg_out[5]).max()), 0.0, TOL["exit_angle"])

    # A static magnetic field does no work. Checked inside each code separately,
    # so a broken integrator shows up even if both codes broke the same way.
    for code, din, dout in (("g4bl", dg_in, dg_out), ("opalx", do_in, do_out)):
        common = np.intersect1d(din["id"], dout["id"])
        i0 = np.searchsorted(din["id"], common)
        i1 = np.searchsorted(dout["id"], common)
        p0 = cmplib.momentum_error({k: v[i0] for k, v in din.items()}, code)
        p1 = cmplib.momentum_error({k: v[i1] for k, v in dout.items()}, code)
        R.check(name, f"|p| change through {code}", float(np.abs(p1 - p0).max()), 0.0,
                TOL["mom_conserved"])

    # Time of flight is deliberately not part of the six coordinates -- they are
    # measured from particle 0 of the same file -- so it is compared on its own.
    R.check(name, "time of flight, particle 0 [s]",
            cmplib.time_of_flight(do_out, "opalx"),
            cmplib.time_of_flight(dg_out, "g4bl"), None,
            note="diagnostic: g4bl prints 6 figures")

    # The transfer matrix uses only the small steps. The large ones are where
    # the map stops being linear, which is the point of having them.
    if len(ids_in) == c["n_pair"] == 19:
        Mg = cmplib.transfer_matrix(Xg_in, Xg_out)
        Mo = cmplib.transfer_matrix(Xo_in, Xo_out)
        diff = np.abs(Mo - Mg)
        floors = cmplib.row_floors(dg_out, cmplib.STEP_SMALL)
        floor = float(floors.max())
        R.check(name, "matrix floor, worst row", floor, 0.0, None,
                note="per row: " + " ".join(f"{f:.1e}" for f in floors))
        # The real test: every element, absolutely, in units of the floor.
        scaled = diff / floors[:, None]
        i = np.unravel_index(int(scaled.argmax()), scaled.shape)
        # On a case with a second run at half the step, this number is step limited
        # and test_convergence tests the extrapolated one instead. Reporting it as a
        # failure here would fail the case for using an affordable step.
        converging = bool(c.get("converge")) and (d / "fine").is_dir()
        R.check(name, "matrix, worst difference / floor", float(scaled.max()), 0.0,
                None if converging else TOL["matrix_floors"],
                note=f"worst |difference| {diff[i]:.3g} at "
                     f"M[{cmplib.COORD_NAMES[i[0]]} <- {cmplib.COORD_NAMES[i[1]]}]"
                     + ("; step limited, see the extrapolated row" if converging else ""))
        # Second check, on elements far enough above the floor that a ratio means
        # something. Everything else is either zero by symmetry -- a static field
        # cannot care what time a particle arrives, so zeta cannot move x -- or
        # small enough that the ratio would just report the floor.
        sig = cmplib.significant(Mg, floor, TOL["matrix_big"])
        if sig.any():
            R.check(name, "matrix, worst relative, large elements",
                    float((diff[sig] / np.abs(Mg)[sig]).max()), 0.0, TOL["matrix_rel"],
                    note=f"{int(sig.sum())} of 36 elements above "
                         f"{TOL['matrix_big']:.0f}x the floor")
        R.check(name, "symplectic residual, g4bl", cmplib.symplectic_residual(Mg), 0.0, None,
                note="diagnostic: finite differences, not the code")
        R.check(name, "symplectic residual, opalx", cmplib.symplectic_residual(Mo), 0.0, None)


# ---------------------------------------------------------------------------
# Step convergence: is what is left the codes, or the integrator?
# ---------------------------------------------------------------------------
def _extrapolate(dh: float, df: float, ratio: float):
    """What is left once the step size is taken out.

    If the error shrinks in proportion to the step -- halving the step halves it,
    a ratio near 2 -- then d(h) = d0 + a*h and d0 = 2*d(h/2) - d(h) removes the a*h
    part. That only works when the ratio really is about 2. When the error falls
    away faster than that the formula over-corrects and can even come out negative,
    which is meaningless. In that case the value at the finer step is already an
    upper bound on what is left, so use it and say which was used.
    """
    if 1.5 <= ratio <= 2.5:
        return abs(2.0 * df - dh), "removed"
    return df, "at the finest step"


def test_convergence(R: Results, c: dict, d: Path) -> None:
    """Run the same 19 particles at two step sizes and extrapolate to zero step.

    Needed because "converged" is a statement about the resolution you can measure
    at, not an absolute. The straight-frame cases cannot resolve below 1e-5 m, so
    DT = 1e-12 s and maxStep = 0.1 mm look converged there. On the bisector geometry
    the same quantity resolves to 1e-7 m, and at that resolution halving both steps
    halves the disagreement -- first order, so what is left is the integrator.

    Two numbers are reported. The ratio d(h) / d(h/2) says the residual really is
    step error: first order gives 2. The Richardson extrapolation
    d0 = 2*d(h/2) - d(h) is the disagreement with the step taken out, and THAT is
    what is tested against the printing floor. Testing d(h) instead would fail a
    case for using an affordable step.
    """
    name = c["name"]
    if not c.get("converge") or not (d / "fine").is_dir():
        return
    z = c["monitors"][-1]
    try:
        dg_h = cmplib.read_bltrack(d / "pair" / f"Z{z}.txt")
        do_h = cmplib.read_monitor(d / "pair" / f"MON_{z}.h5")
        dg_f = cmplib.read_bltrack(d / "fine" / f"Z{z}.txt")
        do_f = cmplib.read_monitor(d / "fine" / f"MON_{z}.h5")
    except Exception as e:
        R.note(name, "step convergence", f"not run ({e})")
        return

    def worst(dg, do, rows):
        Xg, Xo, _, _ = cmplib.matched(dg, do)
        return float(np.abs(Xo[rows] - Xg[rows]).max())

    floor = cmplib.printing_step(np.abs(dg_h["x"]).max() * 1e3) * 1e-3

    # The transfer matrix, extrapolated the same way. Built at both steps from the
    # same entrance plane, so only the exit state differs.
    try:
        zi = c["monitors"][0]
        Xg_i, Xo_i, _, _ = cmplib.matched(
            cmplib.read_bltrack(d / "pair" / f"Z{zi}.txt"),
            cmplib.read_monitor(d / "pair" / f"MON_{zi}.h5"))
        Xg_if, Xo_if, _, _ = cmplib.matched(
            cmplib.read_bltrack(d / "fine" / f"Z{zi}.txt"),
            cmplib.read_monitor(d / "fine" / f"MON_{zi}.h5"))
        Xg_h, Xo_h, _, _ = cmplib.matched(dg_h, do_h)
        Xg_f, Xo_f, _, _ = cmplib.matched(dg_f, do_f)
        floors = cmplib.row_floors(dg_h, cmplib.STEP_SMALL)
        mh = float((np.abs(cmplib.transfer_matrix(Xo_i, Xo_h)
                           - cmplib.transfer_matrix(Xg_i, Xg_h)) / floors[:, None]).max())
        mf = float((np.abs(cmplib.transfer_matrix(Xo_if, Xo_f)
                           - cmplib.transfer_matrix(Xg_if, Xg_f)) / floors[:, None]).max())
        ratio = mh / mf if mf else float("nan")
        R.check(name, "matrix step halving ratio", ratio, 2.0, None,
                note=f"{mh:.3g} floors at the step, {mf:.3g} at half it")
        val, how = _extrapolate(mh, mf, ratio)
        R.check(name, f"matrix, step error {how}", val, 0.0, TOL["matrix_floors"],
                note="in units of the printing floor")
    except Exception as e:
        R.note(name, "matrix convergence", f"not run ({e})")

    for lbl, rows, fl in (("|dx|,|dy| [m]", [0, 2], floor),
                          ("|dx'|,|dy'| [rad]", [1, 3], None)):
        dh = worst(dg_h, do_h, rows)
        df = worst(dg_f, do_f, rows)
        ratio = dh / df if df > 0 else float("nan")
        R.check(name, f"step halving ratio, {lbl}", ratio, 2.0, None,
                note=f"{dh:.3g} at the step, {df:.3g} at half it; 2 means the error "
                     f"shrinks in proportion to the step")
        val, how = _extrapolate(dh, df, ratio)
        R.check(name, f"step error {how}, {lbl}", val, 0.0,
                3.0 * fl if fl else TOL["exit_angle"],
                note=f"printing floor {fl:.1e} m" if fl else "")


# ---------------------------------------------------------------------------
# Stage C: 20000 Gaussian muons
# ---------------------------------------------------------------------------
def test_gauss(R: Results, c: dict, d: Path) -> None:
    name = c["name"]
    if not (d / "gauss").is_dir():
        R.note(name, "gauss stage", "not run")
        return
    free_pos, free_ang = [], []
    inf_pos = inf_ang = 0.0
    inf_planes = []
    stat = next((d / "gauss").glob("*.stat"), None)
    for z in c["monitors"]:
        try:
            dg, do = read_plane(d, "gauss", z)
        except Exception as e:
            R.note(name, f"gauss z={z} mm", f"unreadable ({e})")
            continue
        Xg, Xo, ids, miss = cmplib.matched(dg, do)
        # A momentum-dependent drop would bias every number below it, so the
        # count is tested at every plane, not only at the ends.
        R.check(name, f"gauss particles at z={z} mm", len(ids), c["n_gauss"], 0,
                note=(f"{len(miss)} missing from OPALX, ids {miss[:5]}..." if miss else ""))
        if len(ids) < 100:
            continue
        # per particle, not yet reduced -- the distribution is the point
        vpos = np.abs(Xo[[0, 2]] - Xg[[0, 2]]).max(axis=0)
        vang = np.abs(Xo[[1, 3]] - Xg[[1, 3]]).max(axis=0)
        if cmplib.in_field(z, c["boxes"]):
            inf_pos = max(inf_pos, float(vpos.max()))
            inf_ang = max(inf_ang, float(vang.max()))
            inf_planes.append(z)
        else:
            free_pos.append((z, vpos))
            free_ang.append((z, vang))
        for j, lbl in ((0, "x"), (2, "y")):
            R.check(name, f"gauss z={z} mm, mean {lbl} [m]",
                    float(Xo[j].mean()), float(Xg[j].mean()), TOL["gauss_mean"])
            R.check(name, f"gauss z={z} mm, rms {lbl}",
                    float(Xo[j].std()), float(Xg[j].std()), TOL["gauss_rms"], rel=True)
    # Tested on the median and the 99th percentile, which do not depend on how
    # many particles were drawn, and reported on the maximum, which does.
    # Per plane, never pooled across planes. Pooling puts the entrance plane --
    # where the two codes agree by construction, at 5e-8 m, because nothing has
    # happened yet -- into the same sample as the exit plane, and with two free
    # planes that is half the crossings. The median of the pool is then the
    # entrance value and says nothing about the element. Each plane is reported
    # on its own and the worst plane is what has to pass.
    for lbl, vals, tmed, tp99 in (("|dx|,|dy| [m]", free_pos, TOL["gauss_median_pos"],
                                   TOL["gauss_p99_pos"]),
                                  ("|dx'|,|dy'|", free_ang, TOL["gauss_median_angle"],
                                   TOL["gauss_p99_angle"])):
        if not vals:
            continue
        for z, v in vals:
            R.check(name, f"gauss z={z} mm, median {lbl}", float(np.median(v)), 0.0, None,
                    note=f"{len(v)} particles at this plane")
        zmed = max(vals, key=lambda t: np.median(t[1]))
        zp99 = max(vals, key=lambda t: np.percentile(t[1], 99))
        R.check(name, f"gauss free space, worst plane median {lbl}",
                float(np.median(zmed[1])), 0.0, tmed, note=f"at z={zmed[0]} mm")
        R.check(name, f"gauss free space, worst plane 99th pct {lbl}",
                float(np.percentile(zp99[1], 99)), 0.0, tp99, note=f"at z={zp99[0]} mm")
        allv = np.concatenate([v for _, v in vals])
        R.check(name, f"gauss free space, worst single particle {lbl}",
                float(allv.max()), 0.0, None,
                note="the tail of the same distribution, not a separate effect")
    if inf_planes:
        # Inside the field the two codes report the momentum from whichever
        # tracking step the particle was on, so the angle carries a fraction of a
        # step of sampling offset. Recorded as the length that explains it, not
        # tested: it says where the monitor sampled, not what the field is.
        # The largest field the reference orbit sees at any of these planes, not
        # the middle one: in the bend cases the reference has already left the
        # map sideways by the middle plane and sees nothing there.
        b = max((cmplib.field_at(stat, z / 1000.0) for z in inf_planes),
                default=0.0) if stat else 0.0
        if b > 1e-3:
            off = cmplib.offset_explaining(inf_ang, b)
            note = (f"planes {inf_planes}; at {b:.4g} T that is "
                    f"{off * 1e3:.3f} mm of sampling offset")
        else:
            # A quadrupole is zero on its axis, so there is no reference-orbit
            # field to divide by. The particles that show the difference are off
            # axis, where the field is not zero, but the .stat only carries the
            # reference orbit.
            note = (f"planes {inf_planes}; the reference orbit sees {b:.2g} T here, "
                    f"so no offset can be quoted")
        R.check(name, "gauss, |dx'|,|dy'| inside the field", inf_ang, 0.0, None,
                note=note)
        R.check(name, "gauss, |dx|,|dy| inside the field [m]", inf_pos, 0.0,
                TOL["exit_pos"], note="position IS interpolated onto the plane")


# ---------------------------------------------------------------------------
# A check between two cases rather than between the two codes
# ---------------------------------------------------------------------------
def test_zero_map_changes_nothing(R: Results, man: dict, names) -> None:
    """asr62_d3_group is asr62_dipole plus two placements of a map of zeros.

    asr62shim_280_sm_track.g4blmap is identically zero in all 52726 rows, and
    the second placement only turns it by Y180,Z180. So the two cases must give
    the same exit state, in each code on its own -- if they do not, either the
    file is not all zeros or placing an element changes something it should not.
    Run per code rather than between them, so it stays independent of everything
    else in this file.
    """
    a, b = "asr62_dipole", "asr62_d3_group"
    if not {a, b} <= set(names):
        return
    z = man[a]["monitors"][-1]
    if man[b]["monitors"][-1] != z:
        R.note(b, "zero map changes nothing", "the two cases do not share an exit plane")
        return
    for code, reader, fname in (("g4bl", cmplib.read_bltrack, f"Z{z}.txt"),
                                ("opalx", cmplib.read_monitor, f"MON_{z}.h5")):
        try:
            da = reader(HERE / man[a]["dir"] / "pair" / fname)
            db = reader(HERE / man[b]["dir"] / "pair" / fname)
        except Exception as e:
            R.note(b, f"zero map, {code}", f"not run ({e})")
            continue
        Xa = cmplib.canonical(da, code)
        Xb = cmplib.canonical(db, code)
        n = min(Xa.shape[1], Xb.shape[1])
        R.check(b, f"adding a map of zeros, {code} [m]",
                float(np.abs(Xb[[0, 2], :n] - Xa[[0, 2], :n]).max()), 0.0,
                TOL["entrance"], note="same exit as asr62_dipole, in this code alone")


def main(argv) -> int:
    man = cmplib.load_manifest()
    names = argv[1:] or list(man)
    R = Results()
    for n in names:
        c = man[n]
        d = HERE / c["dir"]
        for stage in (test_field, test_pair, test_convergence, test_gauss):
            try:
                stage(R, c, d)
            except Exception as e:
                R.note(n, stage.__name__, f"raised {type(e).__name__}: {e}")
                traceback.print_exc()
    try:
        test_zero_map_changes_nothing(R, man, names)
    except Exception as e:
        R.note("asr62_d3_group", "zero map check", f"raised {type(e).__name__}: {e}")

    buf = io.StringIO()
    with redirect_stdout(buf):
        R.print_table("g4bl_compare -- muE4 field maps, the same file through both codes")
        ntot = sum(1 for r in R.rows if r["ok"] is not None)
        print(f"\nSUMMARY: {ntot - R.failed}/{ntot} checks passed, {R.failed} failed.")
    print(buf.getvalue())
    (HERE / "results.txt").write_text(buf.getvalue())

    try:
        import plot_tests
        plot_tests.make_all(man, names)
    except Exception:
        # Plotting must never be able to hide a test failure.
        traceback.print_exc()
    return 1 if R.failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
