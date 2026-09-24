#!/usr/bin/env python3
"""Collect every number the report quotes, per case, into report/data.json.

Nothing here computes anything new: it reads the same files run_tests.py reads and
records both codes' values side by side, rather than only their difference.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np

from opalxruns.g4bl import read_map_header

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmplib
from cmplib import MUON_MASS, BG0

OUT = HERE.parent / "report"             # studies/g4bl/report
OUT.mkdir(exist_ok=True)
MAPS = cmplib.MAPS


def map_info(fname):
    header = read_map_header(MAPS / fname)
    if header is None:
        return dict(kind="unknown", file=fname)
    kind, d = header
    if kind == "grid":
        return dict(kind="grid", file=fname,
                    x=[d["X0"], d["X0"] + (d["nX"] - 1) * d["dX"]],
                    y=[d["Y0"], d["Y0"] + (d["nY"] - 1) * d["dY"]],
                    z=[d["Z0"], d["Z0"] + (d["nZ"] - 1) * d["dZ"]],
                    step=[d["dX"], d["dY"], d["dZ"]],
                    n=[int(d["nX"]), int(d["nY"]), int(d["nZ"])],
                    size_mb=round((MAPS / fname).stat().st_size / 1e6, 1))
    return dict(kind="cylinder", file=fname)


def plane(d, sub, z):
    dg = cmplib.read_bltrack(d / sub / f"Z{z}.txt")
    do = cmplib.read_monitor(d / sub / f"MON_{z}.h5")
    Xg, Xo, ids, miss = cmplib.matched(dg, do)
    return dg, do, Xg, Xo, ids, miss


def main():
    man = cmplib.load_manifest()
    out = {}
    for name, c in man.items():
        d = HERE / c["dir"]
        rec = dict(name=name, desc=c["desc"], bisector=bool(c.get("bisector")),
                   zstop=c["zstop"], monitors=c["monitors"],
                   dt_pair=c["dt_finer"] if c.get("converge") else c["dt_fine"],
                   maxstep_pair=c["maxstep_finer"] if c.get("converge") else c["maxstep_fine"],
                   dt_gauss=c["dt_fine"] if c.get("converge") else c["dt_coarse"],
                   maxstep_gauss=c["maxstep_fine"] if c.get("converge") else c["maxstep_coarse"],
                   maps=[], field={}, pair={}, gauss={}, convergence={})
        for m in c["maps"]:
            info = map_info(m["file"])
            info.update(scale=m["scale"], rotation=m.get("rot", ""))
            if c.get("bisector"):
                info["centreline_z_mm"] = m.get("cl_z")
            else:
                info["placed_at_m"] = [v / 1000.0 for v in m["xyz"]]
            rec["maps"].append(info)

        # ---- field, both codes -------------------------------------------
        for k in c["grids"]:
            fg = d / "pair" / f"g4bl_field_{k}.txt"
            fo = d / "pair" / f"opalx_field_{k}.dat"
            if not (fg.exists() and fo.exists()):
                continue
            pos, Bg, Bo = cmplib.match_field(cmplib.read_g4bl_field(fg),
                                             cmplib.read_opalx_field(fo))
            off = ~cmplib.on_face(pos, c["boxes"])
            diff = np.abs(Bo - Bg).max(axis=0)
            rec["field"][k] = dict(
                points=int(Bg.shape[1]), points_off_face=int(off.sum()),
                g4bl_peak_T=float(np.abs(Bg).max()), opalx_peak_T=float(np.abs(Bo).max()),
                worst_diff_T=float(diff[off].max()) if off.any() else None,
                median_diff_T=float(np.median(diff[off])) if off.any() else None,
                tol_T=cmplib.TOL["field_between"] * cmplib.printing_step(float(np.abs(Bg).max())))

        # ---- 19 particles, both codes ------------------------------------
        zi, zo = c["monitors"][0], c["monitors"][-1]
        sub = "fine" if (c.get("converge") and (d / "fine").is_dir()) else "pair"
        try:
            dg_i, do_i, Xg_i, Xo_i, ids_i, _ = plane(d, sub, zi)
            dg_o, do_o, Xg_o, Xo_o, ids_o, miss = plane(d, sub, zo)
            Mg = cmplib.transfer_matrix(Xg_i, Xg_o)
            Mo = cmplib.transfer_matrix(Xo_i, Xo_o)
            floors = cmplib.row_floors(dg_o, cmplib.STEP_SMALL)
            rec["pair"] = dict(
                n=len(ids_o), missing=miss, z_in=zi, z_out=zo,
                entrance_worst_pos=float(np.abs(Xo_i[[0, 2]] - Xg_i[[0, 2]]).max()),
                exit_worst_pos=float(np.abs(Xo_o[[0, 2]] - Xg_o[[0, 2]]).max()),
                exit_worst_ang=float(np.abs(Xo_o[[1, 3]] - Xg_o[[1, 3]]).max()),
                exit_worst_delta=float(np.abs(Xo_o[5] - Xg_o[5]).max()),
                orbit_max_x_g4bl=float(np.abs(Xg_o[0]).max()),
                p_change_g4bl=float(np.abs(cmplib.momentum_error(dg_o, "g4bl")
                                           - cmplib.momentum_error(dg_i, "g4bl")).max()),
                p_change_opalx=float(np.abs(cmplib.momentum_error(do_o, "opalx")
                                            - cmplib.momentum_error(do_i, "opalx")).max()),
                tof_g4bl=cmplib.time_of_flight(dg_o, "g4bl"),
                tof_opalx=cmplib.time_of_flight(do_o, "opalx"),
                M_g4bl=Mg.tolist(), M_opalx=Mo.tolist(), row_floors=floors.tolist(),
                matrix_worst_floors=float((np.abs(Mo - Mg) / floors[:, None]).max()),
                per_particle=[dict(
                    id=int(i),
                    g4bl=[float(v) for v in Xg_o[:, k]],
                    opalx=[float(v) for v in Xo_o[:, k]],
                ) for k, i in enumerate(ids_o)],
            )
        except Exception as e:
            rec["pair"] = dict(error=str(e))

        # ---- step convergence --------------------------------------------
        if c.get("converge") and (d / "fine").is_dir():
            try:
                rows = []
                for lbl, sd, dt, ms in (("coarse", "pair", c["dt_fine"], c["maxstep_fine"]),
                                        ("fine", "fine", c["dt_finer"], c["maxstep_finer"])):
                    _, _, Xg_a, Xo_a, _, _ = plane(d, sd, zi)
                    dgb, _, Xg_b, Xo_b, _, _ = plane(d, sd, zo)
                    fl = cmplib.row_floors(dgb, cmplib.STEP_SMALL)
                    rows.append(dict(
                        label=lbl, dt=dt, maxstep=ms,
                        pos=float(np.abs(Xo_b[[0, 2]] - Xg_b[[0, 2]]).max()),
                        ang=float(np.abs(Xo_b[[1, 3]] - Xg_b[[1, 3]]).max()),
                        matrix=float((np.abs(cmplib.transfer_matrix(Xo_a, Xo_b)
                                             - cmplib.transfer_matrix(Xg_a, Xg_b))
                                      / fl[:, None]).max())))
                rec["convergence"] = dict(steps=rows)
            except Exception as e:
                rec["convergence"] = dict(error=str(e))

        # ---- 20000 particles, both codes ---------------------------------
        if (d / "gauss").is_dir():
            planes = []
            for z in c["monitors"]:
                try:
                    dg, do, Xg, Xo, ids, miss = plane(d, "gauss", z)
                except Exception:
                    continue
                vp = np.abs(Xo[[0, 2]] - Xg[[0, 2]]).max(axis=0)
                va = np.abs(Xo[[1, 3]] - Xg[[1, 3]]).max(axis=0)
                planes.append(dict(
                    z=z, n=len(ids), missing=len(miss),
                    in_field=bool(cmplib.in_field(z, c["boxes"])),
                    mean_x=[float(Xg[0].mean()), float(Xo[0].mean())],
                    mean_y=[float(Xg[2].mean()), float(Xo[2].mean())],
                    rms_x=[float(Xg[0].std()), float(Xo[0].std())],
                    rms_y=[float(Xg[2].std()), float(Xo[2].std())],
                    med_pos=float(np.median(vp)), p99_pos=float(np.percentile(vp, 99)),
                    max_pos=float(vp.max()),
                    med_ang=float(np.median(va)), p99_ang=float(np.percentile(va, 99)),
                    max_ang=float(va.max())))
            rec["gauss"] = dict(n=c["n_gauss"], planes=planes)
        out[name] = rec
        print(f"  {name}")
    (OUT / "data.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT/'data.json'}")


if __name__ == "__main__":
    main()
