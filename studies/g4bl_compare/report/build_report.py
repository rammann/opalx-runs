#!/usr/bin/env python3
"""Build report/index.html from report/data.json and report/figs/."""
from __future__ import annotations
import json, html
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = json.loads((HERE / "data.json").read_text())
FULL = json.loads((HERE / "full_data.json").read_text()) if (HERE / "full_data.json").exists() else []
BZ = json.loads((HERE / "bz_data.json").read_text()) if (HERE / "bz_data.json").exists() else None
SCAN = json.loads((HERE / "scan_data.json").read_text()) if (HERE / "scan_data.json").exists() else None
SPIN = json.loads((HERE / "spin_data.json").read_text()) if (HERE / "spin_data.json").exists() else []

ORDER = ["asr61_dipole", "asr61_300sm", "asr62_dipole", "qsm600_quad",
         "asr61_group", "asr62_d2_group", "asr62_d3_group",
         "asr61_bisector", "asr61_group_bisector"]

TITLES = {
    "asr61_dipole": "ASR61_300d — the first bend's main map",
    "asr61_300sm": "ASR61_300sm — the small map beside it",
    "asr62_dipole": "ASR62 — the second and third bends' main map",
    "qsm600_quad": "QSM600 — the quadrupole",
    "asr61_group": "ASR61 group — main map plus both small ones",
    "asr62_d2_group": "ASR62 group, bend at 7504 mm — turned by 180°",
    "asr62_d3_group": "ASR62 group, bend at 12032 mm — not turned",
    "asr61_bisector": "ASR61_300d at the real muE4 angle",
    "asr61_group_bisector": "ASR61 group at the real muE4 angle",
}
WHY = {
    "asr61_dipole": "The map that bends the beam 40 degrees on the first turn. Sent in "
        "along the map's own axis, so it sweeps right across the map and leaves through "
        "the side. That is not how muE4 uses it, but it crosses more of the map than the "
        "real path does.",
    "asr61_300sm": "A small map covering x = 380 to 630 mm, outside the main map's reach. "
        "A beam on the axis would never touch it, so for this case alone it is moved "
        "sideways until its box sits on the axis. Both codes get the identical move.",
    "asr62_dipole": "The map used by both the second and third bends, at the polarity of "
        "the third. Same straight-on setup as ASR61_300d.",
    "qsm600_quad": "The only map file with nine columns, so it is the only one that "
        "exercises reading an electric field (which is zero everywhere in this file). Run "
        "at the strongest setting muE4 uses for this magnet.",
    "asr61_group": "All three ASR61 placements together, at the spacing muE4 gives them. "
        "The beam crosses from the main map into the small one, which is the hardest thing "
        "in the study for either code to get right.",
    "asr62_d2_group": "The three ASR62 placements of the middle bend. Both small files here "
        "are zero everywhere, so what this really tests is that OPALX handles a map turned "
        "by 180 degrees the same way G4beamline does.",
    "asr62_d3_group": "The same three files with no turn. Its exit must match ASR62 on its "
        "own, because the only thing added is a map of zeros — and it does.",
    "asr61_bisector": "The same main map, but placed the way muE4 actually places it: "
        "between two commands that turn the reference line 20 degrees each, so the beam "
        "enters the map at 20 degrees rather than straight on. This is the only case that "
        "checks OPALX and G4beamline agree on where the magnet is.",
    "asr61_group_bisector": "All three ASR61 placements at the real muE4 angle. The beam "
        "reaches 433 mm across the map, so it does pass through the small maps, and they "
        "move the exit by 7.7 mm.",
}

def esc(s): return html.escape(str(s))

def fmt(v, unit="", sig=4):
    if v is None: return "—"
    if isinstance(v, str): return esc(v)
    a = abs(v)
    if a == 0: return "0"
    if unit == "m":
        if a < 1e-9: return f"{v*1e12:.{sig-1}g} pm"
        if a < 1e-6: return f"{v*1e9:.{sig-1}g} nm"
        if a < 1e-3: return f"{v*1e6:.{sig-1}g} µm"
        return f"{v*1e3:.{sig-1}g} mm"
    if unit == "rad":
        if a < 1e-6: return f"{v*1e9:.{sig-1}g} nrad"
        if a < 1e-3: return f"{v*1e6:.{sig-1}g} µrad"
        return f"{v*1e3:.{sig-1}g} mrad"
    if unit == "T":
        if a < 1e-6: return f"{v*1e9:.{sig-1}g} nT"
        if a < 1e-3: return f"{v*1e6:.{sig-1}g} µT"
        return f"{v:.{sig}g} T"
    if unit == "s":
        return f"{v*1e9:.6f} ns"
    return f"{v:.{sig}g}{(' ' + unit) if unit else ''}"

def row(*cells, head=False):
    t = "th" if head else "td"
    return "<tr>" + "".join(f"<{t}>{c}</{t}>" for c in cells) + "</tr>"

def table(headers, rows, cls="data"):
    h = "<tr>" + "".join(f"<th>{c}</th>" for c in headers) + "</tr>"
    return (f'<div class="scroll"><table class="{cls}"><thead>{h}</thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


def setup_block(c):
    out = []
    for m in c["maps"]:
        if m["kind"] != "grid":
            out.append(row(esc(m["file"]), "cylinder map", "—", "—", fmt(m["scale"])))
            continue
        ext = (f'x {m["x"][0]:.0f}…{m["x"][1]:.0f}, y {m["y"][0]:.0f}…{m["y"][1]:.0f}, '
               f'z {m["z"][0]:.0f}…{m["z"][1]:.0f} mm')
        step = f'{m["step"][0]:.0f} × {m["step"][1]:.0f} × {m["step"][2]:.0f} mm'
        where = (f'centreline z = {m["centreline_z_mm"]:.4f} mm'
                 if c["bisector"] else
                 f'x = {m["placed_at_m"][0]:.3f}, z = {m["placed_at_m"][2]:.3f} m')
        if m["rotation"]:
            where += f' , turned {esc(m["rotation"])}'
        out.append(row(f'<code>{esc(m["file"])}</code>', ext, step,
                       f'{m["n"][0]}×{m["n"][1]}×{m["n"][2]}, {m["size_mb"]} MB',
                       f'{m["scale"]:+.5f}', where))
    return table(["field map file", "extent in its own frame", "grid step",
                  "points, file size", "strength", "placed at"], out)


def field_block(c):
    names = {"grid": "on the map's own grid points",
             "half": "halfway between grid points, in all three directions",
             "yscan": "lines across the full height",
             "edge": "exactly on the top x face (reported, not tested)",
             "frame": "a plane through the magnet, in lab coordinates"}
    rows = []
    for k, f in c["field"].items():
        rows.append(row(names.get(k, k), f'{f["points"]:,}',
                        fmt(f["g4bl_peak_T"], "T"), fmt(f["opalx_peak_T"], "T"),
                        fmt(f["worst_diff_T"], "T"), fmt(f["median_diff_T"], "T"),
                        fmt(f["tol_T"], "T") if k != "edge" else "—"))
    return table(["where the field was sampled", "points", "largest |B|, G4beamline",
                  "largest |B|, OPALX", "worst difference", "typical difference",
                  "allowed"], rows)


def pair_block(c):
    p = c["pair"]
    if "error" in p:
        return f'<p class="note">not available: {esc(p["error"])}</p>'
    r = [
        row("particles recorded", f'{p["n"]} of 19', f'{p["n"]} of 19',
            "none missing" if not p["missing"] else esc(str(p["missing"]))),
        row("difference at the entrance plane", "—", "—", fmt(p["entrance_worst_pos"], "m")),
        row("largest |x| reached at the exit", fmt(p["orbit_max_x_g4bl"], "m"),
            fmt(p["orbit_max_x_g4bl"], "m"), "—"),
        row("worst exit position difference", "—", "—", fmt(p["exit_worst_pos"], "m")),
        row("worst exit angle difference", "—", "—", fmt(p["exit_worst_ang"], "rad")),
        row("change in |p| through the magnet", fmt(p["p_change_g4bl"]),
            fmt(p["p_change_opalx"]), "a magnet cannot change |p|"),
        row("time of flight, first particle", fmt(p["tof_g4bl"], "s"),
            fmt(p["tof_opalx"], "s"), fmt(p["tof_opalx"] - p["tof_g4bl"], "s")),
        row("transfer matrix, worst cell", "—", "—",
            f'{p["matrix_worst_floors"]:.2f} printing steps'),
    ]
    return table(["quantity", "G4beamline", "OPALX", "difference"], r)


def matrices(c):
    p = c["pair"]
    if "error" in p: return ""
    names = ["x", "x'", "y", "y'", "ζ", "δ"]
    def m2t(M, title):
        rows = [row(names[i], *[f'{M[i][j]:+.4f}' for j in range(6)])
                for i in range(6)]
        return (f'<h4>{title}</h4>'
                + table(["out ↓ / in →"] + names, rows, cls="data mono small"))
    return ('<div class="pair-cols">' + m2t(p["M_g4bl"], "Transfer matrix, G4beamline")
            + m2t(p["M_opalx"], "Transfer matrix, OPALX") + '</div>')


def per_particle(c):
    p = c["pair"]
    if "error" in p: return ""
    rows = []
    for e in p["per_particle"]:
        g, o = e["g4bl"], e["opalx"]
        rows.append(row(e["id"],
                        f'{g[0]*1e3:.4f}', f'{o[0]*1e3:.4f}', fmt(o[0]-g[0], "m", 3),
                        f'{g[1]*1e3:.4f}', f'{o[1]*1e3:.4f}', fmt(o[1]-g[1], "rad", 3),
                        fmt(o[2]-g[2], "m", 3), fmt(o[3]-g[3], "rad", 3)))
    return table(["particle", "x G4BL [mm]", "x OPALX [mm]", "Δx",
                  "x' G4BL [mrad]", "x' OPALX [mrad]", "Δx'", "Δy", "Δy'"], rows)


def gauss_block(c):
    g = c.get("gauss")
    if not g: return '<p class="note">not run</p>'
    rows = []
    for pl in g["planes"]:
        tag = ' <span class="chip">inside the field</span>' if pl["in_field"] else ""
        rows.append(row(f'{pl["z"]}{tag}', f'{pl["n"]:,}',
                        f'{pl["mean_x"][0]*1e3:.4f}', f'{pl["mean_x"][1]*1e3:.4f}',
                        fmt(pl["mean_x"][1]-pl["mean_x"][0], "m", 3),
                        f'{pl["rms_x"][0]*1e3:.4f}', f'{pl["rms_x"][1]*1e3:.4f}',
                        fmt(pl["rms_x"][1]-pl["rms_x"][0], "m", 3),
                        fmt(pl["med_pos"], "m", 3), fmt(pl["p99_pos"], "m", 3),
                        fmt(pl["max_pos"], "m", 3)))
    return table(["plane z [mm]", "particles", "mean x G4BL [mm]", "mean x OPALX [mm]",
                  "Δ mean x", "rms x G4BL [mm]", "rms x OPALX [mm]", "Δ rms x",
                  "median per-particle Δ", "99th percentile", "worst"], rows)


def conv_block(c):
    cv = c.get("convergence") or {}
    if not cv.get("steps"): return ""
    rows = [row(f'{s["dt"]} s / {s["maxstep"]} mm', fmt(s["pos"], "m"),
                fmt(s["ang"], "rad"), f'{s["matrix"]:.2f}') for s in cv["steps"]]
    return ('<h4>What happens when the step size is halved</h4>'
            + table(["step size", "worst position difference",
                     "worst angle difference", "matrix, printing steps"], rows))


def section(name):
    c = D[name]
    figs = [(f"figs/{name}_field.png", "The field both codes report, and where they differ"),
            (f"figs/{name}_pair.png", "Nineteen particles, one by one, and the transfer matrix"),
            (f"figs/{name}_gauss.png", "Twenty thousand particles: beam size, beam centre, "
                                       "and how the difference depends on amplitude")]
    figs = [(p, cap) for p, cap in figs if (HERE / p).exists()]
    run = (f'Step size: {c["dt_pair"]} s in OPALX and {c["maxstep_pair"]} mm in G4beamline '
           f'for the 19 particles; {c["dt_gauss"]} s and {c["maxstep_gauss"]} mm for the '
           f'20 000. Tracking stops at {c["zstop"]:.2f} m of path length. '
           f'Recording planes at {", ".join(str(z) for z in c["monitors"])} mm.')
    frame = ("Placed the way muE4 places it, between two <code>cornerarc</code> commands "
             "that each turn the reference line 20 degrees."
             if c["bisector"] else
             "Reference line kept straight, so G4beamline's z in mm divided by 1000 equals "
             "OPALX's path length in m until the field starts to bend the beam.")
    return f"""
<section id="{esc(name)}">
  <h2>{esc(TITLES[name])}</h2>
  <p class="lede">{esc(WHY[name])}</p>
  <h3>Setup</h3>
  {setup_block(c)}
  <p class="note">{frame} {esc(run)}</p>
  <h3>The field, with no tracking</h3>
  <p>Both codes are asked for the field at the same points. No particles are involved, so a
  difference here is a difference in how the file is read or where the magnet is placed.</p>
  {field_block(c)}
  {f'<figure><img src="{figs[0][0]}" alt="{esc(figs[0][1])}"><figcaption>{esc(figs[0][1])}</figcaption></figure>' if figs else ''}
  <h3>Nineteen particles</h3>
  {pair_block(c)}
  {matrices(c)}
  {conv_block(c)}
  {f'<figure><img src="{figs[1][0]}" alt="{esc(figs[1][1])}"><figcaption>{esc(figs[1][1])}</figcaption></figure>' if len(figs) > 1 else ''}
  <details><summary>Every particle, one row each</summary>{per_particle(c)}</details>
  <h3>Twenty thousand particles</h3>
  {gauss_block(c)}
  {f'<figure><img src="{figs[2][0]}" alt="{esc(figs[2][1])}"><figcaption>{esc(figs[2][1])}</figcaption></figure>' if len(figs) > 2 else ''}
</section>"""



def full_line():
    """The whole muE4 line, both codes, 22 recording planes."""
    if not FULL:
        return ""
    r0, rl = FULL[0], FULL[-1]
    # per-plane numbers, everything
    a = [row(f'{r["z"]}', f'{r["ng"]:,}', f'{r["no"]:,}',
             f'{r["mean_x"][0]*1e3:.4f}', f'{r["mean_x"][1]*1e3:.4f}',
             f'{r["rms_x"][0]*1e3:.4f}', f'{r["rms_x"][1]*1e3:.4f}',
             f'{r["rms_y"][0]*1e3:.4f}', f'{r["rms_y"][1]*1e3:.4f}',
             fmt(r["med_pos"], "m", 3), fmt(r["p99_pos"], "m", 3),
             fmt(r["max_pos"], "m", 3)) for r in FULL]
    t_all = table(["plane z [mm]", "particles G4BL", "particles OPALX",
                   "mean x G4BL [mm]", "mean x OPALX [mm]",
                   "rms x G4BL [mm]", "rms x OPALX [mm]",
                   "rms y G4BL [mm]", "rms y OPALX [mm]",
                   "median Δ per particle", "99th percentile", "worst"], a)
    # per-plane numbers, only what the real collimators would pass
    b = [row(f'{r["z"]}', f'{r["kept"]:,}',
             f'{r["cut_mean_x"][0]*1e3:.4f}', f'{r["cut_mean_x"][1]*1e3:.4f}',
             f'{r["cut_rms_x"][0]*1e3:.4f}', f'{r["cut_rms_x"][1]*1e3:.4f}',
             f'{r["cut_rms_y"][0]*1e3:.4f}', f'{r["cut_rms_y"][1]*1e3:.4f}',
             fmt(r["cut_max_pos"], "m", 3)) for r in FULL]
    t_cut = table(["plane z [mm]", "particles kept",
                   "mean x G4BL [mm]", "mean x OPALX [mm]",
                   "rms x G4BL [mm]", "rms x OPALX [mm]",
                   "rms y G4BL [mm]", "rms y OPALX [mm]", "worst Δ"], b)
    planes = "".join(
        f'<figure><img src="figs/full/plane_{r["z"]:05d}.png" loading="lazy" '
        f'alt="Phase space at centreline z = {r["z"]} mm, OPALX above G4beamline">'
        f'<figcaption>Centreline z = {r["z"]} mm. OPALX above, G4beamline below, same axes '
        f'down each column. {r["n"]:,} particles matched one to one.</figcaption></figure>'
        for r in FULL)
    dy_raw = (rl["rms_y"][1] - rl["rms_y"][0]) * 1e6
    dy_cut = (rl["cut_rms_y"][1] - rl["cut_rms_y"][0]) * 1e6
    return f"""
<section id="fullline">
  <h2>The whole muE4 line</h2>
  <p class="lede">All 22 field maps at the positions the <code>cornerarc</code> walk gives
  them — the same positions as the existing <code>lattice.in</code>, digit for digit — with
  22 recording planes along the 19.4 m. 10 000 muons, one file read by both codes. Nothing
  scrapes.</p>

  <h3>Setup</h3>
  {table(["", "value"], [
     row("field maps", "22, at their real muE4 positions and rotations"),
     row("recording planes", "22, roughly every metre, only where the reference line is "
         "straight; the six turns occupy 2578-3387, 7142-7866 and 11670-12394 mm and a "
         "plane inside a turn has no well defined frame"),
     row("particles", "10 000 Gaussian muons, 10 mm by 3 mm, 2 mrad, 28 MeV/c, sample mean "
         "subtracted exactly"),
     row("collimators", "none. The 12 iron objects are vacuum in G4beamline and the 16 "
         "COLLIMATOR elements are left out of the OPALX line"),
     row("step size", "1e-12 s in OPALX, 0.1 mm in G4beamline"),
     row("bending", "108 degrees over three magnets"),
     row("run time", "G4beamline 2970 s, OPALX 534 s, one core each"),
  ])}

  <h3>What came out</h3>
  <p>Every particle reached every plane in both codes, except the last plane where OPALX
  recorded {rl['no']:,} of {rl['ng']:,} — the monitor drop already known from the earlier
  whole-line work. Through the first 12.5 m <strong>not one particle in 10 000 differs by
  more than 100&nbsp;µm</strong>, and the median difference stays between 1 and 5&nbsp;µm
  the whole way.</p>
  {t_all}
  <figure><img src="figs/full/growth.png" alt="How the difference between the codes grows
  along the line"><figcaption>The median stays near a few µm for the whole line. The worst
  single particle climbs steeply after 13 m — that is the halo, explained below, not the
  beam.</figcaption></figure>
  <figure><img src="figs/full/along_line.png" alt="Beam size, beam centre and emittance
  along the line, both codes"><figcaption>Beam size, beam centre and emittance against
  position. G4beamline drawn as a line, OPALX as points on top of it.</figcaption></figure>

  <h3>A halo the real line would never have</h3>
  <p>Taken at face value the beam size in y at the last plane is
  {abs(rl['rms_y'][1]-rl['rms_y'][0])/rl['rms_y'][0]*100:.2f}&nbsp;% apart
  ({dy_raw:+.0f}&nbsp;µm), which looks bad. It is not the beam. About 20 particles out of
  10 000 fly to y = ±300&nbsp;mm, far outside every field map — the ASR61 maps reach
  ±120&nbsp;mm vertically and the ASR62 maps ±110&nbsp;mm. Out there neither code has any
  field, and whether a particle was inside or outside at the boundary differs, so those few
  part company completely; the worst is 645&nbsp;mm apart.</p>
  <p>Those particles exist only because the collimators were opened for this run. The real
  line scrapes them long before. Counting only what fits through the tightest collimator
  muE4 has — 190&nbsp;mm across, 105&nbsp;mm high — the same number becomes
  <strong>{dy_cut:+.1f}&nbsp;µm</strong>, and the beam size in x is
  {abs(rl['cut_rms_x'][1]-rl['cut_rms_x'][0])/rl['cut_rms_x'][0]*100:.3f}&nbsp;% apart.</p>
  <p class="note">Opening the collimators isolated the optics as intended, but it produced
  a halo with no physical meaning that then dominated the beam size. Both sets of numbers
  are given so it is clear which is which.</p>
  {t_cut}

  <h3>Phase space at every plane</h3>
  <p>OPALX above, G4beamline below, with the same axes down each column so the two rows can
  be read against each other. Axis limits come from percentiles rather than the extremes —
  scaling to the few halo particles would squash the beam itself into a dot — and each
  figure says how many points fall outside the view.</p>
  {planes}

  <h3>What is still open</h3>
  <p>The step size has <strong>not</strong> been shown to be small enough on the full line.
  The single magnet cases showed that only halving the step and looking can answer that, and
  that has not been done here. The worst-case differences of 100 to 200&nbsp;µm in the
  middle of the line could be step size or could be real.</p>
</section>"""




def spin_section():
    """Thomas-BMT spin precession: uniform fields, then the whole muE4 line."""
    if not (BZ and SCAN):
        return ""
    f = SCAN["fitted"]
    bz_rows = [row("G4beamline", f'{BZ["g4bl"]["polx"]:.9f}', f'{BZ["g4bl"]["poly"]:.9f}',
                   f'{BZ["g4bl"]["angle"]:.9f}', f'{BZ["g4bl"]["angle"]-BZ["theory"]:+.3e}'),
               row("OPALX", f'{BZ["opalx"]["polx"]:.9f}', f'{BZ["opalx"]["poly"]:.9f}',
                   f'{BZ["opalx"]["angle"]:.9f}', f'{BZ["opalx"]["angle"]-BZ["theory"]:+.3e}'),
               row("closed form", "&mdash;", "&mdash;", f'{BZ["theory"]:.9f}', "&mdash;")]
    step_rows = [row(s["dt"], f'{s["mm"]:.5f}', f'{s["angle"]:.9f}',
                     f'{s["err"]:+.3e}', f'{s["bound"]:.3e}') for s in BZ["steps"]]
    scan_rows = [row(f'{q["p"]:.3f}', f'{q["gamma"]:.4f}', f'{q["theta"]:.7f}',
                     f'{q["g4bl"]*1e3:.5f}', f'{q["opalx"]*1e3:.5f}',
                     f'{q["theory"]*1e3:.5f}', f'{(q["opalx"]-q["g4bl"])*1e3:+.6f}')
                 for q in SCAN["points"]]
    tm = BZ["two_metre"]
    e0 = abs(BZ["steps"][0]["err"])

    full = ""
    if SPIN:
        rl = SPIN[-1]
        fr = [row(f'{r["z"]}', f'{np.degrees(r["theta"]):.2f}&deg;', f'{r["n"]:,}',
                  f'{r["med_raw"]:.6f}', f'{r["med_rot"]:.6f}', f'{r["p99_rot"]:.6f}',
                  f'{r["dg"]*1e3:.5f}', f'{r["do"]*1e3:.5f}') for r in SPIN]
        turns = sorted({round(np.degrees(r["theta"]), 1) for r in SPIN if r["theta"]})
        fcheck = "".join(
            f"<tr><td>plane turned by {t:.0f}&deg;</td>"
            f"<td>{np.median([r['med_raw'] for r in SPIN if round(np.degrees(r['theta']),1)==t]):.6f}</td>"
            f"<td>{2*np.sin(np.radians(t)/2):.6f}</td></tr>" for t in turns)
        full = f"""
  <h3>The whole muE4 line</h3>
  <p>The same 22 field maps and 22 recording planes as the section below, with spin on in
  both codes and the spin starting along the direction of travel.</p>
  <p><strong>OPALX reports the spin in a different frame from the momentum.</strong> A
  <code>MONITOR</code> records position and momentum in its own frame &mdash;
  <code>Monitor::apply</code> runs after the bunch has been transformed into that frame,
  which is why the plane-crossing test is on <code>R(2)</code> &mdash; but it takes the spin
  straight from <code>pc-&gt;Pol</code>, and <code>ParticleContainer::transformBunch</code>
  rotates <code>R</code>, <code>P</code>, <code>E</code> and <code>B</code> while leaving
  <code>Pol</code> alone. G4beamline rotates both into the centreline frame. For a plane
  that is not turned the two frames are the same, which is why none of the uniform-field
  cases above could show it; muE4's planes are turned by up to 108 degrees, so it shows
  there. Turning OPALX's spin into each plane's own frame by that plane's
  <code>THETA</code> makes the two codes agree.</p>
  <p>The difference is not merely large, it is exactly what a frame rotation gives.
  Turning a unit vector by an angle <em>t</em> moves it by 2&nbsp;sin(<em>t</em>/2):</p>
  {table(["", "median difference measured", "2 sin(t/2)"], [fcheck])}
  {table(["plane z [mm]", "how far the frame has turned", "particles",
          "median difference as written", "median after turning OPALX into the plane's frame",
          "99th percentile after turning", "spin ahead of momentum, G4BL [mrad]",
          "OPALX [mrad]"], fr)}
  <p>Turning OPALX's spin into each plane's own frame removes almost all of it — the median
  falls from 0.68 to between 0.001 and 0.012. What is left is a real difference, and it is
  the size the field edges predict: the uniform Bz case showed OPALX carrying up to half a
  step of extra spin rotation at each edge, this line has 22 maps and so 44 edges, and
  44 of those with random signs comes to roughly 2e-4&nbsp;rad.</p>
  <p class="note"><strong>This is not a precise measurement of the anomaly and was never
  going to be.</strong> The signal is under a milliradian and the accumulated edge error is
  a few tenths of one, so they are comparable. The anomaly is measured by the uniform-field
  scan above, to 8 parts per million. What this line shows is that spin survives 22 real
  field maps and 108&deg; of bending with nothing going wrong. Note also that the naive
  <em>G</em>&gamma;&theta; = 2.27&nbsp;mrad does not apply here: the WSX solenoid turns the
  spin about the direction of travel before the beam reaches any bend, so the spin is no
  longer along the momentum when the bending starts.</p>
  <p>One more thing this run shows. The length of the spin vector must stay at 1, and in
  OPALX it drifts to {SPIN[-1]["mago"]:.7f} over the line while G4beamline holds
  {SPIN[-1]["magg"]:.7f}. The rotation itself conserves the length exactly; the storage does
  not — the spin is kept as three <code>float</code>s rather than doubles
  (<code>ParticleContainer.hpp:124</code>), and about 250&nbsp;000 steps of rounding back to
  float accumulate to roughly that.</p>
  <figure><img src="figs/spin/spin_full_line.png" alt="Spin through the whole muE4 line">
  <figcaption>Left: the difference between the two codes' spin, before and after turning
  OPALX's into the plane's frame. Middle: the spin along the line. Right: the anomalous
  precession building up over 108 degrees of bending.</figcaption></figure>
"""

    return f"""
<section id="spin">
  <h2>Spin</h2>
  <p class="lede">Both codes integrate the Thomas-BMT equation and neither had ever been
  checked against the other. In OPALX it is switched on by one attribute,
  <code>POLARIZATION</code> on the <code>BEAM</code>, and no input file in the repository
  used it &mdash; only unit tests, which deliberately do not check which way the spin
  turns.</p>

  <h3>A uniform Bz, where the answer is a formula</h3>
  <p>A muon travelling along a uniform field feels no force at all, so the orbit is a
  straight line and the only thing that changes is the spin. It turns about the field by
  (1&nbsp;+&nbsp;<em>G</em>) times the cyclotron angle. One metre of 0.0933979467&nbsp;T is
  exactly one radian for a 28&nbsp;MeV/c muon.</p>
  {table(["", "spin x", "spin y", "angle [rad]", "minus the closed form"], bz_rows)}
  <p>Both turn the spin the same way, and that way is the one the torque on a magnetic
  moment gives for a positive muon: from +x towards &minus;y. This is the first check of
  the direction in OPALX. The momentum comes out unchanged in both, as it must.</p>

  <h3>Why OPALX sits further from the formula, and why that is not the equation</h3>
  <p>OPALX is a hundred times further from the closed form than G4beamline. Shrinking the
  time step shows what that is: every error sits inside half a step's worth of rotation,
  and that bound falls with the step, but the error jumps about inside it rather than
  falling smoothly.</p>
  {table(["DT [s]", "mm per step", "angle [rad]", "minus the closed form",
          "half a step's rotation"], step_rows)}
  <p>Doubling the length of field settles it. The rotation doubled and the error did not
  &mdash; it went from {e0:.2e} to {abs(tm['err']):.2e}&nbsp;rad. Nothing builds up along
  the field, so the integration itself is right and what is left is made at the field edge.
  Geant4 steps exactly to a geometry boundary; OPALX steps straight through one and carries
  up to half a step of extra rotation.</p>
  <figure><img src="figs/spin/bz_step.png" alt="How the spin error depends on the step and
  on the length of field"><figcaption>Left: the error is bounded by the step but not smooth
  in it. Right: twice the rotation, no more error, so nothing builds up along the
  field.</figcaption></figure>

  <h3>The anomaly itself</h3>
  <p>In a field across the direction of travel the spin turns
  (1&nbsp;+&nbsp;<em>G</em>&gamma;) times faster than the momentum, so after a bend of
  &theta; the spin leads the momentum by <em>G</em>&gamma;&theta;. That is the only term in
  the equation which is not the same one that turns the momentum. Holding the bend angle
  fixed and changing the momentum makes that angle a straight line in &gamma; whose slope
  is <em>G</em>, so the five cases together <em>measure</em> the anomaly instead of
  comparing two numbers once.</p>
  {table(["p [GeV/c]", "&gamma;", "bend [rad]", "G4beamline [mrad]", "OPALX [mrad]",
          "theory [mrad]", "OPALX &minus; G4BL [mrad]"], scan_rows)}
  {table(["", "anomaly recovered from the five points", "the value compiled in", "apart"], [
     row("OPALX", f'{f["opalx"]:.9e}', f'{f["opalx_compiled"]:.9e}',
         f'{(f["opalx"]-f["opalx_compiled"])/f["opalx_compiled"]*100:+.4f} %'),
     row("G4beamline", f'{f["g4bl"]:.9e}', f'{f["g4bl_compiled"]:.9e}',
         f'{(f["g4bl"]-f["g4bl_compiled"])/f["g4bl_compiled"]*100:+.4f} %')])}
  <p>OPALX recovers its own anomaly to 8 parts per million across a factor 28 in &gamma;.
  G4beamline is further off only because it prints six digits, which limits an angle to
  about 0.001&nbsp;mrad &mdash; and that is also where the two codes sit relative to each
  other.</p>
  <figure><img src="figs/spin/g2_scan.png" alt="The anomalous precession against gamma">
  <figcaption>Left: the anomalous precession is a straight line in &gamma; and both codes
  lie on it. Right: both sit inside the floor of what G4beamline can print.</figcaption>
  </figure>
{full}
  <h3>What this cannot cover</h3>
  <p><code>FROMFILE</code> ignores any spin in the particle file and stamps the one
  <code>BEAM</code> vector on every particle, so a beam whose particles start with
  different spins cannot be compared &mdash; which is exactly what the real muE4 surface
  muon file has. <code>POLARIZATION</code> is also rejected for anything but muons, even
  though the electron and proton anomalies are present in the source. Neither code applies
  the force from the magnetic moment: G4beamline documents that as a Geant4 bug and the
  OPALX pusher only turns the spin, so on that point they agree.</p>
</section>"""


def summary_rows():
    rows = []
    for n in ORDER:
        c = D[n]
        p, f = c["pair"], c["field"]
        worst_field = max((v["worst_diff_T"] or 0) for v in f.values()) if f else None
        g = c.get("gauss") or {}
        pl = [x for x in g.get("planes", []) if not x["in_field"]]
        med = max((x["med_pos"] for x in pl), default=None)
        rows.append(row(f'<a href="#{n}">{esc(n)}</a>',
                        f'{len(c["maps"])}',
                        "muE4 angle" if c["bisector"] else "straight",
                        fmt(worst_field, "T"),
                        fmt(p.get("exit_worst_pos"), "m"),
                        fmt(p.get("exit_worst_ang"), "rad"),
                        fmt(med, "m")))
    return table(["case", "maps", "geometry", "worst field difference",
                  "worst exit position Δ, 19 particles", "worst exit angle Δ",
                  "median Δ, 20 000 particles"], rows)


toc = ('<li><a href="#spin"><strong>Spin</strong></a></li>'
       + '<li><a href="#fullline"><strong>The whole line</strong></a></li>'
       + "".join(f'<li><a href="#{n}">{esc(TITLES[n].split(" — ")[0])}</a></li>' for n in ORDER))
body = "".join(section(n) for n in ORDER)

HTML = f"""<title>muE4 Field Map Comparison</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {{
  --ground:#f6f7f9; --surface:#ffffff; --raised:#eef1f5;
  --ink:#14171c; --ink-2:#4d5560; --ink-3:#767f8c;
  --line:#dce1e8; --line-2:#c6cdd7;
  --opalx:#1f5fd1; --g4bl:#c1432d;
  --ok:#1b6b45; --warn:#8a5a00; --bad:#a8291c;
  --sans:"IBM Plex Sans",system-ui,-apple-system,sans-serif;
  --serif:"IBM Plex Serif",Georgia,serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground:#101317; --surface:#171b21; --raised:#1e242c;
    --ink:#e9ecf1; --ink-2:#aeb7c3; --ink-3:#8b95a3;
    --line:#2a313b; --line-2:#3a434f;
    --opalx:#6d9bec; --g4bl:#e07b66;
    --ok:#4fbd8a; --warn:#d9a441; --bad:#e2776a;
  }}
}}
:root[data-theme="dark"] {{
  --ground:#101317; --surface:#171b21; --raised:#1e242c;
  --ink:#e9ecf1; --ink-2:#aeb7c3; --ink-3:#8b95a3;
  --line:#2a313b; --line-2:#3a434f;
  --opalx:#6d9bec; --g4bl:#e07b66;
  --ok:#4fbd8a; --warn:#d9a441; --bad:#e2776a;
}}
body {{ background:var(--ground); color:var(--ink); font-family:var(--sans);
  font-size:15px; line-height:1.6; margin:0; }}
.wrap {{ max-width:1180px; margin:0 auto; padding-inline:16px; padding-block:0 64px;
  display:grid; grid-template-columns:200px minmax(0,1fr); gap:40px; align-items:start; }}
@media (max-width:900px) {{ .wrap {{ grid-template-columns:minmax(0,1fr); gap:0; }}
  nav.toc {{ position:static !important; margin-block:8px 24px; }} }}
header.top {{ max-width:1180px; margin:0 auto; padding-inline:16px; padding-block:48px 28px; }}
h1 {{ font-family:var(--serif); font-size:34px; line-height:1.15; margin:0 0 10px;
  font-weight:600; text-wrap:balance; letter-spacing:-0.01em; }}
.sub {{ color:var(--ink-2); max-width:62ch; margin:0 0 18px; }}
.meta {{ font-family:var(--mono); font-size:12px; color:var(--ink-3);
  display:flex; flex-wrap:wrap; gap:6px 18px; }}
nav.toc {{ position:sticky; top:calc(env(safe-area-inset-top, 0px) + 16px);
  font-size:13px; border-left:2px solid var(--line); padding-left:14px; }}
nav.toc ol {{ list-style:none; margin:0; padding:0; display:flex;
  flex-direction:column; gap:7px; }}
nav.toc a {{ color:var(--ink-2); text-decoration:none; }}
nav.toc a:hover {{ color:var(--opalx); }}
nav.toc .lbl {{ font-family:var(--mono); font-size:10px; letter-spacing:.09em;
  text-transform:uppercase; color:var(--ink-3); margin-bottom:10px; display:block; }}
main {{ min-width:0; }}
section {{ padding-block:34px; border-top:1px solid var(--line); }}
section:first-of-type {{ border-top:none; }}
h2 {{ font-family:var(--serif); font-size:23px; font-weight:600; margin:0 0 8px;
  text-wrap:balance; }}
h3 {{ font-family:var(--mono); font-size:11px; letter-spacing:.1em; text-transform:uppercase;
  color:var(--ink-3); margin:30px 0 10px; font-weight:500; }}
h4 {{ font-size:13px; margin:18px 0 8px; color:var(--ink-2); font-weight:600; }}
p {{ max-width:70ch; }}
.lede {{ color:var(--ink-2); margin:0 0 4px; }}
.note {{ font-size:13px; color:var(--ink-3); border-left:2px solid var(--line-2);
  padding-left:12px; margin:12px 0; }}
.scroll {{ overflow-x:auto; margin:10px 0 4px; }}
table.data {{ border-collapse:collapse; font-size:12.5px; width:100%;
  font-variant-numeric:tabular-nums; }}
table.data th {{ text-align:left; font-family:var(--mono); font-size:10.5px; font-weight:500;
  letter-spacing:.04em; text-transform:uppercase; color:var(--ink-3);
  border-bottom:1px solid var(--line-2); padding:7px 12px 7px 0; white-space:nowrap; }}
table.data td {{ padding:6px 12px 6px 0; border-bottom:1px solid var(--line);
  white-space:nowrap; }}
table.data td:first-child {{ color:var(--ink-2); }}
table.mono td {{ font-family:var(--mono); }}
table.small {{ font-size:11.5px; }}
code {{ font-family:var(--mono); font-size:.9em; background:var(--raised);
  padding:1px 5px; border-radius:3px; }}
figure {{ margin:18px 0 6px; background:var(--surface); border:1px solid var(--line);
  border-radius:6px; padding:12px; }}
figure img {{ display:block; width:100%; max-width:100%; height:auto; }}
figcaption {{ font-size:12px; color:var(--ink-3); margin-top:9px; }}
details {{ margin:14px 0; }}
summary {{ cursor:pointer; font-size:13px; color:var(--opalx); }}
.pair-cols {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(320px,1fr));
  gap:8px 28px; }}
.chip {{ font-family:var(--mono); font-size:9.5px; text-transform:uppercase;
  letter-spacing:.06em; background:var(--raised); color:var(--ink-3);
  padding:1px 6px; border-radius:9px; }}
a {{ color:var(--opalx); }}
:focus-visible {{ outline:2px solid var(--opalx); outline-offset:2px; }}
.key {{ display:flex; gap:20px; flex-wrap:wrap; font-size:12.5px; color:var(--ink-2);
  margin:8px 0 0; }}
.key span::before {{ content:""; display:inline-block; width:10px; height:10px;
  border-radius:2px; margin-right:6px; vertical-align:-1px; }}
.key .o::before {{ background:var(--opalx); }}
.key .g::before {{ background:var(--g4bl); }}
</style>

<header class="top">
  <h1>muE4 field maps: OPALX against G4beamline</h1>
  <p class="sub">The same particles, through the same magnetic field files, in both
  programs. The whole 19.4 m muE4 line, then nine single-magnet cases covering the five grid
  map files it uses. For each: how it is set up, what each program reports, and where they
  differ.</p>
  <div class="meta">
    <span>OPALX branch general-fieldmap-element</span>
    <span>G4beamline 3.08</span>
    <span>mu+ at 28 MeV/c</span>
    <span>vacuum, no material, no decay, no space charge</span>
  </div>
  <p class="key"><span class="g">G4beamline</span><span class="o">OPALX</span></p>
</header>

<div class="wrap">
<nav class="toc"><span class="lbl">Cases</span><ol>{toc}</ol></nav>
<main>
<section id="how">
<h2>How the comparison works</h2>
<p>Each case isolates one magnet. The same particle file is fed to both programs, so they
cannot be given different beams. Both read the <em>same</em> <code>.g4blmap</code> file —
nothing is converted. Three things are compared.</p>
<p><strong>The field, with no tracking.</strong> Both programs are asked for the field at
the same points. A difference here means the file is being read differently, or the magnet
is in a different place.</p>
<p><strong>Nineteen particles.</strong> One on the reference path, then a small step either
way in each of the six coordinates, then a wide step either way in x, y and momentum. The
small steps give a transfer matrix; the wide ones show where the map stops behaving
linearly.</p>
<p><strong>Twenty thousand particles.</strong> A Gaussian bunch, 10 mm across in x and 3 mm
in y, matched particle for particle between the two programs, so the comparison stays per
particle rather than only comparing beam sizes.</p>
<h3>What limits the comparison</h3>
<p>G4beamline writes six digits. When a particle is 1.27 m off centre that is steps of
10&nbsp;µm, and nothing smaller can be seen. The cases at the real muE4 angle keep the beam
within 35 mm of the reference line, so there the limit is 0.1&nbsp;µm — about a hundred
times finer. Every "allowed" figure in the tables below is derived from that limit, not
chosen to make a result pass.</p>
<h3>All nine cases</h3>
{summary_rows()}
<p class="note">Differences are the largest over all particles, which is a harsher measure
than the typical one. The median column gives the typical value for comparison.</p>
</section>
{spin_section()}
{full_line()}
{body}
<section id="limits">
<h2>What this does not cover</h2>
<p><strong>One map file is missing.</strong> muE4 uses six; these cases cover five. The
solenoid file <code>wsx_total.g4blmap</code> is the only one in the older, round format, and
it is the only file OPALX reads with different code. Its only comparison is a pair of older
notebooks that run at step sizes this report shows are too large, and that compare the field
only along the axis. Its off-axis field has never been compared with anything.</p>
<figure><img src="figs/step_convergence.png" alt="How the disagreement falls as the step
size is reduced, for both cases placed at the real muE4 angle"><figcaption>Halving the step
size roughly halves the disagreement, in both cases placed at the real muE4 angle. That is
what says the difference is the step size and not the field map. The lower line has settled;
the upper one, which crosses between two map files, has not.</figcaption></figure>
<p><strong>Two cases have not settled.</strong> <a href="#asr61_group_bisector">ASR61 group
at the real muE4 angle</a> still shows a difference that shrinks every time the step size is
halved, and the shrinking has not stopped. Whether any of what is left is a real difference
between the programs is not yet answered, and four of its checks are left failing rather
than given a looser limit.</p>
<p><strong>Nothing physical beyond the magnet.</strong> No material, no scattering, no
decay, no space charge, no apertures. One particle type, one charge, one momentum. The beam
is 10 mm across where muE4's pipe is 400 mm, so only the middle of each magnet is tested.</p>
<p><strong>Inside the field, the angle is unreliable.</strong> Both programs place a
recording plane's position correctly but report the momentum from whichever step the
particle happened to be on. In a field the momentum is changing, so a recording plane
sitting inside a magnet shows an angle difference of about one step — 0.26 mm worth on the
ASR61 magnet. Planes outside the field do not have this problem.</p>
</section>
</main>
</div>
"""
(HERE / "index.html").write_text(HTML)
print(f"wrote index.html  {len(HTML)/1024:.0f} KB")
