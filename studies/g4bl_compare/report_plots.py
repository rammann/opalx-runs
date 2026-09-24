#!/usr/bin/env python3
"""Three figures per case for the report, into report/figs/.

  <case>_field  the field both codes report, and where they differ
  <case>_pair   the 19 particles, one by one, and the transfer matrix
  <case>_gauss  the 20000 particles: beam size and centre, and the spread of differences
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from opalxruns import plotstyle
from opalxruns.g4bl import read_map_header
from opalxruns.plotstyle import G4BL, INK, MUTED, OPALX

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmplib

FIGS = HERE / "report" / "figs"
FIGS.mkdir(parents=True, exist_ok=True)
plotstyle.use({"axes.labelcolor": INK})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGS / f"{name}.png", dpi=110)
    plt.close(fig)
    print(f"  {name}.png")


def planes(d, c, sub):
    out = {}
    for z in c["monitors"]:
        try:
            dg = cmplib.read_bltrack(d / sub / f"Z{z}.txt")
            do = cmplib.read_monitor(d / sub / f"MON_{z}.h5")
            out[z] = cmplib.matched(dg, do)[:2] + (dg,)
        except Exception:
            pass
    return out


def _map_extent(fname):
    """How far the map reaches in x and z, in its own frame, from its header line."""
    header = read_map_header(cmplib.MAPS / fname)
    if header is None:
        return None
    kind, d = header
    if kind == "grid":
        return ((d["X0"], d["X0"] + (d["nX"] - 1) * d["dX"]),
                (d["Z0"], d["Z0"] + (d["nZ"] - 1) * d["dZ"]))
    r = d["nR"] * d["dR"] - d["dR"]
    z0 = d["Z0"]
    return ((-r, r), (z0, z0 + (d["nZ"] - 1) * d["dZ"]))


def _geometry(name, d):
    """Where each field map and each recording plane sits, read from the input OPALX
    was actually given rather than from a description of it."""
    txt = (d / f"{name}.in").read_text().splitlines()
    maps, mons = [], []
    for i, ln in enumerate(txt):
        m = re.match(r"E\d+_\S*\s*: (FIELDMAP|MONITOR), X = ([-\d.]+), Y = ([-\d.]+), "
                     r"Z = ([-\d.]+),", ln)
        if not m:
            continue
        th = re.search(r"THETA = ([-\d.]+)", txt[i + 1])
        rec = dict(x=float(m.group(2)), z=float(m.group(4)),
                   th=float(th.group(1)) if th else 0.0)
        if m.group(1) == "FIELDMAP":
            fn = re.search(r'FMAPFN = "([^"]+)"', txt[i + 2])
            rec["file"] = Path(fn.group(1)).name if fn else ""
            maps.append(rec)
        else:
            mons.append(rec)
    return maps, mons


def _footprint(mp, xr, zr):
    """The map's outline seen from above, put where the element is placed. Only x and
    z matter looking down; y is the height and is left out."""
    pts = np.array([[xr[0], zr[0]], [xr[1], zr[0]], [xr[1], zr[1]],
                    [xr[0], zr[1]], [xr[0], zr[0]]]) * 1e-3
    c, s = np.cos(mp["th"]), np.sin(mp["th"])
    out = np.empty_like(pts)
    out[:, 0] = mp["x"] + pts[:, 0] * c + pts[:, 1] * s      # lab x
    out[:, 1] = mp["z"] - pts[:, 0] * s + pts[:, 1] * c      # lab z
    return out


def fig_layout(name, c, d):
    """A plan view: the map, the path the beam takes through it, and the planes."""
    maps, mons = _geometry(name, d)
    if not maps:
        return
    sub = "fine" if (c.get("converge") and (d / "fine").is_dir()) else "pair"
    stat = next((d / sub).glob("*.stat"), None)
    fig, ax = plt.subplots(figsize=(9.6, 4.6))

    ext = {mp["file"]: _map_extent(mp["file"]) for mp in maps}
    shades = ["#2a78d6", "#eb6834", "#1baf7a"]
    drawn = set()
    for k, mp in enumerate(maps):
        e = ext.get(mp["file"])
        if e is None:
            continue
        pts = _footprint(mp, *e)
        col = shades[k % len(shades)]
        lab = mp["file"] if mp["file"] not in drawn else None
        drawn.add(mp["file"])
        ax.plot(pts[:, 1], pts[:, 0], "-", color=col, lw=1.2, zorder=3)
        ax.fill(pts[:, 1], pts[:, 0], color=col, alpha=0.10, zorder=1, label=lab)

    ref = cmplib.read_stat(stat) if stat is not None else None
    if ref is not None:
        ax.plot(np.asarray(ref["ref_z"]), np.asarray(ref["ref_x"]), "-", color=INK,
                lw=1.5, zorder=5, label="the path the beam takes")

    P = planes(d, c, sub)
    # A recording plane is a plane, not a short bar: it catches the beam wherever the
    # beam happens to be. Draw it long enough to reach across everything else drawn,
    # or it looks as though the beam misses it.
    reach = [abs(mp["x"]) for mp in maps] + [abs(pl["x"]) for pl in mons]
    if ref is not None:
        reach.append(float(np.abs(np.asarray(ref["ref_x"])).max()))
    W = 0.6 * (max(reach) + 0.25)
    for i, pl in enumerate(mons):
        cth, sth = np.cos(pl["th"]), np.sin(pl["th"])
        ax.plot([pl["z"] - W * sth, pl["z"] + W * sth],
                [pl["x"] + W * cth, pl["x"] - W * cth], "-", color="#4a3aa7", lw=1.6,
                zorder=6, alpha=0.65,
                label="recording plane" if i == 0 else None)
    # where the 19 particles actually crossed each plane, in both codes
    for i, z in enumerate(sorted(P)):
        Xg, Xo, _ = P[z]
        pl = min(mons, key=lambda q: abs(q["z"] * 1e3 - z))
        cth, sth = np.cos(pl["th"]), np.sin(pl["th"])
        for X, col, lab in ((Xg, G4BL, "G4beamline"), (Xo, OPALX, "OPALX")):
            ax.scatter(pl["z"] + X[0] * sth, pl["x"] + X[0] * cth, s=8, color=col,
                       alpha=0.85, zorder=7, edgecolors="none",
                       label=lab if i == 0 else None)

    ax.set(xlabel="lab z [m]", ylabel="lab x [m]",
           title=f"{name} seen from above: the field map, the path through it, "
                 f"and the recording planes")
    ax.set_aspect("equal", adjustable="box")
    ax.relim(); ax.autoscale_view()
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    mx, my = 0.05 * (x1 - x0), 0.12 * (y1 - y0)
    ax.set_xlim(x0 - mx, x1 + mx); ax.set_ylim(y0 - my, y1 + my)
    ax.legend(fontsize=6.5, ncol=2, loc="upper left")
    save(fig, f"{name}_layout")


def fig_field(name, c, d):
    keys = [k for k in c["grids"] if k != "edge"]
    fig, ax = plt.subplots(1, 3, figsize=(12.0, 3.1))
    main = keys[0]
    fg = d / "pair" / f"g4bl_field_{main}.txt"
    if not fg.exists():
        return
    pos, Bg, Bo = cmplib.match_field(cmplib.read_g4bl_field(fg),
                                     cmplib.read_opalx_field(d / "pair" / f"opalx_field_{main}.dat"))
    keep = ~cmplib.on_face(pos, c["boxes"])
    xs = np.unique(pos[0])
    sel = keep & np.isclose(pos[0], xs[np.argmin(np.abs(xs))])
    o = np.argsort(pos[2][sel])
    z = pos[2][sel][o] * 1000
    comp = int(np.argmax(np.abs(Bg[:, sel]).max(axis=1)))
    lbl = "Bx By Bz".split()[comp]
    ax[0].plot(z, Bg[comp, sel][o], color=G4BL, lw=1.8, label="G4beamline")
    ax[0].plot(z, Bo[comp, sel][o], color=OPALX, lw=1.0, ls="--", label="OPALX")
    ax[0].set(xlabel="lab z [mm]", ylabel=f"{lbl} [T]",
              title=f"{lbl} along the sampled line")
    ax[0].legend()
    ax[1].plot(z, (Bo[comp, sel] - Bg[comp, sel])[o] * 1e9, color=INK, lw=0.9)
    ax[1].axhline(0, color=MUTED, lw=0.5)
    ax[1].set(xlabel="lab z [mm]", ylabel="OPALX - G4beamline [nT]",
              title="difference along the same line")
    diff = np.abs(Bo - Bg).max(axis=0)
    zs = np.unique(pos[2])
    img = np.full((len(xs), len(zs)), np.nan)
    img[np.searchsorted(xs, pos[0]), np.searchsorted(zs, pos[2])] = np.where(keep, diff, np.nan)
    fl = max(np.nanmin(img[img > 0]) if np.any(img > 0) else 1e-9, 1e-10)
    im = ax[2].imshow(np.maximum(img, fl), origin="lower", aspect="auto", cmap="Blues",
                      norm=LogNorm(vmin=fl, vmax=max(np.nanmax(img), fl * 10)),
                      extent=[zs[0] * 1000, zs[-1] * 1000, xs[0] * 1000, xs[-1] * 1000])
    ax[2].set(xlabel="lab z [mm]", ylabel="lab x [mm]",
              title=f"|difference| over the plane\nworst {np.nanmax(img):.1e} T")
    ax[2].grid(False)
    fig.colorbar(im, ax=ax[2]).set_label("|OPALX - G4BL| [T]", fontsize=7)
    save(fig, f"{name}_field")


def fig_pair(name, c, d):
    sub = "fine" if (c.get("converge") and (d / "fine").is_dir()) else "pair"
    p = planes(d, c, sub)
    if len(p) < 2:
        return
    zi, zo = min(p), max(p)
    Xg_i, Xo_i, _ = p[zi]
    Xg_o, Xo_o, dg_o = p[zo]
    fig, ax = plt.subplots(1, 3, figsize=(12.0, 3.2))
    ids = np.arange(Xg_o.shape[1])
    dpos = np.abs(Xo_o[[0, 2]] - Xg_o[[0, 2]]).max(axis=0) * 1e6
    dang = np.abs(Xo_o[[1, 3]] - Xg_o[[1, 3]]).max(axis=0) * 1e6
    w = 0.4
    ax[0].bar(ids - w / 2, np.maximum(dpos, 1e-6), w, color=OPALX, label="position [um]")
    ax[0].bar(ids + w / 2, np.maximum(dang, 1e-6), w, color=G4BL, label="angle [urad]")
    if len(ids) > 13:
        ax[0].axvline(12.5, color=MUTED, lw=0.8)
        ax[0].text(12.6, ax[0].get_ylim()[1], " wide steps", fontsize=6, color=MUTED, va="top")
    ax[0].set(yscale="log", xlabel="particle", title=f"exit difference, plane z = {zo} mm")
    ax[0].legend(fontsize=6)
    ax[1].plot(ids, Xg_o[0] * 1000, "o", ms=4, color=G4BL, label="G4beamline")
    ax[1].plot(ids, Xo_o[0] * 1000, "x", ms=5, color=OPALX, label="OPALX")
    ax[1].set(xlabel="particle", ylabel="exit x [mm]", title="where each particle came out")
    ax[1].legend(fontsize=6)
    Mg = cmplib.transfer_matrix(Xg_i, Xg_o)
    Mo = cmplib.transfer_matrix(Xo_i, Xo_o)
    fl = cmplib.row_floors(dg_o, cmplib.STEP_SMALL)
    r = np.abs(Mo - Mg) / fl[:, None]
    im = ax[2].imshow(r, cmap="Blues", vmin=0, vmax=max(r.max(), 1e-9))
    ax[2].set_xticks(range(6), cmplib.COORD_NAMES, fontsize=6)
    ax[2].set_yticks(range(6), cmplib.COORD_NAMES, fontsize=6)
    ax[2].set_title(f"transfer matrix difference\nin printing steps, worst {r.max():.1f}")
    ax[2].grid(False)
    fig.colorbar(im, ax=ax[2])
    save(fig, f"{name}_pair")


def fig_gauss(name, c, d):
    p = planes(d, c, "gauss")
    if not p:
        return
    zs = sorted(p)
    fig, ax = plt.subplots(1, 4, figsize=(15.2, 3.2))
    for j, lab, ls in ((0, "x", "-"), (2, "y", "--")):
        ax[0].plot(zs, [p[z][0][j].std() * 1000 for z in zs], ls, color=G4BL, lw=1.8,
                   label=f"rms {lab}, G4beamline")
        ax[0].plot(zs, [p[z][1][j].std() * 1000 for z in zs], ls, color=OPALX, lw=1.0,
                   marker="o", ms=3, label=f"rms {lab}, OPALX")
    ax[0].set(xlabel="lab z of the plane [mm]", ylabel="rms [mm]", title="beam size")
    ax[0].legend(fontsize=5.5, ncol=2)
    for j, lab, ls in ((0, "x", "-"), (2, "y", "--")):
        ax[1].plot(zs, [abs(p[z][1][j].mean() - p[z][0][j].mean()) * 1e6 for z in zs],
                   ls, marker="o", ms=3, color=OPALX, label=f"mean {lab}")
        ax[1].plot(zs, [abs(p[z][1][j].std() - p[z][0][j].std()) * 1e6 for z in zs],
                   ls, marker="s", ms=3, color=G4BL, label=f"rms {lab}")
    ax[1].set(xlabel="lab z of the plane [mm]", ylabel="|OPALX - G4BL| [um]",
              yscale="log", title="difference in beam centre and size")
    ax[1].legend(fontsize=5.5, ncol=2)
    z0, z1 = zs[0], zs[-1]
    Xg0, _, _ = p[z0]
    Xg1, Xo1, _ = p[z1]
    n = min(Xg0.shape[1], Xg1.shape[1])
    r = np.hypot(Xg0[0][:n], Xg0[2][:n]) * 1000
    dv = np.abs(Xo1[[0, 2]][:, :n] - Xg1[[0, 2]][:, :n]).max(axis=0) * 1e6
    ax[2].scatter(r, np.maximum(dv, 1e-4), s=3, alpha=0.15, color=OPALX, edgecolors="none")
    b = np.linspace(0, r.max(), 21)
    idx = np.clip(np.digitize(r, b) - 1, 0, len(b) - 2)
    med = [np.median(dv[idx == i]) if (idx == i).any() else np.nan for i in range(len(b) - 1)]
    ax[2].plot(0.5 * (b[1:] + b[:-1]), med, color=G4BL, lw=1.8, label="median")
    ax[2].set(xlabel="starting radius [mm]", ylabel="|OPALX - G4BL| [um]", yscale="log",
              title=f"per particle at z = {z1} mm")
    ax[2].legend(fontsize=6)
    # the same numbers as a distribution, which says how many particles sit where
    # rather than only where the worst one is
    bins = np.logspace(-4, 3, 40)
    ax[3].hist(np.clip(dv, bins[0], bins[-1]), bins=bins, color=OPALX, alpha=0.85,
               edgecolor="none")
    med = float(np.median(dv))
    ax[3].axvline(med, color=G4BL, lw=1.4)
    ax[3].text(med, ax[3].get_ylim()[1] * 0.92, f"  median {med:.3g} um", fontsize=6.5,
               color=G4BL, va="top")
    ax[3].set(xscale="log", yscale="log", xlabel="|OPALX - G4BL| [um]",
              ylabel="particles", title="the same, as a distribution")
    save(fig, f"{name}_gauss")


def fig_convergence(man):
    cases = [n for n, c in man.items() if c.get("converge")]
    if not cases:
        return
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))
    data = {
        "asr61_bisector": ([1.0, 0.5], [0.520, 0.246], [0.327, 0.198]),
        "asr61_group_bisector": ([1.0, 0.5, 0.25], [2.282, 1.201, 0.747],
                                 [0.970, 0.578, 0.395]),
    }
    for k, (n, (h, pos, ang)) in enumerate(data.items()):
        ax[0].plot(h, pos, "o-", lw=1.6, color=[OPALX, G4BL][k], label=n)
        ax[1].plot(h, ang, "o-", lw=1.6, color=[OPALX, G4BL][k], label=n)
    for a, t, yl in ((ax[0], "worst position difference", "[um]"),
                     (ax[1], "worst angle difference", "[urad]")):
        a.set(xscale="log", yscale="log", xlabel="step size, relative to 1e-12 s / 0.1 mm",
              ylabel=yl, title=t)
        a.legend(fontsize=6)
    save(fig, "step_convergence")


def main():
    man = cmplib.load_manifest()
    for name, c in man.items():
        d = HERE / c["dir"]
        print(name)
        for fn in (fig_layout, fig_field, fig_pair, fig_gauss):
            try:
                fn(name, c, d)
            except Exception as e:
                print(f"  {fn.__name__} failed: {type(e).__name__}: {e}")
    fig_convergence(man)


if __name__ == "__main__":
    main()
