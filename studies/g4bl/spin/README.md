# spin — spin precession in OPALX and G4beamline against the closed form

| cases | what they check |
|---|---|
| `by_p0028` ... `by_p3094` | the g-2 test: a muon bends through a uniform By at five momenta (0.028 to 3.094 GeV/c) with the bend angle held fixed; the spin leads the momentum by G*gamma*theta, and a straight-line fit in gamma gives the muon anomaly G each code is using |
| `uniform_bz/` | a muon along a uniform Bz: no force, so only the spin turns, by (1+G) times the cyclotron angle. `scan/dt_*.in` repeat the OPALX run at smaller time steps, `bz2m.in` over 2 m of field |

`make_scan.py` writes the `by_p*` cases (and their field maps, which are not tracked);
`uniform_bz/` is written by hand. The `scan/` inputs read the field map and particles of
`uniform_bz/`.

```bash
OPALX=/path/to/opalx ./run_all.sh   # write the cases, run both codes, analyse
./run_all.sh --test-only            # analyse the existing output only
```

Output goes to `output/g4bl/spin/`: one folder per `by_p*` case, and for `uniform_bz/` one per
input (`uniform_bz`, `dt_<dt>`, `bz2m`). `analyse_scan.py` writes `results.txt` and
`scan_data.json`, `analyse_bz.py` writes `bz_results.txt` and `bz_data.json`; the G4BL report
(`../report/`) reads those JSON files. Figures go to `output/g4bl/spin/plots/`.
