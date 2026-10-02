# maps

G4beamline field maps kept in this repo, so that a clone on another machine (daint) can
run without the `g4bl-files` folder.

| file | elements in `studies/g4bl/mue4/daint_scaling/lattice_scaling.in` |
|---|---|
| `wsx_total.g4blmap` | `E05_WSX` |
| `asr61_300d_track.g4blmap` | `E08_ASR61_300d` |
| `asr61_300sm_track.g4blmap` | `E09_ASR61_300sm`, `E10_ASR61_300sm` |
| `asr62shim_280_track.g4blmap` | `E19_ASR62`, `E30_ASR62` |
| `asr62shim_280_sm_track.g4blmap` | `E20_ASR62_sm`, `E21_ASR62_sm`, `E31_ASR62_sm`, `E32_ASR62_sm` |
| `qsm01a_210_track.g4blmap` | the 12 `QSM600` elements |

These are the six maps the muE4 line uses, copied unchanged from
`g4bl-files/muE4/maps/` on 2026-10-02. Three of them are over 50 MB, the size above
which GitHub warns (it refuses files over 100 MB).

Used by `studies/g4bl/mue4/daint_scaling/`. The other muE4 setups still read the maps
from `g4bl-files` (`G4BL_FILES` in `opalxruns/paths.py`).
