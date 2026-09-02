#!/usr/bin/env python3
"""Write the shared particle set for asr61_dipole.

Two files, the same 13 muons in both:
  parts.txt  OPALX FROMFILE      x y z px py pz, positions in m, momenta in beta*gamma
  beam.txt   G4BL #BLTrackFile   x y z Px Py Pz t PDGid EventID ..., mm and MeV/c

Particle 0 is on axis; the rest probe the transverse and momentum acceptance of
the map, staying well inside its -200..+390 mm x and -120..+120 mm y box at the
entrance.
"""

MUON_MASS = 0.1056583755  # GeV
P0 = 0.028                # GeV/c, the muE4 reference momentum
Z0_MM = 0.0               # start plane, global z [mm]

# (x [mm], y [mm], dp/p)
OFFSETS = [
    (0.0, 0.0, 0.0),
    (10.0, 0.0, 0.0),
    (-10.0, 0.0, 0.0),
    (0.0, 10.0, 0.0),
    (0.0, -10.0, 0.0),
    (20.0, 10.0, 0.0),
    (-20.0, -10.0, 0.0),
    (0.0, 0.0, 0.02),
    (0.0, 0.0, -0.02),
    (10.0, 10.0, 0.01),
    (-10.0, 10.0, -0.01),
    (20.0, -20.0, 0.0),
    (-20.0, 20.0, 0.0),
]

with open("parts.txt", "w") as f:
    f.write(f"{len(OFFSETS)}\n")
    f.write("x y z px py pz\n")
    for x, y, dpp in OFFSETS:
        pz = P0 * (1.0 + dpp) / MUON_MASS  # beta*gamma
        f.write(f"{x * 1e-3:.9e} {y * 1e-3:.9e} {Z0_MM * 1e-3:.9e} 0.0 0.0 {pz:.9e}\n")

with open("beam.txt", "w") as f:
    f.write("#BLTrackFile written by make_parts.py\n")
    f.write("#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n")
    f.write("#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n")
    for i, (x, y, dpp) in enumerate(OFFSETS):
        pz = P0 * (1.0 + dpp) * 1000.0  # MeV/c
        f.write(f"{x:.9g} {y:.9g} {Z0_MM:.9g} 0 0 {pz:.9g} 0 -13 {i + 1} 1 0 1\n")

print("wrote parts.txt and beam.txt")
