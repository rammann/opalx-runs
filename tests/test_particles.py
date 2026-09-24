"""Tests for opalxruns.particles.

Run from the repo root:  python -m unittest discover -s tests
"""

import math
import tempfile
import unittest
from pathlib import Path

from opalxruns import particles

EPS = {"x": 1e-4, "xp": 1e-4, "y": 1e-4, "yp": 1e-4, "z": 1e-4, "delta": 1e-3}


class MapParticlesTest(unittest.TestCase):
    def test_reference_plus_a_pair_per_coordinate(self):
        rows = particles.map_particles(0.5, EPS)
        self.assertEqual([label for label, _ in rows],
                         ["ref", "x+", "x-", "xp+", "xp-", "y+", "y-", "yp+", "yp-",
                          "z+", "z-", "delta+", "delta-"])
        self.assertEqual(rows[0][1], [0, 0, 0, 0, 0, 0.5])

    def test_momentum_magnitude_carries_only_the_delta_step(self):
        for label, (x, px, y, py, z, pz) in particles.map_particles(0.5, EPS):
            delta = {"delta+": EPS["delta"], "delta-": -EPS["delta"]}.get(label, 0.0)
            self.assertAlmostEqual(math.sqrt(px * px + py * py + pz * pz), 0.5 * (1 + delta), 15)

    def test_angle_step_sets_the_slope(self):
        rows = dict(particles.map_particles(0.5, EPS))
        x, px, y, py, z, pz = rows["xp+"]
        self.assertAlmostEqual(px / pz, math.tan(EPS["xp"]), 15)


class WritePartsTest(unittest.TestCase):
    def test_count_header_and_row_format(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "parts.txt"
            particles.write_parts(path, [("ref", [0, 0, 0, 0, 0, 0.5])])
            self.assertEqual(path.read_text(),
                             "1\nx px y py z pz\n" + " ".join(f"{c:.12e}" for c in [0, 0, 0, 0, 0, 0.5])
                             + "\n")


if __name__ == "__main__":
    unittest.main()
