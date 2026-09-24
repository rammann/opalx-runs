"""Tests for the G4beamline file helpers in opalxruns.g4bl.

Run from the repo root:  python -m unittest discover -s tests
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from opalxruns import g4bl

TRACK = """#BLTrackFile2 Z1200
#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight
10.0 -2.0 1199.5 1.0 0.5 20.0 4.0 -13 1 1 0 1.0

11.0 -3.0 1200.0 0.0 0.0 25.0 4.1 -13 2 1 0 1.0
"""


class TrackFileTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "Z1200.txt"
        self.path.write_text(TRACK)

    def tearDown(self):
        self.tmp.cleanup()

    def test_rows_skip_comments_and_blank_lines(self):
        a = g4bl.read_track_file(self.path)
        self.assertEqual(a.shape, (2, 12))
        self.assertEqual(a[1, 8], 2.0)

    def test_track_on_plane_drifts_to_the_plane_and_keys_by_eventid_minus_one(self):
        d = g4bl.track_on_plane(self.path, 1200.0)
        self.assertEqual(sorted(d), [0, 1])
        # particle 0 sits 0.5 mm before the plane with x' = 1/20, y' = 0.5/20
        np.testing.assert_allclose(d[0], [10.0 + 0.5 / 20.0, -2.0 + 0.25 / 20.0, 1.0, 0.5, 20.0])
        np.testing.assert_allclose(d[1], [11.0, -3.0, 0.0, 0.0, 25.0])


class MapHeaderTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_grid_header_from_a_written_map(self):
        path = self.dir / "quad.g4blmap"
        axis = [-10.0, 0.0, 10.0]
        g4bl.write_grid_map(path, lambda x, y, z: (y, x, 0.0), axis, axis, [-5.0, 5.0])
        kind, h = g4bl.read_map_header(path)
        self.assertEqual(kind, "grid")
        self.assertEqual((h["X0"], h["nX"], h["dX"], h["nZ"], h["dZ"]), (-10.0, 3.0, 10.0, 2.0, 10.0))

    def test_cylinder_header(self):
        path = self.dir / "sol.g4blmap"
        g4bl.write_cylinder_map(path, lambda r, z: (1.0, 0.0), [0.0, 1.0, 2.0], [0.0, 5.0])
        kind, h = g4bl.read_map_header(path)
        self.assertEqual(kind, "cylinder")
        self.assertEqual((h["nR"], h["dR"], h["nZ"]), (3.0, 1.0, 2.0))

    def test_file_without_a_header(self):
        path = self.dir / "empty.g4blmap"
        path.write_text("param normB=1\n")
        self.assertIsNone(g4bl.read_map_header(path))


if __name__ == "__main__":
    unittest.main()
