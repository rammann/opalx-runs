"""Tests for opalxruns.case.Case, on a small H5Part file written here.

Run from the repo root:  python -m unittest discover -s tests
"""

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from opalxruns.case import Case

BG0 = 2.0


def write_h5(path):
    """Two dumps written out of order, particles stored out of id order."""
    with h5py.File(path, "w") as f:
        for step, spos, shift in ((1, 0.8, 0.1), (0, 0.2, 0.0)):
            g = f.create_group(f"Step#{step}")
            g.attrs["SPOS"] = np.array([spos])
            g["id"] = np.array([2, 0, 1])
            g["x"] = np.array([0.3, 0.1, 0.2]) + shift
            g["y"] = np.array([0.0, 0.0, 0.0])
            g["z"] = np.array([0.01, 0.0, -0.01])
            g["px"] = np.array([0.2, 0.0, 0.1])
            g["py"] = np.array([0.0, 0.0, 0.0])
            g["pz"] = np.array([2.0, 2.0, 2.0])


class CaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        folder = Path(self.tmp.name) / "one_case"
        folder.mkdir()
        write_h5(folder / "one_case.h5")
        self.case = Case(folder, {"h5": "one_case.h5", "stat": "one_case.stat"}, BG0)

    def tearDown(self):
        self.tmp.cleanup()

    def test_name_and_files_come_from_the_folder_and_entry(self):
        self.assertEqual(self.case.name, "one_case")
        self.assertEqual(self.case.h5.name, "one_case.h5")
        self.assertEqual(self.case.stat.name, "one_case.stat")

    def test_dumps_sorted_by_step_number_with_their_path_length(self):
        steps, sp = self.case._steps_spos()
        self.assertEqual(steps, [0, 1])
        np.testing.assert_array_equal(sp, [0.2, 0.8])

    def test_read_plane_sorts_by_id_and_projects_to_z_zero(self):
        X = self.case.read_plane(0)
        xp = np.array([0.0, 0.05, 0.1])
        z = np.array([0.0, -0.01, 0.01])
        np.testing.assert_allclose(X[0], np.array([0.1, 0.2, 0.3]) - xp * z)
        np.testing.assert_allclose(X[1], xp)
        np.testing.assert_allclose(X[4], z)
        np.testing.assert_allclose(X[5], np.sqrt(np.array([0.0, 0.1, 0.2])**2 + 4.0) / BG0 - 1.0)

    def test_read_plane_without_projection(self):
        np.testing.assert_allclose(self.case.read_plane(0, project=False)[0], [0.1, 0.2, 0.3])

    def test_trajectory_follows_one_particle_through_every_dump(self):
        s, x, y = self.case.trajectory(0)
        np.testing.assert_allclose(s, [0.2, 0.8])
        np.testing.assert_allclose(x, [0.1, 0.2])

    def test_momentum_magnitude_per_particle(self):
        np.testing.assert_allclose(self.case._pmag(0), np.sqrt(np.array([0.0, 0.1, 0.2])**2 + 4.0))


if __name__ == "__main__":
    unittest.main()
