"""Tests for the Step#N listing shared by opalx_diagnostics and opalx_run.

Run from the repo root:  python -m unittest discover -s tests
"""

import tempfile
import unittest
from pathlib import Path

import h5py

from opalxruns.opalx_diagnostics import list_h5_steps
from opalxruns.opalx_run import h5_steps


class StepListTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "run.h5"
        with h5py.File(self.path, "w") as f:
            for name in ("Step#10", "Step#2", "Step#0", "Step#extra", "other"):
                f.create_group(name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_numbers_in_numeric_order_and_other_groups_ignored(self):
        self.assertEqual(list_h5_steps(self.path), [0, 2, 10])

    def test_names_in_the_same_order(self):
        self.assertEqual(h5_steps(self.path), ["Step#0", "Step#2", "Step#10"])


if __name__ == "__main__":
    unittest.main()
