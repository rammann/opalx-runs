"""Tests for opalxruns.plotstyle.

Run from the repo root:  python -m unittest discover -s tests
"""

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from opalxruns import plotstyle  # noqa: E402


class UseTest(unittest.TestCase):
    def setUp(self):
        self.saved = dict(plt.rcParams)

    def tearDown(self):
        plt.rcParams.update(self.saved)

    def test_shared_style_is_applied(self):
        plotstyle.use()
        self.assertEqual(plt.rcParams["font.size"], 8)
        self.assertEqual(plt.rcParams["grid.color"], plotstyle.GRID)
        self.assertTrue(plt.rcParams["axes.grid"])

    def test_changes_override_the_shared_style(self):
        plotstyle.use({"axes.titlesize": 10, "grid.color": "#e4e4e4"})
        self.assertEqual(plt.rcParams["axes.titlesize"], 10)
        self.assertEqual(plt.rcParams["grid.color"], "#e4e4e4")


class SaveTest(unittest.TestCase):
    def test_writes_the_file_and_closes_the_figure(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "plots" / "a.png"
            fig = plt.figure()
            plotstyle.save_with_footer(fig, path, "footer text")
            self.assertTrue(path.is_file())
            self.assertFalse(plt.fignum_exists(fig.number))


if __name__ == "__main__":
    unittest.main()
