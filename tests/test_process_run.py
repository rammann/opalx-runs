"""Tests for how opalxruns.process_run picks the run folder.

Run from the repo root:  python -m unittest discover -s tests
"""

import unittest
from pathlib import Path

from opalxruns import paths, process_run


class RunFolderTest(unittest.TestCase):
    def test_case_folder_under_studies_means_its_run_folder_in_output(self):
        case = paths.STUDIES / "features" / "ring" / "square_ring"
        self.assertEqual(process_run.run_folder(case), paths.output_dir(case))

    def test_any_other_folder_is_used_as_it_is(self):
        self.assertEqual(process_run.run_folder(Path("/tmp/somewhere")), Path("/tmp/somewhere"))


if __name__ == "__main__":
    unittest.main()
