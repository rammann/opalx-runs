"""Tests for opalxruns.paths.

Run from the repo root:  python -m unittest discover -s tests
"""

import importlib
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from opalxruns import paths


class RepoTest(unittest.TestCase):
    def test_repo_is_the_folder_that_holds_the_package(self):
        self.assertTrue((paths.REPO / "pyproject.toml").is_file())
        self.assertTrue((paths.REPO / "opalxruns" / "paths.py").is_file())


class G4blFilesTest(unittest.TestCase):
    def tearDown(self):
        importlib.reload(paths)

    def test_default_is_the_folder_next_to_the_repo(self):
        # Generators write this path into tracked input files, so it must keep
        # the spelling they used before: /Users/rammann/Code/OPALX/g4bl-files.
        with mock.patch.dict(os.environ):
            os.environ.pop("G4BL_FILES", None)
            importlib.reload(paths)
            self.assertEqual(paths.G4BL_FILES, (paths.REPO.parent / "g4bl-files").resolve())

    def test_environment_variable_overrides_the_default(self):
        with mock.patch.dict(os.environ, {"G4BL_FILES": "/tmp"}):
            importlib.reload(paths)
            self.assertEqual(paths.G4BL_FILES, Path("/tmp").resolve())


class G4blAppTest(unittest.TestCase):
    def tearDown(self):
        importlib.reload(paths)

    def test_default_is_the_installed_app(self):
        with mock.patch.dict(os.environ):
            os.environ.pop("G4BL_APP", None)
            importlib.reload(paths)
            self.assertEqual(paths.G4BL_APP, Path("/Users/rammann/Code/G4BL/G4beamline-3.08.app"))

    def test_environment_variable_overrides_the_default(self):
        with mock.patch.dict(os.environ, {"G4BL_APP": "/opt/other.app"}):
            importlib.reload(paths)
            self.assertEqual(paths.G4BL_APP, Path("/opt/other.app"))


class OpalxTest(unittest.TestCase):
    def test_unset_variable_is_an_error_that_names_it(self):
        with mock.patch.dict(os.environ):
            os.environ.pop("OPALX", None)
            with self.assertRaisesRegex(RuntimeError, "OPALX"):
                paths.opalx()

    def test_path_that_is_not_an_executable_is_an_error(self):
        with tempfile.NamedTemporaryFile() as f:           # exists, but not executable
            for value in ("/no/such/opalx", f.name):
                with self.subTest(value=value), mock.patch.dict(os.environ, {"OPALX": value}):
                    with self.assertRaisesRegex(RuntimeError, "not an executable"):
                        paths.opalx()

    def test_executable_is_returned(self):
        with tempfile.TemporaryDirectory() as d:
            exe = Path(d) / "opalx"
            exe.write_text("#!/bin/sh\n")
            exe.chmod(0o755)
            with mock.patch.dict(os.environ, {"OPALX": str(exe)}):
                self.assertEqual(paths.opalx(), exe)


if __name__ == "__main__":
    unittest.main()
