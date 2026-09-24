"""Tests for opalxruns.run: run folders under output/, links to inputs, launching.

Run from the repo root:  python -m unittest discover -s tests
"""

import contextlib
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from opalxruns import paths, run


class OutputDirTest(unittest.TestCase):
    def test_output_mirrors_studies(self):
        case = paths.STUDIES / "elements" / "collimator" / "circle"
        self.assertEqual(paths.output_dir(case), paths.OUTPUT / "elements" / "collimator" / "circle")

    def test_folder_outside_studies_is_an_error(self):
        with self.assertRaises(ValueError):
            paths.output_dir(Path("/tmp"))


class PrepareTest(unittest.TestCase):
    """A small repo: studies/ and output/ side by side, plus shared/."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.output = root / "output"
        self.topic = root / "studies" / "g" / "topic"
        self.case = self.topic / "case"
        (self.case / "maps").mkdir(parents=True)
        (self.topic / "lattice.in").write_text('F: FIELDMAP, FMAPFN = "../../../../shared/b.map";\n')
        (root / "shared").mkdir()
        (root / "shared" / "b.map").write_text("")
        (self.case / "maps" / "a.map").write_text("")
        (self.case / "parts.txt").write_text("")
        self.input = self.case / "case.in"
        self.input.write_text('S: SOLENOID, FMAPFN = "maps/a.map";\n'
                              'D: DISTRIBUTION, TYPE = FROMFILE, FNAME = "parts.txt";\n'
                              'CALL, FILE = "../lattice.in";\n')
        self.run_dir = self.output / "g" / "topic" / "case"

    def tearDown(self):
        self.tmp.cleanup()

    def test_links_the_input_and_every_file_it_names(self):
        run.prepare(self.run_dir, [self.input], root=self.output)
        for rel, target in (("case.in", self.input), ("maps/a.map", self.case / "maps" / "a.map"),
                            ("parts.txt", self.case / "parts.txt"),
                            ("../lattice.in", self.topic / "lattice.in")):
            link = Path(os.path.normpath(self.run_dir / rel))
            with self.subTest(rel=rel):
                self.assertTrue(link.is_symlink())
                self.assertEqual(link.resolve(), target.resolve())

    def test_a_path_that_already_reaches_the_file_is_not_linked(self):
        # ../../../../shared/b.map from output/g/topic/case is <root>/shared/b.map
        run.prepare(self.run_dir, [self.input], root=self.output)
        self.assertFalse((self.output / "shared").exists())
        self.assertTrue((self.run_dir / "../../../../shared/b.map").exists())

    def test_old_output_is_removed_unless_kept(self):
        self.run_dir.mkdir(parents=True)
        (self.run_dir / "case.h5").write_text("stale")
        (self.run_dir / "data").mkdir()
        run.prepare(self.run_dir, [self.input], keep=True, root=self.output)
        self.assertTrue((self.run_dir / "case.h5").exists())
        run.prepare(self.run_dir, [self.input], root=self.output)
        self.assertFalse((self.run_dir / "case.h5").exists())
        self.assertFalse((self.run_dir / "data").exists())

    def test_run_folder_outside_output_is_refused(self):
        with self.assertRaisesRegex(ValueError, "outside"):
            run.prepare(self.case / "elsewhere", [self.input], root=self.output)

    def test_missing_file_is_an_error_that_names_it(self):
        (self.case / "parts.txt").unlink()
        with self.assertRaisesRegex(FileNotFoundError, "parts.txt"):
            run.prepare(self.run_dir, [self.input], root=self.output)

    def test_files_can_be_looked_up_from_another_folder(self):
        # an input in a subfolder that names files of the folder above it
        scan = self.case / "scan"
        scan.mkdir()
        (scan / "dt.in").write_text('D: DISTRIBUTION, TYPE = FROMFILE, FNAME = "parts.txt";\n')
        run_dir = self.output / "g" / "topic" / "case" / "dt"
        run.prepare(run_dir, [scan / "dt.in"], ref_dir=self.case, root=self.output)
        self.assertEqual((run_dir / "parts.txt").resolve(), (self.case / "parts.txt").resolve())
        self.assertEqual((run_dir / "dt.in").resolve(), (scan / "dt.in").resolve())


class LaunchTest(unittest.TestCase):
    """The codes are replaced by scripts that record how they were called."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.run_dir = self.dir / "run"
        self.run_dir.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def fake(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('#!/bin/sh\npwd -P > called.txt\necho "$@" >> called.txt\n')
        path.chmod(0o755)
        return path

    def test_opalx_runs_in_the_run_folder_with_the_input_name(self):
        exe = self.fake(self.dir / "opalx")
        with mock.patch.dict(os.environ, {"OPALX": str(exe)}):
            rc = run.run_opalx(self.run_dir, "case.in", args=["--restart", "x.h5"])
        self.assertEqual(rc, 0)
        cwd, argv = (self.run_dir / "called.txt").read_text().splitlines()
        self.assertEqual(Path(cwd), self.run_dir.resolve())
        self.assertEqual(argv, "case.in --info 1 --restart x.h5")
        self.assertTrue((self.run_dir / "run.log").exists())

    def test_g4bl_runs_in_the_run_folder(self):
        app = self.dir / "G4beamline.app"
        self.fake(app / "Contents" / "MacOS" / "g4bl")
        with mock.patch.object(paths, "G4BL_APP", app):
            rc = run.run_g4bl(self.run_dir, "case.g4bl")
        self.assertEqual(rc, 0)
        cwd, argv = (self.run_dir / "called.txt").read_text().splitlines()
        self.assertEqual(Path(cwd), self.run_dir.resolve())
        self.assertEqual(argv, "case.g4bl")


class MainTest(unittest.TestCase):
    """The command line: several inputs share one run folder and run in order."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name).resolve()        # /var is a link on macOS
        self.case = root / "studies" / "case"
        self.case.mkdir(parents=True)
        (self.case / "case.g4bl").write_text("beam ascii filename=beam.txt\n")
        (self.case / "beam.txt").write_text("")
        (self.case / "case.in").write_text("OPTION, PSDUMPFREQ = 1;\n")
        self.run_dir = root / "output" / "case"
        # each fake code appends its name, so the order shows in one file
        self.opalx = root / "opalx"
        self.opalx.write_text('#!/bin/sh\necho opalx "$1" >> order.txt\n')
        self.opalx.chmod(0o755)
        self.app = root / "G4beamline.app"
        g4bl = self.app / "Contents" / "MacOS" / "g4bl"
        g4bl.parent.mkdir(parents=True)
        g4bl.write_text('#!/bin/sh\necho g4bl "$1" >> order.txt\n')
        g4bl.chmod(0o755)
        self.patches = [mock.patch.dict(os.environ, {"OPALX": str(self.opalx)}),
                        mock.patch.object(paths, "G4BL_APP", self.app),
                        mock.patch.object(paths, "OUTPUT", root / "output"),
                        mock.patch.object(paths, "STUDIES", root / "studies")]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()

    def test_inputs_run_in_order_in_one_folder(self):
        with contextlib.redirect_stdout(io.StringIO()):
            rc = run.main([str(self.case / "case.g4bl"), str(self.case / "case.in")])
        self.assertEqual(rc, 0)
        self.assertEqual((self.run_dir / "order.txt").read_text().splitlines(),
                         ["g4bl case.g4bl", "opalx case.in"])
        self.assertTrue((self.run_dir / "beam.txt").is_symlink())

    def test_stops_at_the_first_failure(self):
        self.opalx.write_text("#!/bin/sh\nexit 3\n")
        with contextlib.redirect_stdout(io.StringIO()):
            rc = run.main([str(self.case / "case.in"), str(self.case / "case.g4bl")])
        self.assertNotEqual(rc, 0)
        self.assertFalse((self.run_dir / "order.txt").exists())


class PathsCommandTest(unittest.TestCase):
    def test_prints_the_output_folder(self):
        folder = paths.STUDIES / "elements" / "collimator"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            paths.main([str(folder)])
        self.assertEqual(buf.getvalue().strip(), str(paths.output_dir(folder)))


if __name__ == "__main__":
    unittest.main()
