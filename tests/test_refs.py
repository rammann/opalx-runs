"""Tests for opalxruns.refs: which files an input file names, and which are missing.

Run from the repo root:  python -m unittest discover -s tests
"""

import tempfile
import unittest
from pathlib import Path

from opalxruns import refs

OPALX_INPUT = '''// a comment naming FMAPFN = "in_comment.map"
/* a block comment
   FNAME = "in_block_comment.txt"; */
Q1: SOLENOID, L = 1.0, FMAPFN = "maps/sol.map", KS = 1;
Q2: SOLENOID, FMAPFN="/abs/elsewhere.map";
M1: MONITOR, OUTFN = "MON_1";
D1: DUMPEMFIELDS, FILE_NAME = "field.dat";
Dist1: DISTRIBUTION, TYPE = FROMFILE, fname = "parts.txt";
CALL, FILE = "../lattice.in";
'''

LATTICE = '''// called file; its paths are taken from the working folder too
F1: FIELDMAP, FMAPFN = "maps/quad.map";
'''

G4BL_INPUT = '''# fieldmap B file=in_comment.g4blmap
beam ascii filename=beam.txt
fieldmap BY file=by.g4blmap
fieldmap Q file="/abs/quad.g4blmap" current=1
fieldntuple F filename=g4bl_field_x.txt
'''


class OpalxTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.case = root / "study" / "case"
        (self.case / "maps").mkdir(parents=True)
        (self.case / "case.in").write_text(OPALX_INPUT)
        (root / "study" / "lattice.in").write_text(LATTICE)
        (self.case / "maps" / "sol.map").write_text("")
        (self.case / "parts.txt").write_text("")

    def tearDown(self):
        self.tmp.cleanup()

    def test_names_read_files_follows_call_and_skips_comments_and_outputs(self):
        self.assertEqual(refs.references(self.case / "case.in"),
                         ["maps/sol.map", "/abs/elsewhere.map", "parts.txt", "../lattice.in",
                          "maps/quad.map"])

    def test_missing_are_resolved_against_the_working_folder(self):
        self.assertEqual(refs.missing(self.case / "case.in"),
                         ["/abs/elsewhere.map", "maps/quad.map"])

    def test_other_working_folder(self):
        # the same input run from its parent folder finds nothing it names
        self.assertIn("parts.txt", refs.missing(self.case / "case.in", self.case.parent))


class G4blTest(unittest.TestCase):
    def test_reads_fieldmap_and_beam_files_only(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "case.g4bl"
            path.write_text(G4BL_INPUT)
            (Path(d) / "beam.txt").write_text("")
            self.assertEqual(refs.references(path),
                             ["beam.txt", "by.g4blmap", "/abs/quad.g4blmap"])
            self.assertEqual(refs.missing(path), ["by.g4blmap", "/abs/quad.g4blmap"])


if __name__ == "__main__":
    unittest.main()
