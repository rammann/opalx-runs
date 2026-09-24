"""Tests for opalxruns.results.

Run from the repo root:  python -m unittest discover -s tests
"""

import contextlib
import io
import unittest

from opalxruns.results import Results


def printed(results, title=""):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        results.print_table(title)
    return buf.getvalue().splitlines()


class CheckTest(unittest.TestCase):
    def test_absolute_tolerance(self):
        r = Results()
        self.assertTrue(r.check("1", "inside", 1.05, 1.0, 0.1))
        self.assertFalse(r.check("1", "outside", 1.2, 1.0, 0.1))
        self.assertEqual(r.failed, 1)

    def test_relative_tolerance(self):
        r = Results()
        self.assertTrue(r.check("1", "rel", 101.0, 100.0, 0.02, rel=True))
        self.assertFalse(r.check("1", "rel", 103.0, 100.0, 0.02, rel=True))

    def test_no_tolerance_is_a_diagnostic_that_never_fails(self):
        r = Results()
        self.assertIsNone(r.check("1", "diag", 5.0, 0.0, None))
        self.assertEqual(r.failed, 0)

    def test_at_least_fails_below_the_minimum(self):
        r = Results()
        self.assertTrue(r.at_least("1", "big enough", 2.0, 1.0))
        self.assertFalse(r.at_least("1", "too small", 0.5, 1.0))
        self.assertEqual(r.failed, 1)

    def test_note_adds_a_row_without_a_result(self):
        r = Results()
        r.note("1", "remark", "text")
        self.assertIsNone(r.rows[0]["ok"])


class LayoutTest(unittest.TestCase):
    def rows(self, **layout):
        r = Results(**layout)
        r.check("7", "quantity", 1.5, 1.25, 0.5)
        r.note("7", "remark", "just text")
        return printed(r, "Title")

    def test_default_layout(self):
        out = self.rows()          # "", rule, title, rule, header, dashes, rows...
        self.assertEqual(out[1], "=" * 104)
        self.assertEqual(out[3], "=" * 104)
        self.assertEqual(out[4], f"{'test':4s} {'quantity':34s} {'measured':>13s} "
                                 f"{'expected':>13s} {'|diff|':>10s} {'tol':>9s}  result")
        self.assertEqual(out[5], "-" * len(out[4]))
        self.assertEqual(out[6], f"{'7':>4} {'quantity':34s} {1.5:>13.6g} {1.25:>13.6g} "
                                 f"{0.25:>10.3g} {'5.0e-01':>9s}  PASS")
        self.assertEqual(out[7], f"{'7':>4} {'remark':34s} just text")

    def test_wide_left_aligned_case_column(self):
        out = self.rows(first="case", first_width=16, first_align="<", rule=110)
        self.assertEqual(out[1], "=" * 110)
        self.assertTrue(out[4].startswith(f"{'case':16s} {'quantity':34s} "))
        self.assertTrue(out[6].startswith(f"{'7':16s} {'quantity':34s} "))

    def test_analytic_column_and_narrow_names(self):
        out = self.rows(name_width=30, expected="analytic", rule=100)
        self.assertEqual(out[1], "=" * 100)
        self.assertTrue(out[4].startswith(f"{'test':4s} {'quantity':30s} {'measured':>13s} "
                                          f"{'analytic':>13s} "))
        self.assertTrue(out[6].startswith(f"{'7':>4} {'quantity':30s} "))


if __name__ == "__main__":
    unittest.main()
