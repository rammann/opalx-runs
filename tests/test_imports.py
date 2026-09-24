"""Every module of the package imports, and the old flat module names are gone.

Run from the repo root:  python -m unittest discover -s tests
"""

import importlib
import sys
import unittest

MODULES = [
    "paths",
    "opalx_diagnostics", "opalx_run", "plot_stat", "plot_timing", "plot_monitors",
    "process_run", "particles_to_vtk", "elements_to_vtk",
    "g4bl", "h5", "mue4", "mue4_beam",
    "case", "matrices", "particles", "plotstyle", "refs", "results",
]


class ImportTest(unittest.TestCase):
    def test_every_module_imports(self):
        for name in MODULES:
            with self.subTest(module=name):
                importlib.import_module(f"opalxruns.{name}")

    def test_no_flat_module_names_leak_into_sys_modules(self):
        # The old processing/ scripts put their own folder on sys.path and
        # imported each other by bare name; the package must not do that.
        for name in MODULES:
            importlib.import_module(f"opalxruns.{name}")
        for bare in ("opalx_run", "opalx_diagnostics", "plot_stat", "mue4lib"):
            with self.subTest(module=bare):
                self.assertNotIn(bare, sys.modules)


if __name__ == "__main__":
    unittest.main()
