"""Tests for opalxruns.matrices.

Run from the repo root:  python -m unittest discover -s tests
"""

import unittest

import numpy as np

from opalxruns import matrices


class DriftTest(unittest.TestCase):
    def test_drift_entries(self):
        M = matrices.drift_matrix(2.0, gamma=4.0)
        expected = np.eye(6)
        expected[0, 1] = expected[2, 3] = 2.0
        expected[4, 5] = 2.0 / 16.0
        np.testing.assert_array_equal(M, expected)


class SymplecticTest(unittest.TestCase):
    def test_identity_and_drift_are_symplectic_in_4d_and_6d(self):
        self.assertEqual(matrices.symplectic_residual(np.eye(4)), 0.0)
        self.assertEqual(matrices.symplectic_residual(np.eye(6)), 0.0)
        self.assertEqual(matrices.symplectic_residual(matrices.drift_matrix(1.0, 1.0)), 0.0)

    def test_scaled_identity_is_not(self):
        # (2I)^T J (2I) - J = 3 J, so the largest entry is 3.
        self.assertEqual(matrices.symplectic_residual(2.0 * np.eye(6)), 3.0)


class TransferMatrixTest(unittest.TestCase):
    def test_recovers_a_linear_map_from_plus_minus_pairs(self):
        rng = np.random.default_rng(1)
        M = rng.normal(size=(6, 6))
        eps = 1e-3
        Xin = np.zeros((6, 13))
        for j, (ip, im) in enumerate(matrices.PAIRS):
            Xin[j, ip], Xin[j, im] = +eps, -eps
        Xout = M @ Xin
        np.testing.assert_allclose(matrices.transfer_matrix(Xin, Xout), M, rtol=1e-12)

    def test_centred_differences_are_half_the_pair_difference(self):
        X = np.arange(6 * 13, dtype=float).reshape(6, 13)
        d = matrices.centred_differences(X, matrices.PAIRS)
        for j, (ip, im) in enumerate(matrices.PAIRS):
            np.testing.assert_array_equal(d[:, j], (X[:, ip] - X[:, im]) / 2.0)


if __name__ == "__main__":
    unittest.main()
