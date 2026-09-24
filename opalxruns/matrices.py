"""Transfer matrices in (x, x', y, y', z, delta): drifts, the symplectic check,
and the matrix measured from tracked +/- particle pairs."""

from __future__ import annotations

import numpy as np

# The +/- pairs of the 13-particle set, as column indices: column 0 is the
# reference particle, then one (plus, minus) pair per coordinate.
PAIRS = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12)]


def drift_matrix(L: float, gamma: float) -> np.ndarray:
    M = np.eye(6)
    M[0, 1] = L
    M[2, 3] = L
    M[4, 5] = L / gamma**2
    return M


def symplectic_residual(M: np.ndarray) -> float:
    """max |M^T J M - J| for a 2n x 2n matrix, J being the usual block form with
    one [[0, 1], [-1, 0]] block per coordinate pair."""
    n = M.shape[0]
    J = np.zeros((n, n))
    for i in range(0, n, 2):
        J[i, i + 1] = 1.0
        J[i + 1, i] = -1.0
    return float(np.abs(M.T @ J @ M - J).max())


def centred_differences(X: np.ndarray, pairs=PAIRS) -> np.ndarray:
    """6 x 6: column j is (X[:, plus] - X[:, minus]) / 2 for the j-th pair."""
    d = np.zeros((6, 6))
    for j, (ip, im) in enumerate(pairs):
        d[:, j] = (X[:, ip] - X[:, im]) / 2.0
    return d


def transfer_matrix(Xin: np.ndarray, Xout: np.ndarray, pairs=PAIRS) -> np.ndarray:
    """6 x 6 by centred differences over the +/- pairs, entrance plane to exit.

    Built as (d out) @ inv(d in) rather than assuming the entrance differences
    are exactly the nominal steps, so a plane that sits a little downstream of
    where the particles were made does not bias it.
    """
    return centred_differences(Xout, pairs) @ np.linalg.inv(centred_differences(Xin, pairs))
