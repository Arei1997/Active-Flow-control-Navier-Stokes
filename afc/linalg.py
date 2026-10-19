import numpy as np
from scipy.fft import dctn, idctn
from scipy.linalg import solve_banded


class Tridiagonal:
    """Constant-coefficient tridiagonal system applied along one axis of an array."""

    def __init__(self, diag: np.ndarray, off: float):
        n = len(diag)
        self.n = n
        self.bands = np.zeros((3, n))
        self.bands[0, 1:] = off
        self.bands[1] = diag
        self.bands[2, :-1] = off

    def dense(self) -> np.ndarray:
        return np.diag(self.bands[1]) + np.diag(self.bands[0, 1:], 1) + np.diag(self.bands[2, :-1], -1)

    def solve(self, rhs: np.ndarray, axis: int = 0) -> np.ndarray:
        rhs = np.moveaxis(rhs, axis, 0)
        x = solve_banded((1, 1), self.bands, rhs, check_finite=False)
        return np.moveaxis(x, 0, axis)


class NeumannPoisson:
    """Solves lap(phi) = f with homogeneous Neumann walls on a cell-centred grid.

    The 1D Neumann Laplacian is diagonalised by the DCT-II, which is what the
    eigen-decomposition in assmatFD.m/solveFD.m computes densely. The constant
    null-space mode is removed, giving the zero-mean solution.
    """

    def __init__(self, nx: int, ny: int, hx: float, hy: float):
        lam_x = -4 / hx**2 * np.sin(np.pi * np.arange(nx) / (2 * nx)) ** 2
        lam_y = -4 / hy**2 * np.sin(np.pi * np.arange(ny) / (2 * ny)) ** 2
        denom = lam_x[:, None] + lam_y[None, :]
        denom[0, 0] = 1.0
        self.inv = 1.0 / denom
        self.inv[0, 0] = 0.0

    def solve(self, f: np.ndarray) -> np.ndarray:
        return idctn(dctn(f, type=2, norm="ortho", workers=-1) * self.inv, type=2, norm="ortho", workers=-1)
