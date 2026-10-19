"""Direct-forcing immersed boundary (Pinelli et al. style) with optional RKPM.

Interpolation and spreading are stored as one sparse matrix of discrete delta
weights, W[l, k] = delta_h(x_k - X_l), so that
    interpolate(f) = (W * area) @ f
    spread(F)      = W.T @ (F * epsilon * dS)
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .grid import StaggeredGrid

_WINDOW = np.arange(-5, 6)
_BAND = 5


@dataclass
class Markers:
    x: np.ndarray
    y: np.ndarray
    ds: np.ndarray
    theta: np.ndarray
    n_circle: int
    spacing: float

    @property
    def n(self) -> int:
        return len(self.x)

    @property
    def tail(self) -> slice:
        return slice(self.n_circle, self.n)

    @property
    def hinge(self) -> tuple[float, float]:
        return self.x[self.n_circle - 1], self.y[self.n_circle - 1]

    def with_tail_angle(self, angle: float) -> "Markers":
        k = np.arange(1, self.n - self.n_circle + 1)
        hx, hy = self.hinge
        x, y = self.x.copy(), self.y.copy()
        x[self.tail] = hx + k * self.spacing * np.cos(angle)
        y[self.tail] = hy + k * self.spacing * np.sin(angle)
        return Markers(x, y, self.ds, self.theta, self.n_circle, self.spacing)


def circle_markers(n: int, xc: float, yc: float, rc: float):
    theta = 2 * np.pi / n * np.arange(1, n + 1)
    th_next = np.append(theta[1:], theta[0] + 2 * np.pi)
    th_prev = np.insert(theta[:-1], 0, theta[-1] - 2 * np.pi)
    ds = rc * (0.5 * (th_next + theta) - 0.5 * (theta + th_prev))
    return xc + rc * np.cos(theta), yc + rc * np.sin(theta), ds, theta


def build_markers(grid: StaggeredGrid, xc, yc, rc, angle=0.0, tail=True, n_circle=None) -> Markers:
    """Circle with an optional rigid tail of length ~D hinged at the rear point."""
    spacing = np.hypot(np.diff(grid.xu).min(), np.diff(grid.yu).min())
    if n_circle is None:
        n_circle = int(np.floor(0.95 * 2 * np.pi * rc / spacing))
    x, y, ds, theta = circle_markers(n_circle, xc, yc, rc)
    markers = Markers(x, y, ds, theta, n_circle, spacing)
    if not tail:
        return markers

    n_tail = int(np.ceil(2 * rc / spacing))
    markers = Markers(
        np.append(x, np.zeros(n_tail)),
        np.append(y, np.zeros(n_tail)),
        np.append(ds, np.full(n_tail, spacing)),
        np.append(theta, np.full(n_tail, 2 * np.pi)),
        n_circle,
        spacing,
    )
    return markers.with_tail_angle(angle)


def regular_delta(r: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Roma, Peskin & Berger three-point discrete delta."""
    a = np.abs(r) / h
    near = (1 + np.sqrt(np.clip(1 - 3 * a**2, 0, None))) / (3 * h)
    far = (5 - 3 * a - np.sqrt(np.clip(1 - 3 * (1 - a) ** 2, 0, None))) / (6 * h)
    return np.where(a > 1.5, 0.0, np.where(a > 0.5, far, near))


class ImmersedBoundary:
    """Transfer operators between one staggered Eulerian grid and the markers."""

    def __init__(self, xx: np.ndarray, yy: np.ndarray, markers: Markers, rkpm: bool = False):
        self.shape = (len(xx), len(yy))
        self.markers = markers
        n = markers.n

        ii = np.abs(xx[None, :] - markers.x[:, None]).argmin(axis=1)
        jj = np.abs(yy[None, :] - markers.y[:, None]).argmin(axis=1)
        hk_x = self._kernel_width(xx, ii, markers.x)
        hk_y = self._kernel_width(yy, jj, markers.y)

        I = ii[:, None] + _WINDOW
        J = jj[:, None] + _WINDOW
        rx = xx[I] - markers.x[:, None]
        ry = yy[J] - markers.y[:, None]
        dx = regular_delta(rx, hk_x[:, None])
        dy = regular_delta(ry, hk_y[:, None])
        inside = (np.abs(rx) <= 1.5 * hk_x[:, None])[:, :, None] & (np.abs(ry) <= 1.5 * hk_y[:, None])[:, None, :]

        cell_x = 0.5 * (xx[I + 1] - xx[I - 1])
        cell_y = 0.5 * (yy[J + 1] - yy[J - 1])
        area = cell_x[:, :, None] * cell_y[:, None, :]
        delta = dx[:, :, None] * dy[:, None, :]

        rxx = np.broadcast_to(rx[:, :, None], delta.shape)
        ryy = np.broadcast_to(ry[:, None, :], delta.shape)
        basis = np.stack([np.ones_like(rxx), rxx, ryy, rxx * ryy, rxx**2, ryy**2], axis=-1)
        if rkpm:
            self.b = self._rkpm_coefficients(basis, delta, area, inside, hk_x, hk_y)
        else:
            self.b = np.zeros((n, 6))
            self.b[:, 0] = 1.0
        weights = delta * np.einsum("nijk,nk->nij", basis, self.b)

        rows = np.broadcast_to(np.arange(n)[:, None, None], delta.shape)[inside]
        cols = (I[:, :, None] * self.shape[1] + J[:, None, :])[inside]
        size = (n, self.shape[0] * self.shape[1])
        self.W = sp.csr_matrix((weights[inside], (rows, cols)), shape=size)
        self.W_area = sp.csr_matrix(((weights * area)[inside], (rows, cols)), shape=size)
        self.WT = self.W.T.tocsr()

        if rkpm:
            self.epsilon = self._rkpm_epsilon()
        else:
            self.epsilon = np.full(n, np.hypot(xx[1] - xx[0], yy[1] - yy[0]) / np.sqrt(2))

    @staticmethod
    def _kernel_width(x, idx, xt):
        hp = np.maximum(x[idx + 1] - x[idx], x[idx] - x[idx - 1])
        hm = np.minimum(x[idx + 1] - x[idx], x[idx] - x[idx - 1])
        return 5 / 6 * hp + 1 / 6 * hm + 1 / 9 * np.abs(xt - x[idx])

    @staticmethod
    def _rkpm_coefficients(basis, delta, area, inside, hk_x, hk_y):
        """Polynomial correction so the kernel reproduces moments up to 2nd order."""
        scale = np.stack(
            [np.ones_like(hk_x), 1 / hk_x, 1 / hk_y, 1 / (hk_x * hk_y), 1 / hk_x**2, 1 / hk_y**2], axis=-1
        )
        phi = basis * scale[:, None, None, :]
        w = np.where(inside, area * delta, 0.0)
        moments = np.einsum("nij,nija,nijb->nab", w, phi, phi)
        b = np.zeros((len(hk_x), 6))
        b[:, 0] = 1.0
        rhs = np.zeros(6)
        rhs[0] = 1.0
        for m in np.flatnonzero(np.linalg.det(moments) != 0):
            b[m] = np.linalg.solve(moments[m], rhs)
        return b * scale

    def _rkpm_epsilon(self):
        """Width of the forcing shell from the banded self-interaction system."""
        n = self.markers.n
        gram = (self.W_area @ self.WT).toarray()
        lag = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        band = np.minimum(lag, n - lag) <= _BAND
        system = np.where(band, gram * self.markers.ds[None, :], 0.0)
        return np.linalg.solve(system, np.ones(n))

    def interpolate(self, field: np.ndarray) -> np.ndarray:
        return self.W_area @ field.ravel()

    def spread(self, values: np.ndarray) -> np.ndarray:
        return (self.WT @ (values * self.epsilon * self.markers.ds)).reshape(self.shape)

    def forcing(self, rhs: np.ndarray, target: np.ndarray | None = None):
        """Eulerian force that drives the predicted velocity towards `target` on the markers.

        `rhs` is minus the explicit predictor over dt, so interpolating it gives
        -u~/dt at the markers and the required Lagrangian force is (U_wall - u~)/dt.
        """
        values = self.interpolate(rhs)
        if target is not None:
            values = values + target
        return self.spread(values), values
