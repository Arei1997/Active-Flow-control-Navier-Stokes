"""Implicit viscous step by D'Yakonov approximate factorisation.

Solves (1 - ax Dxx)(1 - ay Dyy) f* = (1 + ax Dxx)(1 + ay Dyy) f^n + dt * rhs
with ax = dt / (2 Re hx^2), as two sweeps of tridiagonal solves.
"""

import numpy as np

from .linalg import Tridiagonal


def _smooth_corners(f: np.ndarray, hx: float, hy: float) -> np.ndarray:
    f = f.copy()
    f[0, 0] = (hx * f[1, 0] + hy * f[0, 1]) / (hx + hy)
    f[-1, 0] = (hx * f[-2, 0] + hy * f[-1, 1]) / (hx + hy)
    f[-1, -1] = (hx * f[-2, -1] + hy * f[-1, -2]) / (hx + hy)
    f[0, -1] = (hx * f[1, -1] + hy * f[0, -2]) / (hx + hy)
    return f


def _second_difference(f: np.ndarray, axis: int) -> np.ndarray:
    if axis == 0:
        return f[2:] - 2 * f[1:-1] + f[:-2]
    return f[:, 2:] - 2 * f[:, 1:-1] + f[:, :-2]


def _interior(f: np.ndarray, axis: int) -> np.ndarray:
    return f[1:-1] if axis == 0 else f[:, 1:-1]


class ViscousSolverU:
    """u: Dirichlet at inlet/outlet (x), zero-gradient ghosts at slip walls (y)."""

    def __init__(self, nxu: int, nyu: int, hx: float, hy: float, dt: float, re: float):
        self.ax = 0.5 * dt / (re * hx**2)
        self.ay = 0.5 * dt / (re * hy**2)
        self.hx, self.hy, self.dt = hx, hy, dt
        self.tx = Tridiagonal(np.full(nxu - 2, 1 + 2 * self.ax), -self.ax)
        dy = np.full(nyu - 2, 1 + 2 * self.ay)
        dy[[0, -1]] = 1 + self.ay
        self.ty = Tridiagonal(dy, -self.ay)

    def solve(self, rhs, u, u_inlet, u_outlet):
        ax, ay = self.ax, self.ay
        u = _smooth_corners(u, self.hx, self.hy)
        g = ay * _second_difference(u, 1) + _interior(u, 1)
        g = ax * _second_difference(g, 0) + _interior(g, 0) + self.dt * rhs

        g[0, :] += ax * (_interior(u_inlet, 0) - ay * _second_difference(u_inlet, 0))
        g[-1, :] += ax * (_interior(u_outlet, 0) - ay * _second_difference(u_outlet, 0))
        g = self.tx.solve(g, axis=0)
        return self.ty.solve(g, axis=1)


class ViscousSolverV:
    """v: reflected ghosts at inlet/outlet (x), Dirichlet at slip walls (y)."""

    def __init__(self, nxv: int, nyv: int, hx: float, hy: float, dt: float, re: float):
        self.ax = 0.5 * dt / (re * hx**2)
        self.ay = 0.5 * dt / (re * hy**2)
        self.hx, self.hy, self.dt = hx, hy, dt
        dx = np.full(nxv - 2, 1 + 2 * self.ax)
        dx[[0, -1]] = 1 + 3 * self.ax
        self.tx = Tridiagonal(dx, -self.ax)
        self.ty = Tridiagonal(np.full(nyv - 2, 1 + 2 * self.ay), -self.ay)

    def solve(self, rhs, v, v_south, v_north, v_inlet, v_outlet):
        ax, ay = self.ax, self.ay
        v = _smooth_corners(v, self.hx, self.hy)
        g = ax * _second_difference(v, 0) + _interior(v, 0)
        g = ay * _second_difference(g, 1) + _interior(g, 1) + self.dt * rhs

        g[:, 0] += ay * (_interior(v_south, 0) - ax * _second_difference(v_south, 0))
        g[:, -1] += ay * (_interior(v_north, 0) - ax * _second_difference(v_north, 0))
        g = self.ty.solve(g, axis=1)

        g[0, :] += 2 * ax * _interior(v_inlet, 0)
        g[-1, :] += 2 * ax * _interior(v_outlet, 0)
        return self.tx.solve(g, axis=0)
