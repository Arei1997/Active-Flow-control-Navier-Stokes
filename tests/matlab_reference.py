"""Literal loop-by-loop translations of MATLAB routines, used as test oracles."""

import numpy as np


def trid_lu(a, e, c):
    n = len(a)
    alfa, beta = np.zeros(n), np.zeros(n)
    alfa[0] = a[0]
    for i in range(1, n):
        beta[i] = e[i] / alfa[i - 1]
        alfa[i] = a[i] - beta[i] * c[i - 1]
    return alfa, beta


def solve_lu(alfa, beta, c, f):
    n = len(f)
    y, u = np.zeros(n), np.zeros(n)
    y[0] = f[0]
    for i in range(1, n):
        y[i] = f[i] - beta[i] * y[i - 1]
    u[-1] = y[-1] / alfa[-1]
    for i in range(n - 2, -1, -1):
        u[i] = (y[i] - c[i] * u[i + 1]) / alfa[i]
    return u


def _factors(n, alpha, ends):
    a = np.full(n, 1 + 2 * alpha)
    a[[0, -1]] = ends
    off = np.full(n, -alpha)
    alfa, beta = trid_lu(a, off, off)
    return alfa, beta, off


def dyakonov_u(f, un, ub, ubs, hx, hy, re, dt):
    m, n = f.shape
    nx, ny = m + 2, n + 2
    ax, ay = 0.5 * dt / (re * hx**2), 0.5 * dt / (re * hy**2)
    lx = _factors(nx - 2, ax, 1 + 2 * ax)
    ly = _factors(ny - 2, ay, 1 + ay)
    un = un.copy()
    un[0, 0] = (hx * un[1, 0] + hy * un[0, 1]) / (hx + hy)
    un[-1, 0] = (hx * un[-2, 0] + hy * un[-1, 1]) / (hx + hy)
    un[-1, -1] = (hx * un[-2, -1] + hy * un[-1, -2]) / (hx + hy)
    un[0, -1] = (hx * un[1, -1] + hy * un[0, -2]) / (hx + hy)
    g1 = np.zeros((nx, ny - 2))
    for i in range(nx):
        for j in range(1, ny - 1):
            g1[i, j - 1] = ay * (un[i, j + 1] - 2 * un[i, j] + un[i, j - 1]) + un[i, j]
    g = np.zeros((nx - 2, ny - 2))
    for i in range(1, nx - 1):
        g[i - 1, :] = ax * (g1[i + 1] - 2 * g1[i] + g1[i - 1]) + g1[i]
    g += dt * f
    for j in range(1, ny - 1):
        left = ub[3, j] - ay * (ub[3, j + 1] - 2 * ub[3, j] + ub[3, j - 1])
        right = ubs[j] - ay * (ubs[j + 1] - 2 * ubs[j] + ubs[j - 1])
        gg = g[:, j - 1].copy()
        gg[0] += ax * left
        gg[-1] += ax * right
        g[:, j - 1] = solve_lu(*lx[:2], lx[2], gg)
    us = np.zeros_like(g)
    for i in range(nx - 2):
        us[i, :] = solve_lu(*ly[:2], ly[2], g[i, :])
    return us


def dyakonov_v(f, un, vb, vbs, hx, hy, re, dt):
    m, n = f.shape
    nx, ny = m + 2, n + 2
    ax, ay = 0.5 * dt / (re * hx**2), 0.5 * dt / (re * hy**2)
    lx = _factors(nx - 2, ax, 1 + 3 * ax)
    ly = _factors(ny - 2, ay, 1 + 2 * ay)
    un = un.copy()
    un[0, 0] = (hx * un[1, 0] + hy * un[0, 1]) / (hx + hy)
    un[-1, 0] = (hx * un[-2, 0] + hy * un[-1, 1]) / (hx + hy)
    un[-1, -1] = (hx * un[-2, -1] + hy * un[-1, -2]) / (hx + hy)
    un[0, -1] = (hx * un[1, -1] + hy * un[0, -2]) / (hx + hy)
    g1 = ax * (un[2:] - 2 * un[1:-1] + un[:-2]) + un[1:-1]
    g = ay * (g1[:, 2:] - 2 * g1[:, 1:-1] + g1[:, :-2]) + g1[:, 1:-1]
    g += dt * f
    for i in range(1, nx - 1):
        bot = vb[0, i] - ax * vb[0, i - 1] + 2 * ax * vb[0, i] - ax * vb[0, i + 1]
        top = vb[2, i] - ax * vb[2, i - 1] + 2 * ax * vb[2, i] - ax * vb[2, i + 1]
        gg = g[i - 1, :].copy()
        gg[0] += ay * bot
        gg[-1] += ay * top
        g[i - 1, :] = solve_lu(*ly[:2], ly[2], gg)
    us = np.zeros_like(g)
    for j in range(ny - 2):
        gg = g[:, j].copy()
        gg[0] += 2 * ax * vb[3, j + 1]
        gg[-1] += 2 * ax * vbs[j + 1]
        us[:, j] = solve_lu(*lx[:2], lx[2], gg)
    return us


def poisson_eig(f, hx, hy):
    """assmatFD.m + solveFD.m: dense eigen-decomposition of the Neumann Laplacian."""
    def neumann(n, h):
        a = (np.diag(np.full(n, -2.0)) + np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)) / h**2
        a[0, 0] = a[-1, -1] = -1 / h**2
        return np.linalg.eigh(a)

    ex, px = neumann(f.shape[0], hx)
    ey, py = neumann(f.shape[1], hy)
    emix = ex[:, None] + ey[None, :]
    zero = np.unravel_index(np.abs(emix).argmin(), emix.shape)
    emix[zero] = 1.0
    g = px.T @ f @ py
    g[zero] = 0.0
    return px @ (g / emix) @ py.T
