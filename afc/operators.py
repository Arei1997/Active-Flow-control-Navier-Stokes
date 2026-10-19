"""Finite-difference operators on the staggered grid."""

import numpy as np


def laplacian(f: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """Five-point Laplacian on the interior of f."""
    return (
        (f[2:, 1:-1] + f[:-2, 1:-1]) / hx**2
        + (f[1:-1, 2:] + f[1:-1, :-2]) / hy**2
        - (2 / hx**2 + 2 / hy**2) * f[1:-1, 1:-1]
    )


def laplacian_neumann(phi: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """Laplacian of a cell-centred field padded with zero-gradient ghosts."""
    return laplacian(np.pad(phi, 1, mode="edge"), hx, hy)


def advection_u(u: np.ndarray, v: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """d(uu)/dx + d(uv)/dy at interior u nodes."""
    nxu = u.shape[0]
    uc = u[1:-1, 1:-1]
    flux_n = (v[2:nxu, 1:] + v[1 : nxu - 1, 1:]) * (u[1:-1, 2:] + uc)
    flux_s = (v[2:nxu, :-1] + v[1 : nxu - 1, :-1]) * (uc + u[1:-1, :-2])
    return 0.25 / hy * (flux_n - flux_s) + (u[2:, 1:-1] ** 2 - u[:-2, 1:-1] ** 2) / (2 * hx)


def advection_v(u: np.ndarray, v: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """d(uv)/dx + d(vv)/dy at interior v nodes."""
    nyv = v.shape[1]
    vc = v[1:-1, 1:-1]
    flux_e = (u[1:, 2:nyv] + u[1:, 1 : nyv - 1]) * (vc + v[2:, 1:-1])
    flux_w = (u[:-1, 2:nyv] + u[:-1, 1 : nyv - 1]) * (vc + v[:-2, 1:-1])
    return 0.25 / hx * (flux_e - flux_w) + (v[1:-1, 2:] ** 2 - v[1:-1, :-2] ** 2) / (2 * hy)


def pressure_gradient(p: np.ndarray, hx: float, hy: float, axis: int) -> np.ndarray:
    h = hx if axis == 0 else hy
    return np.diff(p, axis=axis) / h


def divergence(u: np.ndarray, v: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """Cell-centred divergence, shape (ncx, ncy)."""
    return np.diff(u[:, 1:-1], axis=0) / hx + np.diff(v[1:-1, :], axis=1) / hy


def explicit_rhs(f, lap, nl, lap_old, nl_old, grad_p, dt, re, ab1, ab2):
    """Minus the explicit predictor (f* - f^n)/dt - f^n/dt, as in sumupAB.m."""
    return -f[1:-1, 1:-1] / dt + ab1 * (nl - lap / re) + ab2 * (lap_old / re - nl_old) + grad_p


def convective_outflow_u(u: np.ndarray, dt: float, hx: float, q_in: float) -> np.ndarray:
    """Convective outlet for u, rescaled to conserve the inlet flow rate."""
    uc = u[-2, 1:-1].mean()
    ub = u[-1, :] - dt / hx * uc * (u[-1, :] - u[-2, :])
    return ub * q_in / ub[1:-1].sum()


def convective_outflow_v(u: np.ndarray, v: np.ndarray, dt: float, hx: float) -> np.ndarray:
    uc = u[-2, 1:-1].mean()
    return v[-2, :] - dt / hx * uc * (v[-2, :] - v[-3, :])


def apply_boundary_conditions(us, vs, ub, vb, ub_out, vb_out):
    """Embed predicted interiors and fill boundary/ghost values (setbc.m)."""
    nxu, nyu = us.shape[0] + 2, us.shape[1] + 2
    nxv, nyv = vs.shape[0] + 2, vs.shape[1] + 2
    u = np.zeros((nxu, nyu))
    v = np.zeros((nxv, nyv))
    u[1:-1, 1:-1] = us
    v[1:-1, 1:-1] = vs

    u[0, :] = ub[3, :nyu]
    u[-1, :] = ub_out
    u[1:, 0] = u[1:, 1]
    u[1:, -1] = u[1:, -2]

    v[:, 0] = vb[0, :nxv]
    v[:, -1] = vb[2, :nxv]
    v[0, 1:-1] = 2 * vb[3, 1 : nyv - 1] - v[1, 1:-1]
    v[-1, 1:-1] = 2 * vb_out[1:-1] - v[-2, 1:-1]
    return u, v


def project(phi, u, v, p, dt, re, hx, hy):
    """Velocity correction and rotational pressure update (finalvelocityVK.m)."""
    uf = u.copy()
    vf = v.copy()
    uf[1:-1, :] = u[1:-1, :] - dt / hx * np.diff(np.pad(phi, ((0, 0), (1, 1)), mode="edge"), axis=0)
    vf[:, 1:-1] = v[:, 1:-1] - dt / hy * np.diff(np.pad(phi, ((1, 1), (0, 0)), mode="edge"), axis=1)
    p_new = p + phi - 0.5 * dt / re * laplacian_neumann(phi, hx, hy)
    return uf, vf, p_new
