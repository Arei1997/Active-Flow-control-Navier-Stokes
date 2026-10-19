import numpy as np
import pytest

from afc.linalg import NeumannPoisson, Tridiagonal
from afc.operators import laplacian_neumann
from afc.postprocess import streamfunction
from afc.viscous import ViscousSolverU, ViscousSolverV

from matlab_reference import dyakonov_u, dyakonov_v, poisson_eig

rng = np.random.default_rng(0)


def test_tridiagonal_matches_dense():
    diag = 2 + rng.random(9)
    tri = Tridiagonal(diag, -0.4)
    rhs = rng.random((9, 5))
    np.testing.assert_allclose(tri.solve(rhs, axis=0), np.linalg.solve(tri.dense(), rhs))
    np.testing.assert_allclose(tri.solve(rhs.T, axis=1), np.linalg.solve(tri.dense(), rhs).T)


def test_poisson_matches_matlab_eigen_solver():
    hx, hy = 0.3, 0.2
    f = rng.standard_normal((14, 11))
    f -= f.mean()
    phi = NeumannPoisson(14, 11, hx, hy).solve(f)
    np.testing.assert_allclose(phi, poisson_eig(f, hx, hy), atol=1e-12)
    np.testing.assert_allclose(laplacian_neumann(phi, hx, hy), f, atol=1e-10)
    assert abs(phi.mean()) < 1e-12


@pytest.mark.parametrize("component", ["u", "v"])
def test_dyakonov_matches_matlab_loops(component):
    hx, hy, dt, re = 0.11, 0.09, 0.01, 40.0
    ncx, ncy = 13, 10
    ub = rng.random((4, 16))
    if component == "u":
        nx, ny = ncx + 1, ncy + 2
        un, f, ubs = rng.random((nx, ny)), rng.random((nx - 2, ny - 2)), rng.random(ny)
        ours = ViscousSolverU(nx, ny, hx, hy, dt, re).solve(f, un, ub[3, :ny], ubs)
        ref = dyakonov_u(f, un, ub, ubs, hx, hy, re, dt)
    else:
        nx, ny = ncx + 2, ncy + 1
        un, f, vbs = rng.random((nx, ny)), rng.random((nx - 2, ny - 2)), rng.random(ny)
        ours = ViscousSolverV(nx, ny, hx, hy, dt, re).solve(f, un, ub[0, :nx], ub[2, :nx], ub[3, :ny], vbs)
        ref = dyakonov_v(f, un, ub, vbs, hx, hy, re, dt)
    np.testing.assert_allclose(ours, ref, rtol=1e-12, atol=1e-13)


def test_streamfunction_inverts_laplacian():
    hx, hy = 0.2, 0.15
    omega = rng.standard_normal((12, 9))
    psi = streamfunction(omega, hx, hy)
    lap = (psi[2:, 1:-1] - 2 * psi[1:-1, 1:-1] + psi[:-2, 1:-1]) / hx**2 + (
        psi[1:-1, 2:] - 2 * psi[1:-1, 1:-1] + psi[1:-1, :-2]
    ) / hy**2
    np.testing.assert_allclose(lap, -omega[1:-1, 1:-1], atol=1e-10)
    assert np.all(psi[[0, -1], :] == 0) and np.all(psi[:, [0, -1]] == 0)
