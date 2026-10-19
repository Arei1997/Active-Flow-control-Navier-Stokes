import numpy as np
import pytest

from afc.grid import StaggeredGrid
from afc.immersed import ImmersedBoundary, build_markers, regular_delta


@pytest.fixture
def grid():
    return StaggeredGrid(120, 100, 6.0, 5.0)


def test_regular_delta_moments():
    h = 0.1
    for shift in (0.0, 0.013, 0.05, 0.077):
        x = (np.arange(-4, 5) + shift / h) * h
        w = regular_delta(x, h) * h
        assert w.sum() == pytest.approx(1.0)
        assert (w * x).sum() == pytest.approx(0.0, abs=1e-14)


def test_markers_geometry(grid):
    m = build_markers(grid, 2.0, 2.5, 0.5, angle=np.pi / 4)
    r = np.hypot(m.x[: m.n_circle] - 2.0, m.y[: m.n_circle] - 2.5)
    np.testing.assert_allclose(r, 0.5)
    assert m.ds[: m.n_circle].sum() == pytest.approx(np.pi)
    assert m.hinge == pytest.approx((2.5, 2.5))
    tip = np.hypot(m.x[-1] - 2.5, m.y[-1] - 2.5)
    assert tip == pytest.approx((m.n - m.n_circle) * m.spacing)
    np.testing.assert_allclose(np.arctan2(m.y[m.tail] - 2.5, m.x[m.tail] - 2.5), np.pi / 4)


@pytest.mark.parametrize("rkpm", [False, True])
def test_interpolation_reproduces_polynomials(grid, rkpm):
    m = build_markers(grid, 2.0, 2.5, 0.5, angle=0.3)
    xx, yy = grid.xu[1:-1], grid.yu[1:-1]
    ib = ImmersedBoundary(xx, yy, m, rkpm=rkpm)
    X, Y = np.meshgrid(xx, yy, indexing="ij")
    # without RKPM the kernel is widened to h + |dx|/9, so moments hold only approximately
    tol = 1e-10 if rkpm else 2e-2
    np.testing.assert_allclose(ib.interpolate(np.ones_like(X)), 1.0, atol=tol)
    np.testing.assert_allclose(ib.interpolate(1 + 2 * X - 3 * Y), 1 + 2 * m.x - 3 * m.y, atol=10 * tol)
    if rkpm:
        np.testing.assert_allclose(ib.interpolate(X * Y + X**2), m.x * m.y + m.x**2, atol=1e-10)


def test_spread_is_adjoint_of_interpolation(grid):
    m = build_markers(grid, 2.0, 2.5, 0.5, angle=-0.4)
    ib = ImmersedBoundary(grid.xv[1:-1], grid.yv[1:-1], m)
    field = np.random.default_rng(1).random(ib.shape)
    lag = np.random.default_rng(2).random(m.n)
    lhs = (ib.spread(lag) * field).sum() * grid.hx * grid.hy
    rhs = (ib.interpolate(field) * lag * ib.epsilon * m.ds).sum()
    assert lhs == pytest.approx(rhs)
