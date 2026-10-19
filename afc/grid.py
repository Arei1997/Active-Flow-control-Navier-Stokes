"""Staggered MAC grid.

Pressure lives at cell centres (ncx, ncy). u lives on vertical faces with one
ghost row above and below (ncx+1, ncy+2); v lives on horizontal faces with one
ghost column left and right (ncx+2, ncy+1). Arrays are indexed [i, j] with i
along x, matching the MATLAB layout.

Boundary arrays follow the thesis convention, with rows ordered
    0: south, 1: east, 2: north, 3: west.
"""

from dataclasses import dataclass

import numpy as np

from .config import SimulationConfig


@dataclass
class StaggeredGrid:
    ncx: int
    ncy: int
    lx: float
    ly: float

    def __post_init__(self):
        self.hx = self.lx / self.ncx
        self.hy = self.ly / self.ncy
        self.nxu, self.nyu = self.ncx + 1, self.ncy + 2
        self.nxv, self.nyv = self.ncx + 2, self.ncy + 1
        self.xu = np.arange(self.nxu) * self.hx
        self.yu = -self.hy / 2 + np.arange(self.nyu) * self.hy
        self.xv = -self.hx / 2 + np.arange(self.nxv) * self.hx
        self.yv = np.arange(self.nyv) * self.hy
        self.xp = self.xv[1:-1]
        self.yp = self.yu[1:-1]

    @classmethod
    def from_config(cls, config: SimulationConfig) -> "StaggeredGrid":
        return cls(config.ncx, config.ncy, config.lx, config.ly)


def boundary_velocities(grid: StaggeredGrid, inflow: float = 1.0):
    """Uniform stream: u = inflow on every edge, v = 0 (set_velocity_bc.m)."""
    ub = np.zeros((4, max(grid.nxu, grid.nyu)))
    ub[[0, 2], : grid.nxu] = inflow
    ub[[1, 3], : grid.nyu] = inflow
    vb = np.zeros((4, max(grid.nxv, grid.nyv)))
    return ub, vb


def initial_fields(grid: StaggeredGrid, ub: np.ndarray, vb: np.ndarray):
    """Impulsive start from the uniform stream (initialize.m)."""
    u = np.ones((grid.nxu, grid.nyu))
    v = np.zeros((grid.nxv, grid.nyv))
    p = np.zeros((grid.ncx, grid.ncy))

    u[0, :] = ub[3, : grid.nyu]
    u[-1, :] = ub[3, : grid.nyu]
    v[0, :] = 2 * vb[3, : grid.nyv] - v[1, :]
    v[-1, :] = 2 * vb[1, : grid.nyv] - v[-2, :]
    v[:, 0] = vb[0, : grid.nxv]
    v[:, -1] = vb[2, : grid.nxv]
    return u, v, p
