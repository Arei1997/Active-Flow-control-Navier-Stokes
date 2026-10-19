"""Flow past a cylinder with a hinged, spring-loaded tail.

Fractional-step scheme from NavierStokesENV.m:
  1. explicit predictor (Adams-Bashforth convection, Crank-Nicolson diffusion)
  2. immersed-boundary forcing so the predictor satisfies no-slip on the body
  3. implicit viscous step by D'Yakonov approximate factorisation
  4. pressure projection with a rotational pressure update
  5. tail angle advanced from the fluid torque (torsional spring + inertia)
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy.interpolate import RegularGridInterpolator

from .config import SimulationConfig
from .grid import StaggeredGrid, boundary_velocities, initial_fields
from .immersed import ImmersedBoundary, build_markers
from .linalg import NeumannPoisson
from .operators import (
    advection_u,
    advection_v,
    apply_boundary_conditions,
    convective_outflow_u,
    convective_outflow_v,
    divergence,
    explicit_rhs,
    laplacian,
    pressure_gradient,
    project,
)
from .viscous import ViscousSolverU, ViscousSolverV

_HISTORY_KEYS = ("t", "cd", "cl", "force_x", "force_y", "angle", "angle_dot", "torque", "slip_error")


class Simulation:
    def __init__(self, config: SimulationConfig | None = None):
        self.config = c = config or SimulationConfig()
        self.grid = g = StaggeredGrid.from_config(c)
        self.dt = c.dt
        self.ub, self.vb = boundary_velocities(g, c.inflow)

        self.visc_u = ViscousSolverU(g.nxu, g.nyu, g.hx, g.hy, self.dt, c.re)
        self.visc_v = ViscousSolverV(g.nxv, g.nyv, g.hx, g.hy, self.dt, c.re)
        self.poisson = NeumannPoisson(c.ncx, c.ncy, g.hx, g.hy)

        self.u, self.v, self.p = initial_fields(g, self.ub, self.vb)
        self.phi = np.zeros_like(self.p)
        self.angle = c.initial_angle if c.tail else 0.0
        self.angle_dot = 0.0
        self.torque = 0.0

        self.markers = build_markers(g, c.xc, c.yc, c.radius, self.angle, c.tail, c.n_circle)
        self._build_ib()

        self.t = 0.0
        self.step_count = 0
        self.history = defaultdict(list)
        self._reset_multistep()

    def _reset_multistep(self):
        g = self.grid
        self.q_in = self.u[-1, 1:-1].sum()
        self._first = True
        self._lap_u = np.zeros((g.nxu - 2, g.nyu - 2))
        self._nl_u = np.zeros_like(self._lap_u)
        self._lap_v = np.zeros((g.nxv - 2, g.nyv - 2))
        self._nl_v = np.zeros_like(self._lap_v)
        self._torque_rhs_old = 0.0

    def _build_ib(self):
        g, rkpm = self.grid, self.config.rkpm
        self.ib_u = ImmersedBoundary(g.xu[1:-1], g.yu[1:-1], self.markers, rkpm)
        self.ib_v = ImmersedBoundary(g.xv[1:-1], g.yv[1:-1], self.markers, rkpm)

    @property
    def pivot(self) -> tuple[float, float]:
        if self.config.pivot == "hinge":
            return self.markers.hinge
        return self.config.xc, self.config.yc

    def set_tail_angle(self, angle: float):
        self.angle = angle
        self.markers = self.markers.with_tail_angle(angle)
        self._build_ib()

    def _wall_velocity(self):
        """Rigid-rotation velocity of the tail markers (zero on the cylinder)."""
        m = self.markers
        uw, vw = np.zeros(m.n), np.zeros(m.n)
        if self.config.tail:
            px, py = self.pivot
            uw[m.tail] = -self.angle_dot * (m.y[m.tail] - py)
            vw[m.tail] = self.angle_dot * (m.x[m.tail] - px)
        return uw, vw

    def _coefficients(self):
        if self._first:
            return 1.0, 0.0, 0.0
        if (self.step_count + 1) % self.config.ab_restart_interval == 0:
            return 1.0, 0.0, 1.0
        return 1.5, 0.5, 1.0

    def step(self):
        c, g, dt = self.config, self.grid, self.dt
        hx, hy, re = g.hx, g.hy, c.re
        ab1, ab2, cgp = self._coefficients()
        u, v, p = self.u, self.v, self.p
        m = self.markers
        uw, vw = self._wall_velocity()

        lap = laplacian(u, hx, hy)
        nl = advection_u(u, v, hx, hy)
        gp = cgp * pressure_gradient(p, hx, hy, axis=0)
        rhs = explicit_rhs(u, lap, nl, self._lap_u, self._nl_u, gp, dt, re, ab1, ab2)
        f_x, val_x = self.ib_u.forcing(rhs, uw / dt)
        force_x = hx * hy * f_x.sum()
        u_out = convective_outflow_u(u, dt, hx, self.q_in)
        u_star = self.visc_u.solve(f_x - ab1 * nl + ab2 * self._nl_u - gp, u, self.ub[3, : g.nyu], u_out)
        self._lap_u, self._nl_u = lap, nl

        lap = laplacian(v, hx, hy)
        nl = advection_v(u, v, hx, hy)
        gp = cgp * pressure_gradient(p, hx, hy, axis=1)
        rhs = explicit_rhs(v, lap, nl, self._lap_v, self._nl_v, gp, dt, re, ab1, ab2)
        f_y, val_y = self.ib_v.forcing(rhs, vw / dt)
        force_y = hx * hy * f_y.sum()
        v_out = convective_outflow_v(u, v, dt, hx)
        vb = self.vb
        v_star = self.visc_v.solve(
            f_y - ab1 * nl + ab2 * self._nl_v - gp,
            v,
            vb[0, : g.nxv],
            vb[2, : g.nxv],
            vb[3, : g.nyv],
            v_out,
        )
        self._lap_v, self._nl_v = lap, nl

        u, v = apply_boundary_conditions(u_star, v_star, self.ub, vb, u_out, v_out)
        self.phi = self.poisson.solve(divergence(u, v, hx, hy) / dt)
        self.u, self.v, self.p = project(self.phi, u, v, p, dt, re, hx, hy)

        slip = max(
            np.abs(self.ib_u.interpolate(self.u[1:-1, 1:-1]) - uw).max(),
            np.abs(self.ib_v.interpolate(self.v[1:-1, 1:-1]) - vw).max(),
        )
        if c.tail:
            t = m.tail
            fx_lag = -val_x[t] * self.ib_u.epsilon[t] * m.ds[t]
            fy_lag = -val_y[t] * self.ib_v.epsilon[t] * m.ds[t]
            px, py = self.pivot
            self.torque = np.sum((m.x[t] - px) * fy_lag - (m.y[t] - py) * fx_lag)

        self._record(
            t=self.t,
            cd=2 * force_x / c.diameter,
            cl=2 * force_y / c.diameter,
            force_x=force_x,
            force_y=force_y,
            angle=self.angle,
            angle_dot=self.angle_dot,
            torque=self.torque,
            slip_error=slip,
        )

        if c.tail and c.flap:
            self._advance_flap(ab1, ab2)
        elif c.tail and self.angle_dot != 0.0:
            self.set_tail_angle(self.angle + dt * self.angle_dot)

        self.step_count += 1
        self.t += dt
        self._first = False

    def _advance_flap(self, ab1, ab2):
        """theta'' = -(k theta - T) / I with AB2 for velocity and trapezoid for angle."""
        c = self.config
        rhs = (c.k_spring * self.angle - self.torque) / c.inertia
        if ab2 == 0:
            self._torque_rhs_old = 0.0
        old = self.angle_dot
        self.angle_dot = old - self.dt * (ab1 * rhs - ab2 * self._torque_rhs_old)
        self._torque_rhs_old = rhs
        self.set_tail_angle(self.angle + 0.5 * self.dt * (old + self.angle_dot))

    def _record(self, **values):
        for key, value in values.items():
            self.history[key].append(float(value))

    def run(self, n_steps: int, callback=None, every: int = 20):
        for _ in range(n_steps):
            self.step()
            if callback is not None and self.step_count % every == 0:
                callback(self)

    def history_arrays(self) -> dict[str, np.ndarray]:
        return {k: np.asarray(self.history[k]) for k in _HISTORY_KEYS}

    def probe(self, points: np.ndarray) -> np.ndarray:
        """Velocity (u, v) at arbitrary points by bilinear interpolation, shape (n, 2)."""
        g = self.grid
        iu = RegularGridInterpolator((g.xu, g.yu), self.u)
        iv = RegularGridInterpolator((g.xv, g.yv), self.v)
        return np.column_stack([iu(points), iv(points)])

    def max_divergence(self) -> float:
        return float(np.abs(divergence(self.u, self.v, self.grid.hx, self.grid.hy)).max())

    def save(self, path: str | Path):
        """Checkpoint with full multistep state, so a restart continues exactly."""
        np.savez_compressed(
            path,
            config=json.dumps(self.config.to_dict()),
            u=self.u,
            v=self.v,
            p=self.p,
            phi=self.phi,
            angle=self.angle,
            angle_dot=self.angle_dot,
            torque=self.torque,
            t=self.t,
            step_count=self.step_count,
            q_in=self.q_in,
            first=self._first,
            lap_u=self._lap_u,
            nl_u=self._nl_u,
            lap_v=self._lap_v,
            nl_v=self._nl_v,
            torque_rhs_old=self._torque_rhs_old,
            **{f"history_{k}": v for k, v in self.history_arrays().items()},
        )

    @classmethod
    def load(cls, path: str | Path, **overrides) -> "Simulation":
        data = np.load(path)
        config = SimulationConfig(**{**json.loads(str(data["config"])), **overrides})
        sim = cls(config)
        sim.u, sim.v, sim.p, sim.phi = data["u"], data["v"], data["p"], data["phi"]
        sim.angle_dot, sim.torque = float(data["angle_dot"]), float(data["torque"])
        if config.tail:
            sim.set_tail_angle(float(data["angle"]))
        sim.t, sim.step_count = float(data["t"]), int(data["step_count"])
        sim.q_in, sim._first = float(data["q_in"]), bool(data["first"])
        sim._lap_u, sim._nl_u = data["lap_u"], data["nl_u"]
        sim._lap_v, sim._nl_v = data["lap_v"], data["nl_v"]
        sim._torque_rhs_old = float(data["torque_rhs_old"])
        for key in _HISTORY_KEYS:
            sim.history[key] = list(data[f"history_{key}"])
        return sim

    @classmethod
    def from_matlab(cls, directory: str | Path, config: SimulationConfig | None = None) -> "Simulation":
        """Restart from the .mat files written by the MATLAB code (restart=1, restart_flap=1)."""
        directory = Path(directory)

        def read(name, key):
            return sio.loadmat(directory / f"{name}.mat")[key]

        sim = cls(config)
        sim.u = read("uu", "u").astype(float)
        sim.v = read("vv", "v").astype(float)
        sim.p = read("pp", "p").astype(float)
        if sim.u.shape != (sim.grid.nxu, sim.grid.nyu):
            raise ValueError(f"grid mismatch: u is {sim.u.shape}, config expects {(sim.grid.nxu, sim.grid.nyu)}")
        if (directory / "phi.mat").exists():
            sim.phi = read("phi", "phi").astype(float)
        if sim.config.tail and (directory / "angolo.mat").exists():
            sim.angle_dot = read("angolo_dot", "angolo_dot").item()
            sim.torque = read("TORQUE", "TORQUE").item()
            sim.set_tail_angle(read("angolo", "angolo").item())
        if (directory / "time.mat").exists():
            sim.t = float(read("time", "time").ravel()[-1]) + sim.dt
        sim._reset_multistep()
        return sim
