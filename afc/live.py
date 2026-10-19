"""Live view of a running simulation (vorticity, body and force history), as in movieout.m."""

from pathlib import Path

import numpy as np

from .postprocess import _wake_window, vorticity
from .simulation import Simulation


class LiveView:
    def __init__(self, sim: Simulation, window=None, vmax: float = 6.0, frames_dir=None, show: bool = True):
        import matplotlib

        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        self.plt = plt
        self.show = show
        self.frames_dir = Path(frames_dir) if frames_dir else None
        if self.frames_dir:
            self.frames_dir.mkdir(parents=True, exist_ok=True)
        self.n_frames = 0

        g = sim.grid
        (x0, x1), (y0, y1) = window or _wake_window(sim, upstream=1.5, downstream=10.0, half_height=2.5)
        xo = np.arange(g.ncx + 1) * g.hx
        yo = np.arange(g.ncy + 1) * g.hy
        self.ix = (xo >= x0) & (xo <= x1)
        self.jy = (yo >= y0) & (yo <= y1)
        extent = (xo[self.ix][0], xo[self.ix][-1], yo[self.jy][0], yo[self.jy][-1])

        if show:
            plt.ion()
        self.fig, (self.ax_field, self.ax_force) = plt.subplots(
            2, 1, figsize=(11, 7.5), gridspec_kw={"height_ratios": [2.2, 1]}, constrained_layout=True
        )
        self.image = self.ax_field.imshow(
            np.zeros((self.jy.sum(), self.ix.sum())), origin="lower", extent=extent, cmap="RdBu_r", vmin=-vmax, vmax=vmax
        )
        self.fig.colorbar(self.image, ax=self.ax_field, label=r"vorticity $\omega$", shrink=0.85)
        n = sim.markers.n_circle
        (self.body,) = self.ax_field.fill(sim.markers.x[:n], sim.markers.y[:n], color="0.25")
        (self.tail,) = self.ax_field.plot([], [], "k-", lw=2.5)
        self.ax_field.set_aspect("equal")

        (self.cd_line,) = self.ax_force.plot([], [], lw=1, label="$C_D$")
        self.ax_cl = self.ax_force.twinx()
        (self.cl_line,) = self.ax_cl.plot([], [], lw=1, color="tab:orange", label="$C_L$")
        self.ax_force.set_xlabel("t")
        self.ax_force.set_ylabel("$C_D$", color="tab:blue")
        self.ax_cl.set_ylabel("$C_L$", color="tab:orange")
        self.ax_force.grid(alpha=0.3)

    @property
    def is_open(self) -> bool:
        return not self.show or self.plt.fignum_exists(self.fig.number)

    def update(self, sim: Simulation) -> bool:
        """Redraw; returns False once the window has been closed."""
        if not self.is_open:
            return False
        g, m = sim.grid, sim.markers
        omega = vorticity(sim.u, sim.v, g.hx, g.hy)
        self.image.set_data(omega[np.ix_(self.ix, self.jy)].T)
        self.body.set_xy(np.column_stack([m.x[: m.n_circle], m.y[: m.n_circle]]))
        self.tail.set_data(m.x[m.tail], m.y[m.tail])
        self.ax_field.set_title(
            f"Re = {sim.config.re:g}   t = {sim.t:.2f}   step {sim.step_count}   "
            f"tail angle = {np.degrees(sim.angle):+.1f}°"
        )

        h = sim.history
        if h["t"]:
            t = np.asarray(h["t"])
            cd, cl = np.asarray(h["cd"]), np.asarray(h["cl"])
            self.cd_line.set_data(t, cd)
            self.cl_line.set_data(t, cl)
            # skip the impulsive-start spike when scaling the axes
            recent = slice(100 if len(t) > 200 else 0, None)
            self.ax_force.set_xlim(t[0], max(t[-1], t[0] + 1e-9))
            self.ax_force.set_ylim(*self._limits(cd[recent]))
            self.ax_cl.set_ylim(*self._limits(cl[recent]))

        if self.frames_dir:
            self.n_frames += 1
            self.fig.savefig(self.frames_dir / f"{self.n_frames:05d}.png", dpi=100)
        if self.show:
            self.fig.canvas.draw_idle()
            self.plt.pause(0.001)
        return True

    @staticmethod
    def _limits(values):
        lo, hi = float(values.min()), float(values.max())
        pad = 0.1 * (hi - lo) or 0.05 * abs(hi) or 1e-3
        return lo - pad, hi + pad

    def close(self, keep_open: bool = True):
        if self.show and keep_open and self.is_open:
            self.plt.ioff()
            self.plt.show()
        else:
            self.plt.close(self.fig)
