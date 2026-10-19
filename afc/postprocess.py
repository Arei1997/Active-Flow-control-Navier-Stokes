import numpy as np
from scipy.fft import dstn, idstn

from .simulation import Simulation


def vorticity(u: np.ndarray, v: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """Vorticity at cell corners, shape (ncx+1, ncy+1)."""
    return np.diff(v, axis=0) / hx - np.diff(u, axis=1) / hy


def streamfunction(omega: np.ndarray, hx: float, hy: float) -> np.ndarray:
    """Solve lap(psi) = -omega on interior corners with psi = 0 on the frame (vort.m)."""
    f = -omega[1:-1, 1:-1]
    nx, ny = f.shape
    lam_x = -4 / hx**2 * np.sin(np.pi * np.arange(1, nx + 1) / (2 * (nx + 1))) ** 2
    lam_y = -4 / hy**2 * np.sin(np.pi * np.arange(1, ny + 1) / (2 * (ny + 1))) ** 2
    psi_hat = dstn(f, type=1, norm="ortho") / (lam_x[:, None] + lam_y[None, :])
    psi = np.zeros_like(omega)
    psi[1:-1, 1:-1] = idstn(psi_hat, type=1, norm="ortho")
    return psi


def cell_centred_velocity(u: np.ndarray, v: np.ndarray):
    return 0.5 * (u[1:, 1:-1] + u[:-1, 1:-1]), 0.5 * (v[1:-1, 1:] + v[1:-1, :-1])


def _wake_window(sim: Simulation, upstream=0.75, downstream=7.75, half_height=1.5):
    c = sim.config
    return (c.xc - upstream * c.diameter, c.xc + downstream * c.diameter), (
        c.yc - half_height * c.diameter,
        c.yc + half_height * c.diameter,
    )


def _draw_body(ax, sim: Simulation):
    m = sim.markers
    ax.plot(m.x, m.y, "k.", ms=2)


def plot_fields(sim: Simulation, path=None, window=None):
    """Streamwise velocity, cross-stream velocity, pressure and vorticity around the body."""
    import matplotlib.pyplot as plt

    g = sim.grid
    (x0, x1), (y0, y1) = window or _wake_window(sim)
    omega = vorticity(sim.u, sim.v, g.hx, g.hy)
    xo = np.arange(omega.shape[0]) * g.hx
    yo = np.arange(omega.shape[1]) * g.hy
    lim = 0.8 * np.abs(omega[(xo > x0) & (xo < x1)][:, (yo > y0) & (yo < y1)]).max()

    panels = [
        ("u", g.xu, g.yu, sim.u, "viridis", None),
        ("v", g.xv, g.yv, sim.v, "RdBu_r", None),
        ("p", g.xp, g.yp, sim.p, "RdBu_r", None),
        (r"$\omega$", xo, yo, omega, "RdBu_r", np.linspace(-lim, lim, 101)),
    ]
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 2.6 * len(panels)), constrained_layout=True)
    for ax, (name, x, y, f, cmap, levels) in zip(axes, panels):
        ix = (x >= x0) & (x <= x1)
        jy = (y >= y0) & (y <= y1)
        cs = ax.contourf(x[ix], y[jy], f[np.ix_(ix, jy)].T, levels=levels if levels is not None else 100, cmap=cmap, extend="both")
        _draw_body(ax, sim)
        ax.set_aspect("equal")
        ax.set_title(f"{name}   (Re = {sim.config.re:g}, t = {sim.t:.2f})")
        fig.colorbar(cs, ax=ax, shrink=0.9)
    if path:
        fig.savefig(path, dpi=150)
        plt.close(fig)
    return fig


def plot_streamlines(sim: Simulation, path=None, window=None, density=2.0):
    import matplotlib.pyplot as plt

    g = sim.grid
    (x0, x1), (y0, y1) = window or _wake_window(sim)
    uc, vc = cell_centred_velocity(sim.u, sim.v)
    ix = (g.xp >= x0) & (g.xp <= x1)
    jy = (g.yp >= y0) & (g.yp <= y1)
    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    speed = np.hypot(uc, vc)[np.ix_(ix, jy)].T
    ax.streamplot(g.xp[ix], g.yp[jy], uc[np.ix_(ix, jy)].T, vc[np.ix_(ix, jy)].T, color=speed, density=density, linewidth=0.7)
    ax.fill(sim.markers.x[: sim.markers.n_circle], sim.markers.y[: sim.markers.n_circle], color="0.3")
    ax.plot(sim.markers.x[sim.markers.tail], sim.markers.y[sim.markers.tail], "k-", lw=2)
    ax.set_aspect("equal")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    if path:
        fig.savefig(path, dpi=150)
        plt.close(fig)
    return fig


def plot_history(history: dict, path=None, reference: dict | None = None):
    """Force coefficients and tail motion against time; `reference` overlays another run."""
    import matplotlib.pyplot as plt

    keys = [k for k in ("cd", "cl", "angle", "torque") if k in history and len(history[k])]
    labels = {"cd": "$C_D$", "cl": "$C_L$", "angle": r"$\theta$ [rad]", "torque": "torque"}
    fig, axes = plt.subplots(len(keys), 1, figsize=(9, 2.2 * len(keys)), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, key in zip(axes, keys):
        if reference is not None and key in reference:
            ax.plot(reference["t"], reference[key], color="0.6", lw=1, label="reference")
        ax.plot(history["t"], history[key], lw=1, label="python")
        ax.set_ylabel(labels[key])
        ax.grid(alpha=0.3)
    if reference is not None:
        axes[0].legend(loc="best", frameon=False)
    axes[-1].set_xlabel("t")
    if path:
        fig.savefig(path, dpi=150)
        plt.close(fig)
    return fig


def load_matlab_history(directory) -> dict:
    """Snapshot series saved by the MATLAB driver (time, Cdd, Cll, ango, torques)."""
    from pathlib import Path

    import scipy.io as sio

    directory = Path(directory)
    names = {"t": ("time", "time"), "cd": ("Cdd", "Cdd"), "cl": ("Cll", "Cll"), "angle": ("ango", "ango"), "torque": ("torques", "torques")}
    out = {}
    for key, (file, var) in names.items():
        if (directory / f"{file}.mat").exists():
            out[key] = sio.loadmat(directory / f"{file}.mat")[var].ravel()
    return out
