"""Run the flapping-tail cylinder simulation (Python port of NavierStokesENV.m)."""

import argparse
import math
import time
from pathlib import Path

import numpy as np

from afc import Simulation, SimulationConfig
from afc.live import LiveView
from afc.postprocess import load_matlab_history, plot_fields, plot_history, plot_streamlines


def parse_args():
    defaults = SimulationConfig()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--steps", type=int, default=1000)
    p.add_argument("--re", type=float, default=defaults.re)
    p.add_argument("--ncx", type=int, default=defaults.ncx)
    p.add_argument("--ncy", type=int, default=defaults.ncy)
    p.add_argument("--cfl", type=float, default=defaults.cfl)
    p.add_argument("--angle-deg", type=float, default=math.degrees(defaults.initial_angle), help="initial tail angle")
    p.add_argument("--no-flap", action="store_true", help="hold the tail fixed instead of solving its dynamics")
    p.add_argument("--no-tail", action="store_true", help="bare cylinder")
    p.add_argument("--rkpm", action="store_true", help="RKPM-corrected delta kernel")
    p.add_argument("--pivot", choices=["center", "hinge"], default=defaults.pivot)
    p.add_argument("--k-spring", type=float, default=defaults.k_spring)
    p.add_argument("--restart", type=Path, help="checkpoint .npz, or a directory of MATLAB .mat files")
    p.add_argument("--out", type=Path, default=Path("runs") / time.strftime("%Y%m%d-%H%M%S"))
    p.add_argument("--report-every", type=int, default=20)
    p.add_argument("--checkpoint-every", type=int, default=500)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--live", action="store_true", help="watch the flow in a window while it runs (close it to stop)")
    p.add_argument("--live-every", type=int, default=20, help="steps between live redraws")
    p.add_argument("--vmax", type=float, default=6.0, help="vorticity colour limit for the live view")
    p.add_argument("--frames", type=Path, help="save a PNG of the live view to this folder at every redraw")
    return p.parse_args()


def main():
    args = parse_args()
    config = SimulationConfig(
        ncx=args.ncx,
        ncy=args.ncy,
        re=args.re,
        cfl=args.cfl,
        initial_angle=math.radians(args.angle_deg),
        tail=not args.no_tail,
        flap=not (args.no_flap or args.no_tail),
        rkpm=args.rkpm,
        pivot=args.pivot,
        k_spring=args.k_spring,
    )
    if args.restart is None:
        sim = Simulation(config)
    elif args.restart.is_dir():
        sim = Simulation.from_matlab(args.restart, config)
    else:
        sim = Simulation.load(args.restart)

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"grid {sim.grid.ncx}x{sim.grid.ncy}, Re={sim.config.re:g}, dt={sim.dt:.4g}, markers={sim.markers.n}")
    print(f"{'step':>7} {'t':>8} {'Cd':>9} {'Cl':>9} {'angle':>8} {'max div':>9} {'slip':>9}")

    view = None
    if args.live or args.frames:
        view = LiveView(sim, vmax=args.vmax, frames_dir=args.frames, show=args.live)
        view.update(sim)

    started = time.perf_counter()
    first = sim.step_count
    target = first + args.steps
    try:
        while sim.step_count < target:
            sim.step()
            n = sim.step_count
            if n % args.report_every == 0:
                h = sim.history
                print(
                    f"{n:7d} {h['t'][-1]:8.3f} {h['cd'][-1]:9.4f} {h['cl'][-1]:9.4f} "
                    f"{h['angle'][-1]:8.4f} {sim.max_divergence():9.2e} {h['slip_error'][-1]:9.2e}",
                    flush=True,
                )
            if n % args.checkpoint_every == 0:
                sim.save(args.out / "checkpoint.npz")
            if view and n % args.live_every == 0 and not view.update(sim):
                print("live window closed, stopping")
                break
    except KeyboardInterrupt:
        print("interrupted, saving results")
    done = sim.step_count - first
    elapsed = time.perf_counter() - started
    print(f"{done} steps in {elapsed:.1f}s ({1000 * elapsed / max(done, 1):.1f} ms/step)")

    sim.save(args.out / "checkpoint.npz")
    np.savez(args.out / "history.npz", **sim.history_arrays())
    if not args.no_plots:
        reference = load_matlab_history(args.restart) if args.restart and args.restart.is_dir() else None
        plot_fields(sim, args.out / "fields.png")
        plot_streamlines(sim, args.out / "streamlines.png")
        plot_history(sim.history_arrays(), args.out / "history.png", reference=reference)
    print(f"results written to {args.out}")
    if view:
        view.close()


if __name__ == "__main__":
    main()
