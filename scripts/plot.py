"""Plot fields and force history from a checkpoint (Python port of postpro.m)."""

import argparse
from pathlib import Path

from afc import Simulation
from afc.postprocess import load_matlab_history, plot_fields, plot_history, plot_streamlines


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, help="checkpoint .npz, or a directory of MATLAB .mat files")
    p.add_argument("--out", type=Path, help="output directory (default: next to the source)")
    p.add_argument("--reference", type=Path, help="MATLAB directory whose Cdd/Cll/ango history is overlaid")
    args = p.parse_args()

    if args.source.is_dir():
        sim = Simulation.from_matlab(args.source)
        history = load_matlab_history(args.source)
        out = args.out or args.source / "plots"
    else:
        sim = Simulation.load(args.source)
        history = sim.history_arrays()
        out = args.out or args.source.parent
    out.mkdir(parents=True, exist_ok=True)

    plot_fields(sim, out / "fields.png")
    plot_streamlines(sim, out / "streamlines.png")
    if len(history.get("t", [])):
        reference = load_matlab_history(args.reference) if args.reference else None
        plot_history(history, out / "history.png", reference=reference)
    print(f"plots written to {out}")


if __name__ == "__main__":
    main()
