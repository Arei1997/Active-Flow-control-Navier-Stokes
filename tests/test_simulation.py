from pathlib import Path

import numpy as np
import pytest

from afc import Simulation, SimulationConfig
from afc.postprocess import load_matlab_history

MATLAB_DIR = Path(__file__).resolve().parents[1] / "matlab"
COARSE = dict(ncx=100, ncy=84)


def test_projection_is_divergence_free():
    sim = Simulation(SimulationConfig(**COARSE))
    sim.run(5)
    assert sim.max_divergence() < 1e-10
    assert np.all(np.isfinite(sim.u)) and np.all(np.isfinite(sim.v))


def test_checkpoint_restart_is_exact(tmp_path):
    reference = Simulation(SimulationConfig(**COARSE))
    reference.run(12)

    sim = Simulation(SimulationConfig(**COARSE))
    sim.run(7)
    sim.save(tmp_path / "state.npz")
    restarted = Simulation.load(tmp_path / "state.npz")
    restarted.run(5)

    np.testing.assert_array_equal(restarted.u, reference.u)
    np.testing.assert_array_equal(restarted.history["cd"], reference.history["cd"])
    assert restarted.angle == reference.angle


def test_fixed_tail_does_not_move():
    sim = Simulation(SimulationConfig(**COARSE, flap=False, initial_angle=0.2))
    sim.run(3)
    assert sim.angle == 0.2
    assert set(sim.history["angle_dot"]) == {0.0}


@pytest.mark.skipif(not (MATLAB_DIR / "Cdd.mat").exists(), reason="MATLAB reference output not present")
def test_matches_matlab_reference_run():
    """The saved MATLAB run is NavierStokesENV.m with Re=150 and an initial angle of pi/3."""
    ref = load_matlab_history(MATLAB_DIR)
    sim = Simulation(SimulationConfig(initial_angle=np.pi / 3))
    sim.run(40)
    h = sim.history_arrays()
    snapshots = np.array([19, 39])
    for key in ("t", "cd", "cl", "angle"):
        np.testing.assert_allclose(h[key][snapshots], ref[key][:2], rtol=1e-10)
    # MATLAB stores the torque before updating it, so its snapshot lags by one step
    np.testing.assert_allclose(h["torque"][snapshots - 1], ref["torque"][:2], rtol=1e-10)
