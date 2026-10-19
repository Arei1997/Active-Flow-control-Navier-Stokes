import numpy as np
import pytest

gym = pytest.importorskip("gymnasium")

from afc import SimulationConfig
from afc.env import CylinderFlapEnv, drag_for_angle


def test_env_step_and_ramp():
    env = CylinderFlapEnv(SimulationConfig(ncx=100, ncy=84), steps_per_action=4, max_actions=2)
    obs, _ = env.reset(seed=0)
    assert obs.shape == env.observation_space.shape

    start = env.sim.angle
    obs, reward, terminated, truncated, info = env.step(np.array([0.9], dtype=np.float32))
    assert np.isfinite(reward) and not terminated and not truncated
    max_step = env.max_rate * env.steps_per_action * env.sim.dt
    assert env.sim.angle == pytest.approx(min(0.9 * env.max_angle, start + max_step))
    assert np.all(np.diff(env.sim.history["angle"]) > 0)

    *_, truncated, _ = env.step(np.array([-1.0], dtype=np.float32))
    assert truncated


def test_drag_for_angle_runs():
    cd = drag_for_angle(0.3, n_steps=3, config=SimulationConfig(ncx=100, ncy=84))
    assert np.isfinite(cd) and cd < 0
