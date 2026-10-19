"""Reinforcement-learning interface: the agent sets the tail angle to reduce drag.

Forces follow the thesis sign convention: `cd` is the force the body exerts on
the fluid, so it is negative for a drag-producing body and drag = -cd.
"""

from dataclasses import replace
from pathlib import Path

import gymnasium as gym
import numpy as np
from gymnasium import spaces

from .config import SimulationConfig
from .simulation import Simulation


def drag_for_angle(angle: float, n_steps: int = 50, config: SimulationConfig | None = None) -> float:
    """Port of solver.m: impulsive start with the tail held at `angle`, returns Cd."""
    config = replace(config or SimulationConfig(), initial_angle=angle, flap=False, tail=True)
    sim = Simulation(config)
    sim.run(n_steps)
    return sim.history["cd"][-1]


def default_probes(config: SimulationConfig) -> np.ndarray:
    d = config.diameter
    xs = config.xc + d * np.array([0.75, 1.5, 2.5, 3.5])
    ys = config.yc + d * np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    return np.array([(x, y) for x in xs for y in ys])


class CylinderFlapEnv(gym.Env):
    """Each action is a target tail angle, reached by a linear ramp over `steps_per_action` CFD steps.

    The ramp rate is capped at `max_rate` (rad per convective time D/U) so the tail
    tip cannot outrun the flow and break the CFL limit.

    Observation: [Cd, Cl, angle / max_angle, u and v at the wake probes].
    Reward: mean drag reduction relative to `baseline_drag`, minus a lift penalty.
    """

    metadata = {"render_modes": []}

    def __init__(
        self,
        config: SimulationConfig | None = None,
        initial_state: str | Path | None = None,
        steps_per_action: int = 50,
        max_actions: int = 100,
        max_angle: float = np.pi / 3,
        max_rate: float = 1.0,
        lift_penalty: float = 0.2,
        baseline_drag: float = 0.0,
        probes: np.ndarray | None = None,
        warmup_steps: int = 0,
    ):
        super().__init__()
        self.config = replace(config or SimulationConfig(), flap=False, tail=True)
        self.initial_state = Path(initial_state) if initial_state else None
        self.steps_per_action = steps_per_action
        self.max_actions = max_actions
        self.max_angle = max_angle
        self.max_rate = max_rate
        self.lift_penalty = lift_penalty
        self.baseline_drag = baseline_drag
        self.probes = default_probes(self.config) if probes is None else np.asarray(probes)
        self.warmup_steps = warmup_steps

        n_obs = 3 + 2 * len(self.probes)
        self.observation_space = spaces.Box(-np.inf, np.inf, shape=(n_obs,), dtype=np.float32)
        self.action_space = spaces.Box(-1.0, 1.0, shape=(1,), dtype=np.float32)
        self.sim: Simulation | None = None
        self._n_actions = 0

    def _fresh_simulation(self) -> Simulation:
        if self.initial_state is None:
            sim = Simulation(self.config)
            sim.run(self.warmup_steps)
        elif self.initial_state.is_dir():
            sim = Simulation.from_matlab(self.initial_state, self.config)
        else:
            sim = Simulation.load(self.initial_state, flap=False, tail=True)
        sim.angle_dot = 0.0
        return sim

    def _observation(self) -> np.ndarray:
        h = self.sim.history
        cd = h["cd"][-1] if h["cd"] else 0.0
        cl = h["cl"][-1] if h["cl"] else 0.0
        vel = self.sim.probe(self.probes).ravel()
        return np.concatenate([[cd, cl, self.sim.angle / self.max_angle], vel]).astype(np.float32)

    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        self.sim = self._fresh_simulation()
        self._n_actions = 0
        return self._observation(), {}

    def step(self, action):
        sim = self.sim
        duration = self.steps_per_action * sim.dt
        target = float(np.clip(action, -1, 1)[0]) * self.max_angle
        target = float(np.clip(target, sim.angle - self.max_rate * duration, sim.angle + self.max_rate * duration))
        sim.angle_dot = (target - sim.angle) / duration
        start = len(sim.history["cd"])
        sim.run(self.steps_per_action)
        sim.angle_dot = 0.0
        sim.set_tail_angle(target)

        drag = -np.mean(sim.history["cd"][start:])
        lift = np.mean(sim.history["cl"][start:])
        reward = self.baseline_drag - drag - self.lift_penalty * abs(lift)
        self._n_actions += 1

        terminated = not np.isfinite(drag) or np.abs(sim.u).max() > 10 * sim.config.inflow
        truncated = self._n_actions >= self.max_actions
        info = {"drag": drag, "lift": lift, "angle": target, "t": sim.t}
        return self._observation(), float(reward), terminated, truncated, info
