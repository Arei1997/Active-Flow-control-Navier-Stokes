"""Train a PPO agent to choose the tail angle that minimises drag.

Requires the optional RL dependencies: pip install -e ".[rl]"
"""

import argparse
import importlib.util
from pathlib import Path

from afc import SimulationConfig
from afc.env import CylinderFlapEnv


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--episodes", type=int, default=100, help="MaxEpisodes in the MATLAB agent setup")
    p.add_argument("--episode-length", type=int, default=50, help="actions per episode (MaxStepsPerEpisode)")
    p.add_argument("--steps-per-action", type=int, default=50, help="CFD time steps between actions")
    p.add_argument("--initial-state", type=Path, default=Path("matlab"), help="developed flow to start each episode from")
    p.add_argument("--ncx", type=int, default=500)
    p.add_argument("--ncy", type=int, default=420)
    p.add_argument("--out", type=Path, default=Path("runs/ppo"))
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    from stable_baselines3 import PPO
    from stable_baselines3.common.callbacks import CheckpointCallback
    from stable_baselines3.common.monitor import Monitor

    args.out.mkdir(parents=True, exist_ok=True)
    initial = args.initial_state if args.initial_state and args.initial_state.exists() else None
    env = Monitor(
        CylinderFlapEnv(
            SimulationConfig(ncx=args.ncx, ncy=args.ncy),
            initial_state=initial,
            steps_per_action=args.steps_per_action,
            max_actions=args.episode_length,
        ),
        str(args.out / "monitor"),
    )
    model = PPO(
        "MlpPolicy",
        env,
        n_steps=min(512, args.episode_length * 4),
        batch_size=64,
        learning_rate=3e-4,
        gamma=0.99,
        seed=args.seed,
        verbose=1,
        tensorboard_log=str(args.out / "tb") if importlib.util.find_spec("tensorboard") else None,
    )
    model.learn(
        total_timesteps=args.episodes * args.episode_length,
        callback=CheckpointCallback(save_freq=1000, save_path=str(args.out / "checkpoints"), name_prefix="ppo"),
    )
    model.save(args.out / "ppo_final")
    print(f"model saved to {args.out / 'ppo_final.zip'}")


if __name__ == "__main__":
    main()
