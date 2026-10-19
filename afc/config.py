from dataclasses import asdict, dataclass
import math


@dataclass
class SimulationConfig:
    """Physical and numerical parameters (defaults match NavierStokesENV.m)."""

    ncx: int = 500
    ncy: int = 420
    lx: float = 16.5
    ly: float = 13.86

    re: float = 150.0
    diameter: float = 1.0
    xc: float = 3.003
    yc: float = 6.930

    cfl: float = 0.2
    inflow: float = 1.0

    tail: bool = True
    flap: bool = True
    rkpm: bool = False
    n_circle: int | None = None

    initial_angle: float = math.pi / 6
    k_spring: float = 0.5
    rho_beam: float = 2.0
    # "center" reproduces the thesis code (tail velocity and torque taken about the
    # cylinder centre); "hinge" uses the rear stagnation point the tail rotates about.
    pivot: str = "center"

    ab_restart_interval: int = 200

    @property
    def radius(self) -> float:
        return self.diameter / 2

    @property
    def hx(self) -> float:
        return self.lx / self.ncx

    @property
    def hy(self) -> float:
        return self.ly / self.ncy

    @property
    def dt(self) -> float:
        return self.cfl * self.hx / self.inflow

    @property
    def inertia(self) -> float:
        length = self.diameter
        centre = length
        delta_rho = self.rho_beam - 1
        return delta_rho * length**3 / 3 + delta_rho * length * centre * (centre - length)

    def to_dict(self) -> dict:
        return asdict(self)
