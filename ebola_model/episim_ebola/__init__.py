"""EpiSim Ebola model package. See MODEL.md."""
from .model import (EbolaParams, Response, simulate_deterministic,
                    simulate_stochastic, summarise, reproduction_numbers)

__version__ = "0.1.0"
