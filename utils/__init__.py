"""The run's random generator.

Model       : one module-level generator shared by every draw in a run.
Equation    : not applicable.
Source      : not applicable; this project's own.
Convention  : seeding is global rather than per-call, which is what makes a given
              ``-z`` reproduce a run byte for byte. ``random`` returns a uniform
              variate in [0, 1) and ``gamma`` the gamma sampler used by the
              continuous-gamma path.
Invariants  : the same seed produces the same sequence of draws, and therefore
              the same output bytes. Because the generator is global, the package
              is not thread-safe; a concurrent run must use a separate process.
"""

import time

import numpy as np

_rng = np.random.default_rng()

__all__ = ["set_seed", "create_seed", "random", "gamma"]


def set_seed(seed):
    global _rng
    _rng = np.random.default_rng(seed)


def create_seed():
    """Pick a seed when the user does not supply one."""
    return int(time.time() * 1e6) & 0xFFFFFFFF


def random():
    return _rng.random()


def gamma(shape):
    return _rng.gamma(shape, 1.0)
