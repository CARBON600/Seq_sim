"""Module-level PRNG wrapper around numpy.random.default_rng.
Seeding is global: a fixed -z seed reproduces a run byte-for-byte. Because the
generator is module-global, concurrent runs belong in separate processes.
input: a seed (set_seed)
output: random(), gamma(shape)
"""

import time

import numpy as np

_rng = np.random.default_rng()

__all__ = ["set_seed", "create_seed", "random", "gamma"]


###function definition

def set_seed(seed):
    global _rng
    _rng = np.random.default_rng(seed)


def create_seed():
    """Pick a time-based 32-bit seed when the user does not supply one."""
    return int(time.time() * 1e6) & 0xFFFFFFFF


def random():
    return _rng.random()


def gamma(shape):
    return _rng.gamma(shape, 1.0)
