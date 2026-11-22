"""seq_sim: Monte Carlo simulation of molecular sequence evolution."""

from ._version import __version__
from .seq_sim import main

VERSION = __version__

__all__ = ["main", "VERSION", "__version__"]
