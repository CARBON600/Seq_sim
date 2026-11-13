"""Pytest configuration: make the package importable whatever the checkout is called.

The package root is this repository's root (``pyproject.toml`` maps "." to
``seq_sim_v2``), so ``import seq_sim_v2`` normally needs a directory of that name
on ``sys.path``. A GitHub zip unpacks to ``Seq_sim-main/`` and a clone may be
called anything, so when the package is not already importable the repository is
registered under its own name from here.
"""

import importlib.util
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(_ROOT))

if importlib.util.find_spec("seq_sim_v2") is None:
    _spec = importlib.util.spec_from_file_location(
        "seq_sim_v2",
        os.path.join(_ROOT, "__init__.py"),
        submodule_search_locations=[_ROOT],
    )
    _module = importlib.util.module_from_spec(_spec)
    sys.modules["seq_sim_v2"] = _module
    _spec.loader.exec_module(_module)
