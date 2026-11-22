"""Pytest configuration: make src/ importable whatever the checkout is called.

The CLI tests run python -m seq_sim in a subprocess, so src/ is exported
through PYTHONPATH and inherited by those children.
"""

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SRC = os.path.join(_ROOT, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
_TESTS = os.path.join(_ROOT, "tests")
if _TESTS not in sys.path:
    sys.path.insert(0, _TESTS)

_OLD = os.environ.get("PYTHONPATH", "")
os.environ["PYTHONPATH"] = _SRC + (os.pathsep + _OLD if _OLD else "")
