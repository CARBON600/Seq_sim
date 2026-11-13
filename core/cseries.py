"""C-series empirical profile mixtures (Si Quang, Gascuel & Lartillot 2008).

Loads the bundled C10-C60 amino-acid equilibrium frequency profiles and weights
from ``data/cseries/`` (transcribed from IQ-TREE 2's ``model/modelmixture.cpp``
in ``ARNDCQEGHILKMFPSTWYV`` order) and validates row and weight normalisation.
"""

from importlib.resources import files

import numpy as np

AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
CSERIES = ("C10", "C20", "C30", "C40", "C50", "C60")

_DATA_DIR = files("seq_sim_v2") / "data" / "cseries"


def load_cseries(name):
    """Load one C-series model.

    Returns ``(freqs, weights)``: ``freqs`` is a ``(C, 20)`` float64 array of
    profile equilibrium frequencies (each row sums to 1) and ``weights`` is a
    length-``C`` array of profile weights (sums to 1).
    """
    name = name.upper()
    if name not in CSERIES:
        raise ValueError(
            f"unknown C-series model {name!r}; expected one of {CSERIES}")
    freqs = _read_matrix(_DATA_DIR / f"{name}.freq.dat")
    weights = _read_vector(_DATA_DIR / f"{name}.weight.dat")
    n = freqs.shape[0]
    if n != weights.shape[0]:
        raise ValueError(f"{name}: {n} profiles but {weights.shape[0]} weights")
    if freqs.shape[1] != len(AA_ORDER):
        raise ValueError(f"{name}: {freqs.shape[1]} columns, expected 20")
    for i in range(n):
        s = float(freqs[i].sum())
        if abs(s - 1.0) > 1e-6:
            raise ValueError(f"{name} profile {i}: frequencies sum to {s}")
    w = float(weights.sum())
    if abs(w - 1.0) > 1e-6:
        raise ValueError(f"{name}: weights sum to {w}")
    return freqs, weights


def _read_matrix(path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    if not rows:
        raise ValueError(f"empty data file: {path}")
    return np.array(rows, dtype=np.float64)


def _read_vector(path):
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values.extend(float(x) for x in line.split())
    if not values:
        raise ValueError(f"empty data file: {path}")
    return np.array(values, dtype=np.float64)
