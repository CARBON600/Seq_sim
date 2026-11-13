"""GFmix site and branch heterogeneity (McCarthy, Susko, Harada & Roger 2026).

Implements per-site C-series profile sampling and per-branch GARP/FYMINK
compositional rescaling (Eqs. 1.4 and 1.7 of the McCarthy et al. bioRxiv
preprint). On heterogeneous branches, GARP (G, A, R, P) frequencies are scaled
by ``gamma_G``, FYMINK (F, Y, M, I, N, K) by ``gamma_F``, and the remaining ten
amino acids by 1, then renormalised to sum to 1 while keeping the model's
exchangeabilities unchanged.
"""

import numpy as np

from .aa_models import AA_ORDER, AAModel
from .cseries import load_cseries

GARP = "GARP"
FYMINK = "FYMINK"

_GARP_IDX = np.array([AA_ORDER.index(a) for a in GARP], dtype=np.int64)
_FYMINK_IDX = np.array([AA_ORDER.index(a) for a in FYMINK], dtype=np.int64)


def branch_transform(freqs, gamma_g, gamma_f):
    """Apply Eq. 1.4 to one profile frequency vector.

    ``freqs`` is a length-20 profile summing to 1. Returns a length-20 vector
    that again sums to 1: the GARP class is scaled by ``gamma_g``, the FYMINK
    class by ``gamma_f`` and the remaining class by 1, followed by
    renormalisation via the factor ``mu``.
    """
    freqs = np.asarray(freqs, dtype=np.float64)
    sum_g = float(freqs[_GARP_IDX].sum())
    sum_f = float(freqs[_FYMINK_IDX].sum())
    sum_o = 1.0 - sum_g - sum_f
    mu = 1.0 / (gamma_g * sum_g + gamma_f * sum_f + sum_o)
    out = freqs * mu
    out[_GARP_IDX] *= gamma_g
    out[_FYMINK_IDX] *= gamma_f
    return out


def assign_site_profiles(num_sites, weights, rng):
    """Draw a profile index for every site from ``weights``.

    ``rng`` is a callable returning a uniform variate in ``[0, 1)``. Draws are
    independent across sites and are redrawn once per replicate, matching the
    per-replicate semantics of the gamma category assignment.
    """
    cum = np.cumsum(weights)
    out = np.zeros(num_sites, dtype=np.int64)
    n = len(weights)
    for i in range(num_sites):
        r = rng()
        idx = int(np.searchsorted(cum, r, side="right"))
        out[i] = idx if idx < n else n - 1
    return out


def build_profile_models(cseries_name, exchangeability, gamma_by_node):
    """Build the AA models for every (branch, profile) pair.

    ``exchangeability`` is the 190-element rate vector shared by all profiles
    (the ``-m`` model's exchangeabilities, e.g. LG). ``gamma_by_node`` maps a
    node identifier to its ``(gamma_g, gamma_f)`` pair; nodes absent from the
    map are homogeneous (gamma 1, 1) and share the base models.

    Returns ``(base_models, per_node_models, weights)``: ``base_models`` is a
    length-``C`` list of AA models on the raw C-series profiles;
    ``per_node_models`` maps each heterogeneous node to a length-``C`` list of
    models on the branch-transformed profiles; ``weights`` are the C-series
    built-in profile weights.
    """
    freqs, weights = load_cseries(cseries_name)
    n = freqs.shape[0]
    # Rows carry ~1e-10 rounding from the 10-significant-digit transcription;
    # renormalise so each profile sums to exactly 1 before building Q.
    freqs = freqs / freqs.sum(axis=1, keepdims=True)
    base = [AAModel(freqs[c], exchangeability) for c in range(n)]
    per_node = {}
    for node_id, (gg, gf) in gamma_by_node.items():
        per_node[node_id] = [
            AAModel(branch_transform(freqs[c], gg, gf), exchangeability)
            for c in range(n)
        ]
    return base, per_node, weights
