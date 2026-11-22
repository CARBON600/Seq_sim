"""GFmix site and branch heterogeneity (McCarthy, Susko, Harada & Roger 2026).
Per-site C-series profile sampling plus per-branch GARP/FYMINK compositional
rescaling (Eqs. 1.4 and 1.7 of the McCarthy et al. bioRxiv preprint). On a
heterogeneous branch, GARP (G, A, R, P) frequencies are scaled by gamma_G,
FYMINK (F, Y, M, I, N, K) by gamma_F and the other ten amino acids by 1, then
renormalised; the exchangeabilities are unchanged.
input: C-series name or an explicit weight vector, per-node (gamma_G, gamma_F)
output: base models, per-node models, weights, and a per-site profile index
"""

import numpy as np

from .aa_models import AA_ORDER, AAModel
from .cseries import load_cseries

GARP = "GARP"
FYMINK = "FYMINK"

_GARP_IDX = np.array([AA_ORDER.index(a) for a in GARP], dtype=np.int64)
_FYMINK_IDX = np.array([AA_ORDER.index(a) for a in FYMINK], dtype=np.int64)


###function definition

def branch_transform(freqs, gamma_g, gamma_f):
    """Apply Eq. 1.4 to one profile frequency vector.

    freqs is a length-20 profile summing to 1. Returns a length-20 vector
    that again sums to 1: the GARP class is scaled by gamma_g, the FYMINK
    class by gamma_f and the remaining class by 1, followed by
    renormalisation via the factor mu.
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


def normalise_weights(weights):
    """Return weights divided by their sum.

    The division does not make the running total reach 1.0 exactly;
    assign_site_profiles clips a draw above the final cumulative weight onto
    the last class, so a rounding tail cannot fall through.
    """
    w = np.asarray(weights, dtype=np.float64)
    total = float(w.sum())
    if not (total > 0.0):
        raise ValueError("profile weights must sum to a positive value")
    return w / total


def assign_site_profiles(num_sites, weights, rng):
    """Draw a profile index for every site from weights.

    rng is a callable returning a uniform variate in [0, 1). Draws are
    independent across sites and are redrawn once per replicate, matching the
    per-replicate semantics of the gamma category assignment.
    """
    # Build the CDF over the positive weights only. A class with weight 0 can
    # then never be drawn, even when the weights sum to slightly less than 1 in
    # floating point (the tail is absorbed by the last positive class instead
    # of leaking onto a profile that was assigned no probability).
    w = np.asarray(weights, dtype=np.float64)
    pos = np.flatnonzero(w > 0.0)
    if pos.size == 0:
        raise ValueError("profile weights must contain a positive entry")
    cum = np.cumsum(w[pos])
    cum[-1] = 1.0
    k = pos.size
    out = np.zeros(num_sites, dtype=np.int64)
    for i in range(num_sites):
        r = rng()
        idx = int(np.searchsorted(cum, r, side="right"))
        out[i] = pos[idx if idx < k else k - 1]
    return out


def build_profile_models(cseries_name, exchangeability, gamma_by_node):
    """Build the AA models for every (branch, profile) pair.

    exchangeability is the 190-element rate vector shared by all profiles (the
    -m model's exchangeabilities, e.g. LG). gamma_by_node maps a node to its
    (gamma_g, gamma_f) pair; nodes absent from the map are homogeneous and
    share the base models. Returns (base_models, per_node_models, weights):
    length-C model lists for the raw and the branch-transformed profiles, plus
    the C-series weights.
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
    weights = weights / weights.sum()
    return base, per_node, weights
