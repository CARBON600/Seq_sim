"""Site-rate variation: discrete gamma on top of the C-series mixture."""

import numpy as np
from scipy.stats import chi2

from seq_sim.core.aa_models import empirical_rates
from seq_sim.core.evolve import evolve_sequences
from seq_sim.core.gamma import discrete_gamma
from seq_sim.core.profile_mixture import assign_site_profiles, build_profile_models
from seq_sim.core.rate_heterogeneity import assign_categories, setup_rates
from seq_sim.tree.tree_parser import allocate_sequences, parse_trees
from seq_sim.utils import random as rnd
from seq_sim.utils import set_seed

NEWICK = "((A:0.10,B:0.15)H1:0.20,(C:0.12,D:0.08)H2:0.20)ROOT;"
ALPHA, K, NSITES, SEED = 0.5, 4, 6000, 424242
GAMMA_H1 = (3.0, 0.4)


def _raw(model, t):
    """Un-cumulate the row-cumulative probabilities set_matrix stores."""
    n = model.num_states
    buf = np.zeros(n * n)
    model.set_matrix(buf, t)
    m = buf.reshape(n, n)
    out = np.empty((n, n))
    out[:, 0] = m[:, 0]
    out[:, 1:] = np.diff(m, axis=1)
    return out


def test_cseries_with_discrete_gamma_matches_exact_distribution():
    set_seed(SEED)
    tree = parse_trees(NEWICK)[0][2]
    base, per_node, weights = build_profile_models(
        "C10", empirical_rates("LG"), {tree.root.branch1: GAMMA_H1})
    _, rates = discrete_gamma(ALPHA, K)
    tips = {n.name: n for n in tree.node_list if n.tip_no >= 0}
    h1, h2 = tree.root.branch1, tree.root.branch2

    exact = np.zeros((20, 20, 20, 20))
    for c in range(len(base)):
        for k in range(K):
            def mat(ch):
                model = per_node.get(ch, base)[c]
                return _raw(model, ch.length0 * ch.branch_rate * rates[k])
            u = np.einsum("rh,ha,hb->rab", mat(h1), mat(tips["A"]),
                          mat(tips["B"]))
            v = np.einsum("rh,hc,hd->rcd", mat(h2), mat(tips["C"]),
                          mat(tips["D"]))
            prior = per_node.get(tree.root, base)[c].freq
            exact += (weights[c] / K) * np.einsum("r,rab,rcd->abcd",
                                                 prior, u, v)
    assert abs(exact.sum() - weights.sum()) < 1e-9

    site_profile = assign_site_profiles(NSITES, weights, rnd)
    rate_info = setup_rates(NSITES, gamma_shape=ALPHA, gamma_cats=K)
    assign_categories(rate_info, NSITES, ALPHA, K)
    allocate_sequences(tree, NSITES)
    evolve_sequences(tree, 0, NSITES, 1.0, {}, rate_info, None,
                     {"base_models": base, "per_node": per_node,
                      "site_profile": site_profile})

    flat = np.zeros(NSITES, dtype=np.int64)
    for name in "ABCD":
        flat = flat * 20 + np.asarray(tips[name].sequence, dtype=np.int64)
    counts = np.bincount(flat, minlength=20 ** 4).astype(float)
    exp = (exact * NSITES).reshape(-1)
    # Cells below the expected-count floor are pooled into one bin, so the
    # kept bins still sum to the sample size and m bins have m - 1 dof.
    keep = exp >= 5.0
    obs = np.append(counts[keep], counts[~keep].sum())
    expd = np.append(exp[keep], exp[~keep].sum())
    stat = float(((obs - expd) ** 2 / expd).sum())
    dof = int(obs.size) - 1
    p = float(chi2.sf(stat, dof))
    assert dof > 20, dof
    assert p > 0.001, (stat, dof, p)
