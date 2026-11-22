"""Step A: does the GFmix-style simulator sample from its own model?

Exact side
    Contraction over the whole pattern space (20^4 = 160000 patterns) using the
    very P matrices the sampler draws from:

        P(x_A,x_B,x_C,x_D) = sum_c w_c sum_k (1/K) sum_root pi_c(root)
                             U(root,x_A,x_B) * V(root,x_C,x_D)
        U(r,a,b) = sum_h P(root->h) P(h->a) P(h->b)
        V(r,c,d) = sum_h P(root->h) P(h->c) P(h->d)

Sample side
    evolve_sequences() on the same tree with the same rate_info / profile_info,
    then count tip patterns.

Both sides read the same matrices, so a pass says the sampler and the data
plumbing are faithful to the model the code defines.  It says nothing about
whether those matrices are the right ones: that is Step B (gfmixf) and Step C
(p4).

Two negative controls are run at the end: the same comparison with the root-edge
gamma dropped from the simulated side, and then with the branch gamma dropped.
If the test had no power, both controls would also pass.
"""
import os
import sys

import numpy as np
from scipy.stats import chi2

SRC = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "src")
sys.path.insert(0, SRC)

from seq_sim.core.aa_models import empirical_rates
from seq_sim.core.evolve import evolve_sequences
from seq_sim.core.gamma import discrete_gamma
from seq_sim.core.profile_mixture import assign_site_profiles, build_profile_models
from seq_sim.core.rate_heterogeneity import assign_categories, setup_rates
from seq_sim.tree.tree_parser import allocate_sequences, parse_trees
from seq_sim.utils import random as rnd
from seq_sim.utils import set_seed

NEWICK = "((A:0.10,B:0.15)H1:0.20,(C:0.12,D:0.08)H2:0.20)ROOT;"
CSERIES = "C10"
GAMMA_H1 = (3.0, 0.4)
GAMMA_ROOT = (2.0, 0.5)
ALPHA = 0.5
K = 2
NSITES = 200_000
SEED = 20261121
AA = "ARNDCQEGHILKMFPSTWYV"


def tip_map(tree):
    return {node.name: node for node in tree.node_list if node.tip_no >= 0}


def raw_matrix(model, t):
    """Recover the un-cumulated P(t) from set_matrix's cumulative rows."""
    n = model.num_states
    buf = np.zeros(n * n)
    model.set_matrix(buf, t)
    m = buf.reshape(n, n)
    out = np.empty((n, n))
    out[:, 0] = m[:, 0]
    out[:, 1:] = np.diff(m, axis=1)
    return out


def exact_distribution(tree, per_node, base_models, weights, cat_rates):
    tips = tip_map(tree)
    root, h1, h2 = tree.root, tree.root.branch1, tree.root.branch2
    total = np.zeros((20, 20, 20, 20))
    for c in range(len(base_models)):
        for k in range(K):
            t = cat_rates[k]

            def mat(child):
                model = per_node.get(child, base_models)[c]
                return raw_matrix(model, child.length0 * child.branch_rate * t)

            U = np.einsum("rh,ha,hb->rab", mat(h1), mat(tips["A"]),
                          mat(tips["B"]))
            V = np.einsum("rh,hc,hd->rcd", mat(h2), mat(tips["C"]),
                          mat(tips["D"]))
            prior = per_node.get(root, base_models)[c].freq
            total += (weights[c] / K) * np.einsum("r,rab,rcd->abcd",
                                                 prior, U, V)
    return total


def simulate(tree, per_node, base_models, weights, cat_rates, n_sites):
    site_profile = assign_site_profiles(n_sites, weights, rnd)
    rate_info = setup_rates(n_sites, gamma_shape=ALPHA, gamma_cats=K)
    assign_categories(rate_info, n_sites, ALPHA, K)
    allocate_sequences(tree, n_sites)
    evolve_sequences(tree, 0, n_sites, 1.0, {}, rate_info, None,
                     {"base_models": base_models, "per_node": per_node,
                      "site_profile": site_profile})
    tips = tip_map(tree)
    return [np.asarray(tips[name].sequence, dtype=np.int64) for name in "ABCD"]


def compare(name, exact, seqs, n_sites):
    flat = ((seqs[0] * 20 + seqs[1]) * 20 + seqs[2]) * 20 + seqs[3]
    counts = np.bincount(flat, minlength=20 ** 4).astype(np.float64)
    obs = counts.reshape(20, 20, 20, 20) / n_sites
    diff = np.abs(obs - exact)
    where = np.unravel_index(diff.argmax(), diff.shape)
    exp4 = exact * n_sites
    # Pool the cells below the expected-count floor into one bin, so the bins
    # still sum to the sample size and m bins have m - 1 dof.
    keep = (exp4 >= 5.0).reshape(-1)
    o = np.append(counts[keep], counts[~keep].sum())
    e = np.append(exp4.reshape(-1)[keep], exp4.reshape(-1)[~keep].sum())
    stat = float(((o - e) ** 2 / e).sum())
    dof = int(o.size) - 1
    p = chi2.sf(stat, dof)
    print("  %-18s max|d|=%.6f at %-4s  mean|d|=%.7f  X2=%9.1f dof=%4d  p=%.3g"
          % (name, diff.max(), "".join(AA[i] for i in where), diff.mean(),
             stat, dof, p))
    return stat, p


def main():
    set_seed(SEED)
    tree = parse_trees(NEWICK)[0][2]
    assert tree.rooted, "test tree must be rooted (bifurcating at the root)"
    base_models, per_node, weights = build_profile_models(
        CSERIES, empirical_rates("LG"),
        {tree.root: GAMMA_ROOT, tree.root.branch1: GAMMA_H1})
    _, cat_rates = discrete_gamma(ALPHA, K)

    print("tree           : %s" % NEWICK)
    print("C-series       : %s (%d profiles)" % (CSERIES, len(base_models)))
    print("branch gamma   : edge into H1 = %s ; root edge = %s" % (GAMMA_H1, GAMMA_ROOT))
    print("discrete gamma : alpha=%s K=%d -> rates %s" % (ALPHA, K, np.round(cat_rates, 6)))
    print("sites          : %d, seed %d" % (NSITES, SEED))
    print()

    wsum = float(np.sum(weights))
    print("C-series weight sum = %.13f" % wsum)
    wn = weights / wsum
    print("  (the sampler's CDF is pinned at 1.0, so it uses w/sum(w); the exact")
    print("   side below uses the same normalisation; the sum is printed as a check)")
    print()

    exact = exact_distribution(tree, per_node, base_models, wn, cat_rates)
    mass = float(exact.sum())
    print("self-check  sum of all 20^4 probabilities = %.15f  (|dev| = %.2e)"
          % (mass, abs(mass - 1.0)))
    print()

    print("MAIN RUN")
    seqs = simulate(tree, per_node, base_models, wn, cat_rates, NSITES)
    compare("simulator", exact, seqs, NSITES)

    print()
    print("NEGATIVE CONTROL (root-edge gamma dropped on the simulated side)")
    set_seed(SEED)
    tree2 = parse_trees(NEWICK)[0][2]
    base2, per_node2, _ = build_profile_models(
        CSERIES, empirical_rates("LG"),
        {tree2.root: (1.0, 1.0), tree2.root.branch1: GAMMA_H1})
    seqs2 = simulate(tree2, per_node2, base2, wn, cat_rates, NSITES)
    compare("no root gamma", exact, seqs2, NSITES)

    print()
    print("NEGATIVE CONTROL (branch gamma dropped on the simulated side)")
    set_seed(SEED)
    tree3 = parse_trees(NEWICK)[0][2]
    base3, per_node3, _ = build_profile_models(
        CSERIES, empirical_rates("LG"),
        {tree3.root: GAMMA_ROOT, tree3.root.branch1: (1.0, 1.0)})
    seqs3 = simulate(tree3, per_node3, base3, wn, cat_rates, NSITES)
    compare("no branch gamma", exact, seqs3, NSITES)

    print()
    print("tip marginals (main run):")
    flat = ((seqs[0] * 20 + seqs[1]) * 20 + seqs[2]) * 20 + seqs[3]
    counts = np.bincount(flat, minlength=20 ** 4).astype(np.float64)
    obs = counts.reshape(20, 20, 20, 20) / NSITES
    for j, name in enumerate("ABCD"):
        axes = tuple(i for i in range(4) if i != j)
        ex_m, ob_m = exact.sum(axis=axes), obs.sum(axis=axes)
        se = np.sqrt(np.maximum(ex_m * (1.0 - ex_m) / NSITES, 1e-300))
        print("  %s  max|obs-exact| = %.6f   max|z| = %.2f"
              % (name, np.abs(ob_m - ex_m).max(), np.abs((ob_m - ex_m) / se).max()))

    print()
    print("10 most probable patterns (main run):")
    print("  pat   exact      observed    z")
    flat_ex = exact.reshape(-1)
    for p in np.argsort(flat_ex)[-10:][::-1]:
        a, b, c, d = np.unravel_index(p, (20, 20, 20, 20))
        pat = AA[a] + AA[b] + AA[c] + AA[d]
        e_p = float(flat_ex[p])
        o_p = counts[p] / NSITES
        z = (o_p - e_p) / np.sqrt(max(e_p * (1.0 - e_p) / NSITES, 1e-300))
        print("  %s  %.6f   %.6f   %+.2f" % (pat, e_p, o_p, z))


if __name__ == "__main__":
    main()