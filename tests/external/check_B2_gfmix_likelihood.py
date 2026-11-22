"""Step B2: our site-pattern likelihood vs gfmixf's, at fixed parameters.

gfmixf runs with -n 0, so it evaluates the likelihood without optimising.
Both sides get the same tree, C10 table, weights, alpha and per-edge
GARP/FYMINK multipliers.

  A  no branch/root gamma                -> mixture weights, LG, C10 rows, gamma
  B  every branch (3.0, 0.4), root (2.0, 0.5)  -> adds the transform end to end
  C  row 0 = (5.0, 1.0), rest 1          -> identifies which edge row 0 is
  D  same config at 4 sites and at 1000 sites -> is the residual accumulation?
"""
import os
import subprocess
import sys

import numpy as np

from _env import directory

W = directory('SEQ_SIM_EXT_WORK', 'a writable working directory')
GM = directory('GFMIX_DIR', 'the directory holding qdrv2 and gfmix_lnl_align')

SRC = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "src")
sys.path.insert(0, SRC)
from seq_sim.core.aa_models import AAModel, empirical_rates
from seq_sim.core.gamma import discrete_gamma
from seq_sim.core.profile_mixture import branch_transform, build_profile_models
from seq_sim.tree.tree_parser import parse_trees

G = GM
DRV = G + "/qdrv2"
BIN = G + "/gfmix_lnl_align"
D = os.path.join(W, 'b2')
NEWICK = "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);"
ALPHA = 0.5
K = 4
PATTERNS = [("A", "A", "A", "A"), ("A", "G", "A", "G"),
            ("G", "G", "G", "G"), ("A", "R", "N", "D")]
AA = "ARNDCQEGHILKMFPSTWYV"


def write_seq(reps):
    seqs = [[], [], [], []]
    for pat in PATTERNS:
        for _ in range(reps):
            for t in range(4):
                seqs[t].append(pat[t])
    with open(D + "/seq.phy", "w") as fh:
        fh.write("4 %d\n" % (len(PATTERNS) * reps))
        for t, name in enumerate("ABCD"):
            fh.write("%-10s%s\n" % (name, "".join(seqs[t])))


def setup():
    os.makedirs(D, exist_ok=True)
    with open(D + "/tree.nwk", "w") as fh:
        fh.write(NEWICK + "\n")
    base, per_node, weights = build_profile_models("C10", empirical_rates("LG"), {})
    with open(D + "/freq.dat", "w") as fh:
        for m in base:
            fh.write(" ".join("%.17g" % x for x in m.freq) + "\n")
    with open(D + "/wt.dat", "w") as fh:
        fh.write(" ".join("%.17g" % x for x in weights) + "\n")
    return base, weights


def write_mfile(path, rows):
    with open(path, "w") as fh:
        for g in rows:
            fh.write("%.17g %.17g\n" % g)


def run_ref(mfile):
    args = [BIN, "-s", D + "/seq.phy", "-t", D + "/tree.nwk", "-c", "10",
            "-f", D + "/freq.dat", "-w", D + "/wt.dat", "-e", "9",
            "-a", repr(ALPHA), "-n", "0", "-m", mfile]
    out = subprocess.run(args, capture_output=True, check=True)
    return float(out.stdout.decode().split()[0])


def model_for(pi, g):
    return AAModel(branch_transform(np.asarray(pi, dtype=float), g[0], g[1]),
                   empirical_rates("LG"))


def rawP(model, t):
    buf = np.zeros(400)
    model.set_matrix(buf, t)
    cum = buf.reshape(20, 20)
    P = np.empty((20, 20))
    P[:, 0] = cum[:, 0]
    P[:, 1:] = np.diff(cum, axis=1)
    return P


def our_lnl(tree, base, weights, branch_g, root_g, reps, target=None, target_g=None):
    root = tree.root
    _, cat_rates = discrete_gamma(ALPHA, K)
    cache = {}

    def model_of(c, g):
        key = (c, g)
        if key not in cache:
            cache[key] = model_for(base[c].freq, g)
        return cache[key]

    root_models = [model_of(c, root_g) for c in range(len(base))]
    total = 0.0
    for pat in PATTERNS:
        obs = {name: AA.index(ch) for name, ch in zip("ABCD", pat)}
        lik = 0.0
        for c in range(len(base)):
            for k in range(K):
                def vec(node):
                    if node.tip_no >= 0:
                        v = np.zeros(20)
                        v[obs[node.name]] = 1.0
                        return v
                    L = np.ones(20)
                    for ch in (node.branch1, node.branch2):
                        g = target_g if (target is not None and ch is target) else branch_g
                        L = L * (rawP(model_of(c, g), ch.length0 * ch.branch_rate * cat_rates[k])
                                 @ vec(ch))
                    return L
                lik += (weights[c] / K) * float(np.asarray(root_models[c].freq) @ vec(root))
        total += reps * np.log(lik)
    return total


def main():
    base, weights = setup()
    tree = parse_trees(NEWICK)[0][2]

    _, cat_rates = discrete_gamma(ALPHA, K)
    rr = [float(l.split()[0]) for l in subprocess.run(
        [DRV, "rates"], input=("%r %d" % (ALPHA, K)).encode(),
        capture_output=True, check=True).stdout.decode().splitlines()]
    print("gamma rates ours   : %s" % np.round(cat_rates, 10))
    print("gamma rates gfmixf : %s" % np.round(rr, 10))
    print("max |diff| = %.3e" % np.abs(np.array(rr) - cat_rates).max())
    print()

    ones = (1.0, 1.0)
    write_mfile(D + "/m_ones.dat", [ones] * 7)
    write_mfile(D + "/m_unif.dat", [(3.0, 0.4)] * 6 + [(2.0, 0.5)])

    write_seq(250)
    print("  config (1000 sites)                        ref lnL            our lnL          diff")
    for tag, mf, bg, rg in (("A no gamma", D + "/m_ones.dat", ones, ones),
                            ("B branch+root gamma", D + "/m_unif.dat", (3.0, 0.4), (2.0, 0.5))):
        ref = run_ref(mf)
        ours = our_lnl(tree, base, weights, bg, rg, 250)
        print("  %-38s %.10f  %.10f  %+.2e" % (tag, ref, ours, ref - ours))

    print()
    print("  C probe: row 0 = (5.0, 1.0), other rows = 1")
    write_mfile(D + "/m_probe.dat", [(5.0, 1.0)] + [ones] * 6)
    ref = run_ref(D + "/m_probe.dat")
    print("    ref lnL = %.10f" % ref)
    cands = [tree.root.branch1, tree.root.branch2,
             tree.tips[0], tree.tips[1], tree.tips[2], tree.tips[3]]
    labels = ["root->H1", "root->H2", "root->A", "root->B", "root->C", "root->D"]
    best = None
    for lab, ch in zip(labels, cands):
        val = our_lnl(tree, base, weights, ones, ones, 250, ch, (5.0, 1.0))
        d = abs(val - ref)
        print("    if row0 = %-8s our lnL = %.10f   diff = %+.3e" % (lab, val, val - ref))
        if best is None or d < best[0]:
            best = (d, lab)
    print("    best match: %s (|diff| %.2e)" % (best[1], best[0]))

    print()
    print("  D residual scaling")
    for reps in (1, 10, 100, 250):
        write_seq(reps)
        a = run_ref(D + "/m_ones.dat")
        b = our_lnl(tree, base, weights, ones, ones, reps)
        print("    %4d sites  lnL=%14.6f  diff=%+.3e   diff/site=%+.3e"
              % (4 * reps, a, a - b, (a - b) / (4 * reps)))


if __name__ == "__main__":
    main()