"""Step B1: our per-(edge, class) P(t) vs gfmixf's own frca/Sfr2WLam/Mij_nod.

qdrv2 is compiled from the unmodified gfmixff.c and calls, in main()'s order:

    mu = frca(nchar, ca, cia, fr);     # pi' = gamma*pi / sum(gamma*pi)
    Sfr2WLam(nchar, S, fr, W, Lam);    # Q from pi' and S, normalised to mean 1
    Mij_nod(nchar, t, fr, W, Lam, M);  # P(t) = exp(Q t)

So one comparison covers the composition transform, the rate-matrix
normalisation and the matrix exponential.
"""
import os
import subprocess
import sys

import numpy as np

from _env import directory

GM = directory('GFMIX_DIR', 'the directory holding qdrv2 and gfmix_lnl_align')

SRC = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "src")
sys.path.insert(0, SRC)
from seq_sim.core.aa_models import AAModel, empirical_rates
from seq_sim.core.profile_mixture import branch_transform, build_profile_models

DRIVER = GM + "/qdrv2"


def c10_profile0():
    base, per_node, weights = build_profile_models("C10", empirical_rates("LG"), {})
    return np.asarray(base[0].freq, dtype=float)


def ref_P(pi, g, f, t):
    inp = (" ".join("%.17g" % x for x in pi) + " "
           + "%.17g %.17g" % (g, f) + " %.17g" % t + "\n")
    out = subprocess.run([DRIVER], input=inp.encode(), capture_output=True,
                         check=True)
    rows = [r.split() for r in out.stdout.decode().splitlines()]
    return np.array([[float(x) for x in r] for r in rows])


def our_P(pi, g, f, t):
    freqs = branch_transform(np.array(pi, dtype=float), g, f)
    model = AAModel(freqs, empirical_rates("LG"))
    buf = np.zeros(400)
    model.set_matrix(buf, t)
    cum = buf.reshape(20, 20)
    P = np.empty((20, 20))
    P[:, 0] = cum[:, 0]
    P[:, 1:] = np.diff(cum, axis=1)
    return P, freqs


def main():
    pi0 = c10_profile0()
    print("C10 profile 0, sum = %.12f" % pi0.sum())
    print()
    print("  gammas          t     max|ref-ours|   ref row-sum dev   rel frobenius")
    worst = 0.0
    for (g, f) in [(1.0, 1.0), (3.0, 0.4), (8.0, 0.1), (0.3, 2.5)]:
        for t in (0.1, 0.5, 2.0):
            R = ref_P(pi0, g, f, t)
            O, _ = our_P(pi0, g, f, t)
            d = np.abs(R - O).max()
            worst = max(worst, d)
            rs = np.abs(R.sum(axis=1) - 1.0).max()
            rel = np.linalg.norm(R - O) / np.linalg.norm(R)
            print("  (%5.2f,%5.2f)  %5.2f    %.3e        %.2e          %.3e"
                  % (g, f, t, d, rs, rel))
    print()
    print("worst |ref - ours| over the grid = %.3e" % worst)


if __name__ == "__main__":
    main()