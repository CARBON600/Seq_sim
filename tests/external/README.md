# External cross-checks (run by hand; pytest does not collect them)

These scripts answer a question the in-repo tests cannot: does the simulator
sample from the *same* distribution that the reference programs fit? They need
external software, so they are never run by `pytest` (file names start with
`check_`, which pytest does not collect).

## Environment

No directory is guessed. Each script reads what it needs from the environment
and stops with a message naming the variable if it is unset.

| variable | meaning | used by |
|---|---|---|
| `SEQ_SIM_EXT_WORK` | writable dir; a `b2/` subdir is created for inputs | B1, B2, C1, C2, `p4_*_side.py` |
| `GFMIX_DIR` | holds `qdrv2` (built from `gfmix_qdrv.c` + the original `gfmixff.c`) and `gfmix_lnl_align` | B1, B2 |
| `P4_DIR` | the built p4-phylogenetics checkout | C1, C2 |

    export SEQ_SIM_EXT_WORK=$HOME/seq-sim-ext
    export GFMIX_DIR=$HOME/gfmix
    export P4_DIR=$HOME/p4

## Scripts

| script | compares | passes when |
|---|---|---|
| `check_A_sampler.py` | simulated pattern frequencies vs the exact pattern distribution implied by the same P matrices (all 20^4 patterns) | chi-square p > 0.01, and two negative controls blow up |
| `check_B1_gfmix_matrices.py` | our `P(t)` vs `frca`+`Sfr2WLam`+`Mij_nod` from the unmodified `gfmixff.c` | max abs diff < 1e-12 |
| `check_B2_gfmix_likelihood.py` | our pruning lnL vs `gfmix_lnl_align -n 0` (fixed parameters, no optimisation) | |diff|/site < 1e-8 after weights are normalised |
| `check_C1_p4_hky.py` | our HKY(`-t` ts/tv) vs p4 `newRMatrix(spec='2p', val=kappa)` | |diff| < 1e-6 on 1000 sites |
| `check_C2_p4_clade.py` | p4 `setModelComponentOnNode(clade=1)` vs our per-edge composition | exactly one edge subset matches |

## Measured results (pre-release tree, 2026-11)

* A: X2 = 1697.6, df = 1742, p = 0.773; controls X2 = 26006 (no root gamma) and 5885 (no branch gamma).
* B1: worst |ref - ours| = 2.132e-14 over 4 gamma pairs x 3 branch lengths.
* B2: -1.295e-09 per site, +5.1e-12 once the weights are normalised; the residual
  equals ln(sum of C10 weights) = ln(1.0000000013), i.e. the shipped weight table
  does not sum to exactly one. The normalised column is the one the current tree
  reproduces. `-m` row order was identified by a probe: rows 0-3 = edges into
  A, B, C, D; row 4 = root->H2; row 5 = root->H1; row 6 = root prior.
* C1: three (pi, kappa) pairs, |diff| <= 2.7e-12. Same tstv=2.0 needs kappa=4.0 at
  equal frequencies but kappa=6.25 at (0.4,0.1,0.1,0.4) -- the mapping is not a
  constant, so the comparison is built from Q; the model name is not enough.
* C2: p4 `clade=1` on node H1 equals our transformed composition applied to
  H1 + A + B (diff 6.7e-09); applying it to H1 alone is off by 303 nats.
  seq_sim's `[f=N]` does NOT inherit to descendants.

## Not covered

* GFmix `+F` (a free frequency class fitted internally by gfmixf) and `+I`
  (absent from the GFmix paper/handbook): not implemented.
* p4 `newComp(spec='empirical')` path; only `spec='specified'` is exercised.
* The p4 snapshot here is the old API (`calcLogLike`, `Tree.newComp`); p4 1.x is not exercised.

## Building the GFmix driver

    gcc -O0 -fopenmp -ffunction-sections -fdata-sections \
        -o qdrv2 gfmix_qdrv.c gfmixff.c -Wl,--gc-sections -lm

Build the driver from source: an object file produced by a different toolchain
will not link, and compiling here also drops the Fortran-only helpers.
