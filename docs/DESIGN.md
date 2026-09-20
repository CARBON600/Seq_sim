# Design notes

A Python reimplementation of Seq-Gen's Monte Carlo simulation of sequence
evolution (Rambaut & Grassly 1997): nucleotide and amino-acid sequences, standard
substitution models, rate heterogeneity, site- and branch-heterogeneous models,
node-heterogeneous models, NEXUS input.

## Layout

    seq_sim_v2/
      seq_sim.py              CLI, validation, run pipeline
      manual.py               the --manual text
      core/
        nuc_models.py         JC69, HKY, F84, GTR
        aa_models.py          JTT, WAG, PAM, BLOSUM62, LG, GENERAL + their data
        models.py             model factory and per-node registry
        gamma.py              Yang (1994) discrete gamma
        rate_heterogeneity.py codon / gamma / invariable-site rate model
        cseries.py            C10-C60 empirical profile mixtures
        profile_mixture.py    GFmix site and branch heterogeneity
        evolve.py             root sampling, traversal, mutation loop
        _validators.py        shared argument checks
      tree/tree_parser.py     TNode/TTree and the Newick parser
      output/
        writers.py            PHYLIP, relaxed PHYLIP, NEXUS, FASTA
        nexus_reader.py       NEXUS preprocessing for BEAST/MrBayes trees
      utils/                  the seeded generator
      data/cseries/           C10-C60 frequencies and weights

## Decisions worth knowing

**Transition probabilities.** `set_matrix` fills a flat array with row-wise
*cumulative* P(t) for every state; `set_vector` does the same for one state.
Cumulative rows let the sampler scan without renormalising. JC69, HKY and F84
have closed forms in exp(-beta t); GTR and the AA models build Q, scale it to
mean rate 1, symmetrise as S = sqrt(pi) Q / sqrt(pi) and eigendecompose once, so
any branch length is a spectral update.

**Rate heterogeneity.** Codon-position rates, discrete gamma, continuous gamma
and invariable sites all live in one `rate_info` dict. Codon beats discrete gamma
beats continuous gamma; the invariable mask is independent of all three. An
invariable site is not simply skipped: the remaining sites are scaled by 1/(1-p)
so the mean rate stays 1, which is what Seq-Gen does and what the analytic
expectations assume.

**Node heterogeneity.** The branch leading *into* a node may carry `[f=N]` and
`[r=N]`, selecting a composition vector and a rate matrix from the registry built
out of `--alt-freqs`, `--alt-rates` and `--alt-tstv`. Index 0 is always the model
built from the command line. Internal node labels are stored, so a branch can
also be addressed by name from `--gfmix-branch-hetero`. The idea is Foster
(2004); the `(comp_idx, rmat_idx)` registry keying is this project's own, written
from the paper. P4, Foster's reference implementation of the idea, is GPL v2 and
nothing here is derived from it. An index the registry does not contain is
rejected at use time with the available keys listed, so an annotation typo cannot
silently fall back to another model.

**Site heterogeneity.** A site draws one C-series profile by weight, and the draw
is redrawn once per replicate, matching the per-replicate semantics of the gamma
category assignment. Each profile becomes its own amino-acid model, so a run with
C profiles carries C rate matrices per branch.

**Branch heterogeneity.** This changes the *composition*, not the rate. On a
designated branch the profile frequencies are rescaled (GARP by `gamma_G`, FYMINK
by `gamma_F`, the remaining ten by 1) and renormalised so the profile still sums
to one; the exchangeabilities are untouched. `mu` is exactly the normalising
constant implied by that sum-to-one requirement. Every profile on the branch uses
the same `(gamma_G, gamma_F)`. A node with both gammas equal to 1 is skipped, so
a homogeneous entry is a no-op.

Because no exchangeability is altered, a branch-heterogeneous run has the same
overall rate scale as a homogeneous one; the effect is a shift in equilibrium
composition. A length-scale difference between a simulated tree and an estimate
of it therefore comes from the estimator, not from the simulator's own
bookkeeping. `--normalize-branch-rates`, which rescales per-branch rate
multipliers to a length-weighted mean of 1, is a separate switch and is off by
default.

**Eq. 1.8.** The weighted-average constraint `sum_c w_c mu_c = 1` is deliberately
not imposed, and the source says it is absent from the maximum likelihood
implementations. `branch_transform` returns `mu` as the sum-to-one normaliser and
nothing else. `docs/REFERENCES.md` records the reasoning; a simulator that does
impose the constraint differs from this one by a scale factor that must not be
read as bias.

**Reproducibility.** Every draw goes through one module-level generator, so a
given `-z` reproduces a run byte for byte. `tests/test_golden.py` locks three
output hashes; they only move if the order or arithmetic of the draws changes.

## Known limits

- No UDM profile mixtures. The C-series profiles are supported; UDM is not.
- No indel simulation, no codon models, no MTREV/CPREV45/MTART/HIVB, no +H
  heterotachy, no Gillespie sampler (matrix exponentiation is exact for a fixed
  tree, so this is a cost choice, not a correctness one). Seq-Gen 1.3.5 does
  carry those matrices, in `source/aamodels.c`, if a port is ever wanted.
- Branch heterogeneity covers the interior of the tree. The source includes the
  root edge in its parameter vector; whether a root-edge gamma behaves as
  intended here has not been tested.
- The interaction between the C-series profile mixture and per-site rate
  variation (`-a`, `-g`, `-i`) is not exercised by the test suite.
- `-wa` with several partitions writes the first partition only unless `-os`.
- The generator in `utils` is global, so the package is not thread-safe.

## References

Full bibliography with verification tiers: `docs/REFERENCES.md`.

Algorithms: Rambaut & Grassly (1997) Seq-Gen; Jukes & Cantor (1969) JC69;
Hasegawa, Kishino & Yano (1985) HKY; Felsenstein (1984) F84; Tavare (1986) GTR;
Yang (1994) discrete gamma; Foster (2004) compositional heterogeneity and its
reference implementation P4 (https://github.com/pgfoster/p4-phylogenetics);
Kosiol & Goldman (2005) DCMUT corrections; McCarthy, Susko, Harada & Roger (2026)
site-and-branch heterogeneity.

Empirical matrices: Jones, Taylor & Thornton (1992) JTT; Whelan & Goldman (2001)
WAG; Dayhoff, Schwartz & Orcutt (1978) PAM; Henikoff & Henikoff (1992) BLOSUM62;
Le & Gascuel (2008) LG; Si Quang, Gascuel & Lartillot (2008) C10-C60.
