# Design notes

`seq_sim_v2` is a Python implementation of Monte Carlo molecular sequence
evolution along phylogenetic trees, following the simulation architecture of
Seq-Gen (Rambaut & Grassly 1997) and extending it with node-heterogeneous
composition/rate registries (Foster 2004) and GFmix site-and-branch
heterogeneity (McCarthy et al. 2026).

## Package layout

    seq_sim_v2/
      seq_sim.py              CLI, argument validation, simulation pipeline
      manual.py               --manual documentation text
      core/
        nuc_models.py         JC69, HKY, F84, GTR
        aa_models.py          JTT, WAG, PAM, BLOSUM62, LG, GENERAL
        models.py             model factory and (comp_idx, rmat_idx) registry
        gamma.py              Yang (1994) discrete and continuous gamma
        rate_heterogeneity.py codon, gamma, and invariable-site rate setup
        cseries.py            C10-C60 empirical profile mixtures
        profile_mixture.py    GFmix site and branch heterogeneity
        evolve.py             root sampling, pre-order traversal, mutation loop
        _validators.py        shared parameter validation helpers
      tree/tree_parser.py     TNode/TTree structures and Newick parser
      output/
        writers.py            PHYLIP, relaxed PHYLIP, NEXUS, FASTA writers
        nexus_reader.py       NEXUS tree block extractor
      utils/                  module-level seeded PRNG wrapper
      data/cseries/           C10-C60 frequency and weight tables

## Key implementation details

**Transition probabilities.** `set_matrix` fills a flat array with row-wise
cumulative transition probabilities P(t) for all states, while `set_vector`
computes the cumulative row for a single parent state. Using cumulative rows
allows the per-site sampler to draw child states with a linear scan. JC69, HKY,
and F84 use closed-form analytic expressions; GTR and amino-acid models build
Q = R diag(pi), scale it to mean rate 1, symmetrise it as
S = sqrt(pi) Q / sqrt(pi), and eigendecompose S once per model instance. The
symmetrisation is checked to 1e-12, and a violation raises FloatingPointError
rather than producing a plausible wrong answer.

**Rate heterogeneity.** Codon-position rates, discrete gamma, continuous gamma,
and invariable sites are configured in a single `rate_info` dictionary.
Precedence follows codon rates > discrete gamma > continuous gamma, with the
invariable-site mask applied independently. When invariable sites (`-i p`) are
enabled, variable sites are scaled by 1 / (1 - p) so the overall mean
substitution rate across all sites remains 1.

**Node heterogeneity.** A branch leading into a node can carry `[f=N]` and
`[r=N]` annotations selecting a composition vector and an exchangeability matrix
from the `(comp_idx, rmat_idx)` registry built from `--alt-freqs`,
`--alt-rates`, and `--alt-tstv`. Index 0 corresponds to the primary CLI model.
Unregistered index pairs raise a `ValueError` listing the available registry
keys, so an annotation typo cannot silently fall back to another model. The idea
is Foster (2004); the `(comp_idx, rmat_idx)` keying and the annotation syntax are
this project's own, written from the paper. P4, Foster's reference
implementation of the idea, is GPL v2 and nothing here is derived from it.

**Site and branch heterogeneity (GFmix).** Under `--gfmix-cseries`, each site
draws one C-series profile according to the mixture weights at the start of each
replicate. Under `--gfmix-branch-hetero`, designated branches rescale each
profile's equilibrium frequencies by gamma_G on GARP (`G, A, R, P`), gamma_F on
FYMINK (`F, Y, M, I, N, K`), and 1 on the remaining ten amino acids, then
renormalise the profile to sum to 1 (McCarthy et al. 2026, Eqs. 1.4 and 1.7).
Because only the equilibrium composition changes while the exchangeabilities
remain fixed and each profile matrix is normalised to mean rate 1, branch
lengths keep their standard substitutions-per-site scale.

**Eq. 1.8.** The weighted-average constraint sum_c w_c mu_c = 1 is deliberately
not imposed, matching the maximum likelihood implementations (see
`docs/REFERENCES.md`). `branch_transform` returns mu as the sum-to-one
normaliser and nothing else. A simulator that does impose the constraint
differs from this one by a scale factor that must not be read as bias.

**Reproducibility.** All random draws route through the module-level generator
in `utils`, ensuring that a fixed `-z` seed and command line reproduce identical
output bytes. `tests/test_golden.py` checks SHA-256 digests across
representative nucleotide and protein runs.

## Current limitations

- Supports C-series empirical mixtures (`C10`-`C60`), but not UDM mixtures.
- Does not simulate indels, codon-level 61-state matrices, covarion (`+H`)
  heterotachy, or mitochondrial/viral empirical matrices (`MTREV`, `CPREV45`,
  `MTART`, `HIVB`). Seq-Gen 1.3.5 does carry those matrices, in
  `source/aamodels.c`, if a port is ever wanted.
- No Gillespie sampler: matrix exponentiation is exact for a fixed tree, so
  this is a cost choice, not a correctness choice.
- Branch heterogeneity covers the interior of the tree. The source includes the
  root edge in its parameter vector; whether a root-edge gamma behaves as
  intended here has not been tested.
- The interaction between the C-series profile mixture and per-site rate
  variation (`-a`, `-g`, `-i`) is not exercised by the test suite.
- `-wa` (ancestral sequence output) with multi-partition trees writes only the
  first partition unless `-os` (separate files) is active.
- Because the PRNG in `utils` is module-global, concurrent simulations should
  run in separate processes rather than threads.
