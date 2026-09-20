# seq_sim_v2

Monte Carlo simulation of molecular sequence evolution along phylogenetic trees.

Status: 1.4.1, research release. No CI, no package index; install from a checkout.

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [What it does](#what-it-does)
- [Site-and-branch heterogeneity (GFmix)](#site-and-branch-heterogeneity-gfmix)
- [Command-line reference](#command-line-reference)
- [Input formats](#input-formats)
- [Output formats](#output-formats)
- [Data provenance](#data-provenance)
- [Reproducibility](#reproducibility)
- [Tests](#tests)
- [Verification against reference implementations](#verification-against-reference-implementations)
- [Citing this software](#citing-this-software)
- [Contributing](#contributing)
- [License](#license)

## Installation

Requires Python >= 3.9, NumPy >= 1.21, SciPy >= 1.7.

```bash
pip install .
# or, for development (adds pytest):
pip install -e ".[test]"
```

This installs the `seq-sim` console script and the `python -m seq_sim_v2` entry
point.

## Quick start

`-m` is required; the tree is the positional argument. Output goes to a directory
named `output/` unless `-y` is given.

```bash
# 1000 bp under HKY with ts/tv = 4.0, FASTA (default)
python -m seq_sim_v2 -m HKY -l 1000 -t 4.0 tree.nwk

# 500 amino-acid sites under WAG
python -m seq_sim_v2 -m WAG -l 500 tree.nwk

# GTR with four discrete gamma categories, shape 0.5
python -m seq_sim_v2 -m GTR -a 0.5 -g 4 tree.nwk

# Normalise total tree height to 0.5 substitutions/site (needs a rooted tree)
python -m seq_sim_v2 -m HKY -d 0.5 tree.nwk

# C10 profile mixture with branch heterogeneity, 3 replicates, PHYLIP
python -m seq_sim_v2 -m LG -l 500 -n 3 --gfmix-cseries C10 \
    --gfmix-branch-hetero gammas.txt -o p -y run1 tree.nwk

# Full manual
python -m seq_sim_v2 --manual
```

## What it does

seq_sim_v2 evolves nucleotide or amino-acid sequences down a phylogenetic tree
and writes the result in a standard alignment format. Every draw comes from one
seeded generator, so a given `-z` reproduces a run byte for byte.

Substitution models are JC69, HKY, F84 and GTR for nucleotides and JTT, WAG, PAM,
BLOSUM (BLOSUM62), LG and GENERAL for amino acids; the five named empirical
matrices agree with Seq-Gen 1.3.5 entry by entry and in the same order, to
`max|difference| = 0.000e+00`. Per-site rate variation covers discrete gamma
(`-a` with `-g`), continuous gamma (`-a` alone), codon-position rates (`-c`) and
invariable sites (`-i`), with codon rates taking precedence over discrete gamma,
which takes precedence over continuous gamma, and the invariable mask independent
of all three. Site heterogeneity means a C10-C60 profile mixture
(`--gfmix-cseries`) with one profile drawn per site by weight, and branch
heterogeneity means per-branch `gamma_G` / `gamma_F` multipliers on amino-acid
profiles (`--gfmix-branch-hetero`). Per-node models come from `[f=N]` and `[r=N]`
annotations on branches, which select a composition vector and a rate matrix from
a registry built out of `--alt-freqs`, `--alt-rates` and `--alt-tstv` (Foster
2004; see Acknowledgements).

Input is Newick, or a NEXUS tree from BEAST, MrBayes or IQ-TREE, with partitioned
data, an inserted text block per dataset, and an optional supplied ancestral root
sequence. Output is PHYLIP (strict and relaxed), NEXUS or FASTA, either one file
per dataset or one combined file.

## Site-and-branch heterogeneity (GFmix)

This is the model of McCarthy et al. (2026). It has two independent axes.

**Site heterogeneity.** Each site draws one profile from the C-series mixture
(Si Quang et al. 2008) with probability equal to that profile's weight. With `C`
profiles there are `C` amino-acid models per branch, and the site-to-profile
assignment is redrawn once per replicate.

**Branch heterogeneity.** On designated branches the profile frequencies are
rescaled before the rate matrix is built. Write `pi^(c)` for a profile, `G` for
the GARP set (G, A, R, P), `F` for FYMINK (F, Y, M, I, N, K) and `O` for the
remaining ten amino acids. The rescaled profile is

```
pi^(ce)_j = mu^(ce) * gamma_G * pi^(c)_j     j in G
          = mu^(ce) * gamma_F * pi^(c)_j     j in F
          = mu^(ce) * pi^(c)_j               otherwise
```

where `mu^(ce)` is fixed by the requirement that the profile still sums to one:

```
mu^(ce) = 1 / ( gamma_G * sum_G pi + gamma_F * sum_F pi + sum_O pi )
```

Two properties matter when reading results. The rescaled vector is a
**composition**, not a rate multiplier: the exchangeabilities are untouched, so a
transformed branch has the same overall rate scale as any other branch. And
`gamma_G` / `gamma_F` are shared by every profile on that branch.

The branch file takes one line per branch, `name gamma_G gamma_F`. Branches that
are not listed are homogeneous.

```
taxon_02  2.0  0.5
het0      0.5  2.0
```

Names resolve against the tree, so an internal node works if the tree labels it.
A node with both gammas equal to 1 is skipped, which makes a homogeneous entry a
no-op.

```bash
# site heterogeneity only
python -m seq_sim_v2 -m LG -l 1000 --gfmix-cseries C30 -o p -y c30 tree.nwk

# custom weights plus a per-site assignment log
python -m seq_sim_v2 -m LG -l 1000 --gfmix-cseries C60 \
    --gfmix-weights "$(cat weights60.txt)" --gfmix-log -y c60log tree.nwk

# site and branch heterogeneity together
python -m seq_sim_v2 -m LG -l 1000 --gfmix-cseries C10 \
    --gfmix-branch-hetero gammas.txt -o p -y both tree.nwk
```

`--gfmix-cseries` requires an amino-acid model and is mutually exclusive with
`--alt-freqs`, `--alt-rates` and `--alt-tstv`. `--gfmix-weights` overrides the
built-in weights and must sum to 1. `--gfmix-log` writes the per-site profile
assignment, the weights and the branch gammas to `<prefix>.gfmix.log`.

Equation numbers (1.4, 1.7) follow the preprint version of McCarthy et al., which
is the version this code was written against; see
[docs/REFERENCES.md](docs/REFERENCES.md).

## Command-line reference

`python -m seq_sim_v2 --manual` prints the full manual, including worked
examples. The table below is the short form.

**Model and tree**

| Option | Default | Meaning |
|---|---|---|
| `treefile` | stdin | Newick tree file (positional) |
| `-m, --model` | required | JC69, HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, LG, GENERAL |
| `-l, --length` | 1000 | Sequence length in sites |
| `-n, --num-datasets` | 1 | Number of datasets per tree |
| `-s, --branch-scale` | 0.0 | Branch length scaling factor |
| `-d, --tree-scale` | 0.0 | Normalise total tree height (root to tip) to N substitutions/site; needs a rooted tree; mutually exclusive with `-s` |
| `-br, --branch-rates` | off | Treat `[value]` after a branch length as a per-branch rate multiplier |
| `--normalize-branch-rates` | off | Scale branch rate multipliers to weighted mean 1.0 |
| `-z, --seed` | none | Random number seed |
| `-q, --quiet` | off | Suppress informational output |

**Per-node model registry**

| Option | Default | Meaning |
|---|---|---|
| `--alt-freqs` | none | Alternative composition vectors, comma separated; repeat for several. Nucleotide order A,C,G,T |
| `--alt-rates` | none | Alternative exchangeability vectors; GTR takes 6 values (AC, AG, AT, CG, CT, GT); repeat for several |
| `--alt-tstv` | none | Alternative ts/tv ratios for HKY and F84; repeat for several |
| `--alt-freqs-file` | none | File of alternative composition vectors, one per line |
| `--alt-rates-file` | none | File of alternative exchangeability vectors, one per line |

**Site and branch heterogeneity (GFmix)**

| Option | Default | Meaning |
|---|---|---|
| `--gfmix-cseries` | none | C-series profile mixture (C10, C20, C30, C40, C50, C60) for amino-acid site heterogeneity; mutually exclusive with the `--alt-*` options |
| `--gfmix-weights` | built-in | Per-profile weights, comma or space separated, summing to 1 |
| `--gfmix-branch-hetero` | none | File mapping branch names to gamma_G and gamma_F, one `name gamma_G gamma_F` per line |
| `--gfmix-log` | off | Write the per-site profile assignment, weights and branch gammas |

**Per-site rates and composition**

| Option | Default | Meaning |
|---|---|---|
| `-c, --codon-rates` | none | Codon position rates `r1 r2 r3` |
| `-a, --gamma` | none | Gamma shape (alpha). With `-g`, discrete; alone, continuous |
| `-g, --gamma-cats` | none | Number of discrete gamma categories, 2 to 32 |
| `-i, --invariable` | none | Proportion of invariable sites, 0 to 1 |
| `-f, --freqs` | model | State frequencies: 4 for nucleotides, 20 for amino acids, or `e` for equal |
| `-t, --tstv` | equal | Transition/transversion ratio; HKY and F84 only |
| `-r, --rates` | model | Rate matrix values: 6 for GTR, 190 for amino acids |
| `-wr, --write-rates` | off | Write the relative rate of each site to stderr |

**Output**

| Option | Default | Meaning |
|---|---|---|
| `-o, --format` | `f` | `p` PHYLIP, `r` relaxed PHYLIP, `n` NEXUS, `f` FASTA; append `s` for one file per dataset (FASTA is always per-dataset) |
| `-y, --output-prefix` | `output` | Output directory or file prefix |
| `-u, --suffix` | derived | Output file suffix or extension |
| `-wa, --write-ancestral` | off | Write ancestral sequences for internal nodes |
| `--ancestor` | none | Ancestral root sequence; used instead of drawing one |
| `-x, --insert-file` | none | Text file inserted after each dataset |
| `--version` | | Print the version and exit |
| `--manual` | | Print the manual and exit |

## Input formats

**Newick.** Branch lengths, internal labels, and the annotations `[f=N]` and
`[r=N]` selecting registry entries. Quoted names are supported. Internal labels
are stored, so a branch can be named in `--gfmix-branch-hetero`.

```
((taxon_01:0.05,taxon_02:0.05)het0:0.10,(taxon_03:0.05,taxon_04:0.05):0.10);
```

**NEXUS.** BEAST, MrBayes and IQ-TREE tree blocks. The `TRANSLATE` map and the
`[&R]` / `[&U]` root flags are honoured; bracket comments that are not `f=` or
`r=` annotations are stripped.

## Output formats

| Flag | Format | Notes |
|---|---|---|
| `-o p` | PHYLIP strict | Names truncated to the strict field width |
| `-o r` | PHYLIP relaxed | Full names preserved |
| `-o n` | NEXUS | Includes ancestral sequences under `-wa` |
| `-o f` | FASTA | The default |

Appending `s` (`-o ps`) writes one file per dataset instead of one combined file.

## Data provenance

The C-series profile frequencies and weights in `data/cseries/` are from
Si Quang, Gascuel & Lartillot (2008), transcribed from IQ-TREE 2's
`model/modelmixture.cpp` (the blocks labelled "CAT-CXX profile mixture model of
Le, Gascuel & Lartillot (2008)"). The amino-acid column order in every
`.freq.dat` file is `A R N D C Q E G H I L K M F P S T W Y V`.

See [data/cseries/SOURCES.md](data/cseries/SOURCES.md) for the cross-validation
against the GFmix distribution and the one archived discrepancy in the C50 file,
and [docs/REFERENCES.md](docs/REFERENCES.md) for the full bibliography with
verification tiers.

## Reproducibility

A fixed `-z` and a fixed tree reproduce a run byte for byte.
`tests/test_golden.py` locks three output hashes, so any change to the order or
the arithmetic of the random draws fails the suite. The generator is
module-level and global; the package is not thread-safe.

## Tests

```bash
pytest
```

## Verification against reference implementations

The method is stated with each number so it can be repeated. Sources and
locators are in [docs/REFERENCES.md](docs/REFERENCES.md).

- Empirical amino-acid matrices: JTT, WAG, PAM, BLOSUM62 and LG are identical to
  Seq-Gen 1.3.5's `source/aamodels.c`, the 190 exchangeabilities and the 20
  frequencies agreeing entry by entry, in the same order, to
  `max|difference| = 0.000e+00`.
- C-series profiles: C10-C60 agree with IQ-TREE 2's transcription to
  `max|difference| = 0.000e+00`, and with the GFmix distribution to about
  `5e-07`, consistent with the precision at which GFmix stores its tables. The
  GFmix `C50.aafreq.dat` does not match and is documented as an archived
  discrepancy in `data/cseries/SOURCES.md`; the IQ-TREE values are used
  throughout.
- Equilibrium transform: the rescaled profile sums to 1 to
  `|sum - 1| <= 8.882e-16` over 6 C-series sets and 25 `(gamma_G, gamma_F)` pairs,
  and `gamma_G = gamma_F = 1` is the identity to `1.110e-16`.

## Citing this software

See [CITATION.cff](CITATION.cff). The site-and-branch heterogeneity model is
McCarthy et al. (2026) and the C-series profiles are Si Quang et al. (2008);
cite both if you use those features.

## Contributing

Issues and pull requests are welcome. Two constraints on any change:

1. Simulation output is byte-stable for a fixed seed and options within a
   version. If a change moves an output hash, `tests/test_golden.py` fails by
   design. Regenerate the hashes only when the arithmetic change is intended, and
   say so in the change description.
2. The empirical model data are reference data. Changes to `data/cseries/*.dat`
   or to the matrices in `core/aa_models.py` need a source, not a fit.

Run `pytest` before opening a pull request. The manual in `manual.py` is part of
the interface: a new option is not complete until `--manual` documents it.

## Acknowledgements

The simulation engine reimplements the design of Seq-Gen (Rambaut & Grassly
1997) and contains no Seq-Gen source code. The discrete gamma follows Yang (1994)
as implemented in PAML; no PAML code is used, and every probability computation
goes through SciPy's public API. The node-heterogeneous model is inspired by
Foster (2004) and P4, the reference implementation of that idea. It is an
independent implementation written from the paper: no P4 code is used, and P4 is
distributed under the GNU GPL v2 while this package is MIT. The C-series data were
cross-checked against the `gfmixf` distribution by E. Susko.

## License

MIT License, Copyright (c) 2024-2026 seq_sim_v2 authors. See [LICENSE](LICENSE)
for the full text.
