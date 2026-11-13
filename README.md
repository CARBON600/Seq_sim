# seq_sim_v2

Monte Carlo simulation of molecular sequence evolution along phylogenetic trees.

Version: 1.4.1

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [Overview](#overview)
- [Site-and-branch heterogeneity (GFmix)](#site-and-branch-heterogeneity-gfmix)
- [Command-line reference](#command-line-reference)
- [Input formats](#input-formats)
- [Output formats](#output-formats)
- [Data provenance](#data-provenance)
- [Reproducibility and tests](#reproducibility-and-tests)
- [Verification against reference implementations](#verification-against-reference-implementations)
- [Citing this software](#citing-this-software)
- [Contributing](#contributing)
- [Acknowledgements](#acknowledgements)
- [License](#license)

## Installation

Requires Python >= 3.9, NumPy >= 1.21, and SciPy >= 1.7.

```bash
pip install .
# or, for development (includes pytest):
pip install -e ".[test]"
```

This installs the `seq-sim` console script and the `python -m seq_sim_v2` entry
point.

## Quick start

`-m` is required; the tree file is the positional argument (defaults to `stdin`).
Output files are written to `output/` unless `-y` is specified.

```bash
# 1000 bp under HKY with ts/tv = 4.0, FASTA (default)
python -m seq_sim_v2 -m HKY -l 1000 -t 4.0 tree.nwk

# 500 amino-acid sites under WAG
python -m seq_sim_v2 -m WAG -l 500 tree.nwk

# GTR with four discrete gamma categories, shape 0.5
python -m seq_sim_v2 -m GTR -a 0.5 -g 4 tree.nwk

# Normalise total tree height to 0.5 substitutions/site (requires a rooted tree)
python -m seq_sim_v2 -m HKY -d 0.5 tree.nwk

# C10 profile mixture with branch heterogeneity, 3 replicates, PHYLIP
python -m seq_sim_v2 -m LG -l 500 -n 3 --gfmix-cseries C10 \
    --gfmix-branch-hetero gammas.txt -o p -y run1 tree.nwk

# Full command-line manual
python -m seq_sim_v2 --manual
```

## Overview

`seq_sim_v2` evolves nucleotide or amino-acid sequences along a phylogenetic
tree and writes the simulated alignments in standard formats. Every random draw
comes from a single seeded generator, so a fixed `-z` seed reproduces a run
byte-for-byte.

Supported substitution models are `JC69`, `HKY`, `F84`, and `GTR` (`REV`) for
nucleotides, and `JTT`, `WAG`, `PAM`, `BLOSUM` (BLOSUM62), `LG`, and `GENERAL`
for amino acids; the five built-in empirical amino-acid matrices use the exact
exchangeabilities and equilibrium frequencies from Seq-Gen 1.3.5
(`source/aamodels.c`). Per-site rate variation supports discrete gamma (`-a`
with `-g`), continuous gamma (`-a` alone), codon-position rates (`-c`), and
invariable sites (`-i`). Site heterogeneity uses the C10-C60 empirical profile
mixtures (`--gfmix-cseries`), and branch heterogeneity applies per-branch
`gamma_G` / `gamma_F` compositional multipliers (`--gfmix-branch-hetero`).
Branch-specific substitution models can also be assigned via `[f=N]` and `[r=N]`
Newick annotations backed by `--alt-freqs`, `--alt-rates`, and `--alt-tstv`
(Foster 2004).

Input trees can be in Newick format or NEXUS tree blocks (BEAST, MrBayes,
IQ-TREE), with optional partition prefixes and user-supplied ancestral root
sequences. Output formats include strict PHYLIP, relaxed PHYLIP, NEXUS, and
FASTA.

## Site-and-branch heterogeneity (GFmix)

This implements the site-and-branch heterogeneity model of McCarthy et al.
(2026) along two axes:

**Site heterogeneity.** Each site draws one profile from the C-series mixture
(Si Quang et al. 2008) according to the profile weights. With `C` profiles, `C`
rate matrices are built per branch, and site-to-profile assignments are redrawn
for each replicate.

**Branch heterogeneity.** On designated branches, profile frequencies are
rescaled before constructing the rate matrix. Writing `pi^(c)` for profile `c`,
`G` for the GARP amino-acid class (`G, A, R, P`), `F` for FYMINK
(`F, Y, M, I, N, K`), and `O` for the remaining ten amino acids, the transformed
frequencies are:

```
pi^(ce)_j = mu^(ce) * gamma_G * pi^(c)_j     j in G
          = mu^(ce) * gamma_F * pi^(c)_j     j in F
          = mu^(ce) * pi^(c)_j               otherwise
```

where `mu^(ce)` normalises the rescaled profile to sum to 1:

```
mu^(ce) = 1 / ( gamma_G * sum_G pi + gamma_F * sum_F pi + sum_O pi )
```

Because only the equilibrium frequencies are transformed while exchangeabilities
remain unchanged and each rate matrix is normalised to mean rate 1, transformed
branches retain the standard branch-length scale. All profiles on a given branch
share the same `(gamma_G, gamma_F)` pair.

The `--gfmix-branch-hetero` file specifies one branch per line as
`name gamma_G gamma_F` (unlisted branches remain homogeneous):

```
taxon_02  2.0  0.5
het0      0.5  2.0
```

Both tip names and labelled internal nodes in the Newick tree can be referenced.

```bash
# Site heterogeneity only
python -m seq_sim_v2 -m LG -l 1000 --gfmix-cseries C30 -o p -y c30 tree.nwk

# Custom weights plus a per-site assignment log
python -m seq_sim_v2 -m LG -l 1000 --gfmix-cseries C60 \
    --gfmix-weights "$(cat weights60.txt)" --gfmix-log -y c60log tree.nwk

# Site and branch heterogeneity together
python -m seq_sim_v2 -m LG -l 1000 --gfmix-cseries C10 \
    --gfmix-branch-hetero gammas.txt -o p -y both tree.nwk
```

`--gfmix-cseries` requires an amino-acid model and is mutually exclusive with
`--alt-freqs`, `--alt-rates`, and `--alt-tstv`. See
[docs/REFERENCES.md](docs/REFERENCES.md) for literature details.

## Command-line reference

Run `python -m seq_sim_v2 --manual` for the complete manual with examples.

**Model and tree**

| Option | Default | Meaning |
|---|---|---|
| `treefile` | stdin | Newick or NEXUS tree file (positional) |
| `-m, --model` | required | JC69, HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, LG, GENERAL |
| `-l, --length` | 1000 | Sequence length in sites |
| `-n, --num-datasets` | 1 | Number of datasets per tree |
| `-s, --branch-scale` | 0.0 | Branch length scaling factor |
| `-d, --tree-scale` | 0.0 | Normalise total tree height (root to tip) to N substitutions/site; requires a rooted tree; mutually exclusive with `-s` |
| `-br, --branch-rates` | off | Treat `[value]` after a branch length as a per-branch rate multiplier |
| `--normalize-branch-rates` | off | Scale branch rate multipliers to weighted mean 1.0 |
| `-z, --seed` | none | Random number seed |
| `-q, --quiet` | off | Suppress informational output |

**Per-node model registry**

| Option | Default | Meaning |
|---|---|---|
| `--alt-freqs` | none | Alternative composition vectors, comma-separated; repeat for multiple (`A,C,G,T` order for DNA) |
| `--alt-rates` | none | Alternative exchangeability vectors; 6 values for GTR (`AC, AG, AT, CG, CT, GT`), 190 for AA |
| `--alt-tstv` | none | Alternative ts/tv ratios for HKY and F84; repeat for multiple |
| `--alt-freqs-file` | none | File of alternative composition vectors, one per line |
| `--alt-rates-file` | none | File of alternative exchangeability vectors, one per line |

**Site and branch heterogeneity (GFmix)**

| Option | Default | Meaning |
|---|---|---|
| `--gfmix-cseries` | none | C-series profile mixture (`C10`–`C60`) for amino-acid site heterogeneity; mutually exclusive with `--alt-*` |
| `--gfmix-weights` | built-in | Per-profile weights, comma or space separated, summing to 1 |
| `--gfmix-branch-hetero` | none | File mapping branch names to `gamma_G` and `gamma_F` (`name gamma_G gamma_F` per line) |
| `--gfmix-log` | off | Write per-site profile assignments, weights, and branch gammas to `<prefix>.gfmix.log` |

**Per-site rates and composition**

| Option | Default | Meaning |
|---|---|---|
| `-c, --codon-rates` | none | Codon position rates `r1 r2 r3` |
| `-a, --gamma` | none | Gamma shape (alpha); discrete when paired with `-g`, continuous when used alone |
| `-g, --gamma-cats` | none | Number of discrete gamma categories (2 to 32) |
| `-i, --invariable` | none | Proportion of invariable sites (0 to 1) |
| `-f, --freqs` | model | Equilibrium frequencies: 4 for nucleotides, 20 for amino acids, or `e` for equal |
| `-t, --tstv` | equal | Transition/transversion ratio (HKY and F84 only) |
| `-r, --rates` | model | Exchangeability values: 6 for GTR, 190 for amino acids |
| `-wr, --write-rates` | off | Write relative substitution rate of each site to stderr |

**Output**

| Option | Default | Meaning |
|---|---|---|
| `-o, --format` | `f` | `p` PHYLIP, `r` relaxed PHYLIP, `n` NEXUS, `f` FASTA; append `s` for separate files per replicate |
| `-y, --output-prefix` | `output` | Output directory or file prefix |
| `-u, --suffix` | derived | Output file extension |
| `-wa, --write-ancestral` | off | Write ancestral sequences for internal nodes |
| `--ancestor` | none | Fixed ancestral root sequence string |
| `-x, --insert-file` | none | Text file inserted after each dataset |
| `--version` | | Print version and exit |
| `--manual` | | Print full manual and exit |

## Input formats

**Newick.** Supports branch lengths, quoted taxon names, internal node labels,
and `[f=N]` / `[r=N]` model annotations:

```
((taxon_01:0.05,taxon_02:0.05)het0:0.10,(taxon_03:0.05,taxon_04:0.05):0.10);
```

**NEXUS.** Parses `BEGIN TREES` blocks from BEAST, MrBayes, and IQ-TREE,
preserving `TRANSLATE` tables, `[&R]` / `[&U]` rooting flags, and `[f=N]` /
`[r=N]` branch annotations while stripping other metadata comments.

## Output formats

| Flag | Format | Notes |
|---|---|---|
| `-o p` | PHYLIP strict | Names padded/truncated to 10 characters |
| `-o r` | PHYLIP relaxed | Full taxon names preserved |
| `-o n` | NEXUS | Emits `DATA` blocks (including ancestral nodes under `-wa`) |
| `-o f` | FASTA | Default format (always one file per dataset) |

Appending `s` (`-o ps`, `-o rs`, `-o ns`) writes one file per replicate instead
of a single concatenated file.

## Data provenance

The C-series profile frequencies and weights in `data/cseries/` (Si Quang,
Gascuel & Lartillot 2008) are transcribed from IQ-TREE 2's
`model/modelmixture.cpp` in the amino-acid order
`A R N D C Q E G H I L K M F P S T W Y V`. See
[data/cseries/SOURCES.md](data/cseries/SOURCES.md) for cross-validation details
against the `gfmixf` distribution and [docs/REFERENCES.md](docs/REFERENCES.md)
for the bibliography.

## Reproducibility and tests

Run the test suite with `pytest`:

```bash
pytest
```

A fixed `-z` seed and input tree reproduce output files byte-for-byte.
`tests/test_golden.py` checks SHA-256 digests against reference runs to catch
any unintended changes in PRNG draw order or numerical arithmetic.

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

See [CITATION.cff](CITATION.cff). If you use the GFmix site-and-branch
heterogeneity options, please also cite McCarthy et al. (2026) and Si Quang et
al. (2008).

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

MIT License. See [LICENSE](LICENSE) for full text.
