# seq_sim

Simulates nucleotide and amino-acid sequences along a phylogenetic tree, with
the models and command-line options of Seq-Gen. Partitioned datasets, discrete
gamma, per-branch rate multipliers and C-series profile mixtures are supported.

```bash
pip install .
seq-sim -m HKY -l 1000 -t 4.0 examples/tree.nwk
```

Sequences go to `output/output_1.fasta` unless `-y` gives another prefix.
Requires Python 3.9+ and NumPy 1.21+ / SciPy 1.7+.

## Installation

```bash
pip install .                 # `seq-sim` and `python -m seq_sim`
pip install -e ".[test]"      # development: adds pytest and ruff
```

Tests run from the repository root with `pytest`; `tests/conftest.py` puts
`src/` on `sys.path`, so the checkout directory can be named anything.

## Quick start

`-m` is required and the tree file is positional (stdin when omitted).

```bash
# 1000 bp under HKY, ts/tv = 4.0, FASTA (default)
seq-sim -m HKY -l 1000 -t 4.0 tree.nwk

# 500 amino-acid sites under WAG
seq-sim -m WAG -l 500 tree.nwk

# GTR with four discrete gamma categories, shape 0.5
seq-sim -m GTR -a 0.5 -g 4 tree.nwk

# Tree height normalised to 0.5 substitutions/site (rooted tree)
seq-sim -m HKY -d 0.5 tree.nwk

# C10 profile mixture with branch heterogeneity, 3 replicates, PHYLIP
seq-sim -m LG -l 500 -n 3 --gfmix-cseries C10 \
    --gfmix-branch-hetero gammas.txt -o p -y run1 tree.nwk

# Full manual
seq-sim --manual
```

Runnable inputs are in `examples/`.

## Options

Models: `JC69`, `HKY`, `F84`, `GTR` (`REV`); `JTT`, `WAG`, `PAM`, `BLOSUM`
(BLOSUM62), `LG`, `GENERAL`.

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
| `--normalise-branch-rates` | off | Scale branch rate multipliers to weighted mean 1.0 |
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
| `-o, --format` | f | `p` PHYLIP, `r` relaxed PHYLIP, `n` NEXUS, `f` FASTA; append `s` for separate files per replicate |
| `-y, --output-prefix` | output | Output directory or file prefix |
| `-u, --suffix` | derived | Output file extension |
| `-wa, --write-ancestral` | off | Write ancestral sequences for internal nodes |
| `--ancestor` | none | Fixed ancestral root sequence string |
| `-x, --insert-file` | none | Text file inserted after each dataset |
| `--version` | off | Print version and exit |
| `--manual` | off | Print full manual and exit |

## Input and output

Trees are Newick (branch lengths, quoted names, internal labels, `[f=N]` and
`[r=N]` annotations) or NEXUS `BEGIN TREES` blocks from BEAST, MrBayes and
IQ-TREE, including `TRANSLATE` tables and `[&R]` / `[&U]` flags.

| Flag | Format | Notes |
|---|---|---|
| `-o p` | PHYLIP strict | Names padded/truncated to 10 characters |
| `-o r` | PHYLIP relaxed | Full taxon names preserved |
| `-o n` | NEXUS | Emits `DATA` blocks (including ancestral nodes under `-wa`) |
| `-o f` | FASTA | Default format (always one file per dataset) |

Appending `s` (`-o ps`, `-o rs`, `-o ns`) writes one file per replicate instead
of a single concatenated file.

The C-series data in `src/seq_sim/data/cseries/` come from IQ-TREE 2 and were
compared with the GFmix distribution; see
[SOURCES.md](src/seq_sim/data/cseries/SOURCES.md). Citations and measured bounds
are in [docs/REFERENCES.md](docs/REFERENCES.md).

## Notes

- `-f` sets the distribution the root is drawn from. For length L, state j
  appears L*pi_j times on average (variance L*pi_j*(1-pi_j)), so the realised
  count is not fixed. Use `--ancestor` to pin the root string.
- Under `--gfmix-cseries`, `--ancestor` overwrites the root state after the
  profile for that site has been drawn.
- `--normalise-branch-rates` works per dataset: each dataset is scaled by its
  own length-weighted mean.
- Codon rates restart at every partition: a 1+2 split of one CDS gives position
  classes 1 | 1,2.
- `--gfmix-cseries` loads the C10-C60 tables only; `-f` adds no free-frequency
  class.
- `[f=N]` and `[r=N]` apply to the single edge they are written on; they do not
  inherit to descendants.
- Options a model cannot use are reported on stderr; `-q` suppresses that.

## Not implemented

- Indels and covarion models.
- Mitochondrial and viral amino-acid matrices (`MTREV`, `CPREV45`, `MTART`,
  `HIVB`); the bundled sets are JTT, WAG, PAM, BLOSUM62 and LG.
- `--gfmix-cseries` with nucleotide or codon models; it requires an
  amino-acid model.
- Thread-safe simulation: the generator is module-global, so concurrent runs
  need separate processes.

## Citation and contact

Citing information is in [CITATION.cff](CITATION.cff). Work that uses the
site-and-branch heterogeneity options should also cite McCarthy et al. (2026)
and Si Quang et al. (2008). Issues and pull requests are welcome.

## License

MIT. See [LICENSE](LICENSE). Third-party notices, including Seq-Gen's BSD
3-clause licence, are in [NOTICE](NOTICE).
