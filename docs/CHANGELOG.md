# Changelog

All notable changes to `seq_sim_v2`. Simulation output is deterministic for a
fixed random seed (`-z`) and option set within each minor release.

## 1.4.1

Added GFmix site-and-branch heterogeneity for amino-acid models and support for
named internal nodes in Newick trees.

### Added
- `core/cseries.py` and `data/cseries/`: C10-C60 empirical amino-acid profile
  mixtures and weights (Si Quang et al. 2008), transcribed from IQ-TREE 2
  (`model/modelmixture.cpp`) and cross-checked against the GFmix distribution
  (see `data/cseries/SOURCES.md`).
- `core/profile_mixture.py`: `branch_transform` implementing the GARP/FYMINK
  compositional rescaling (McCarthy et al. 2026), per-site profile sampling,
  and per-branch profile model construction.
- `tests/test_gfmix.py`: unit tests covering C-series loading, branch frequency
  transformation, site profile sampling, and internal node gamma assignment.
- CLI options `--gfmix-cseries`, `--gfmix-weights`, `--gfmix-branch-hetero`,
  and `--gfmix-log`.

### Changed
- `core/evolve.py`: added per-site profile evolution path for root sampling and
  branch mutation.
- `core/aa_models.py`: added `empirical_rates()` helper to share exchangeability
  matrices across profile models.
- `tree/tree_parser.py`: store internal node labels on `node.name` so internal
  branches can be targeted by `--gfmix-branch-hetero`.

### Fixed
- `-` as the tree-file argument now reads the tree from stdin, as `--manual` has
  documented since 1.3.0. Previously it reached `open("-")` and raised
  `FileNotFoundError`.

## 1.3.1

Maintenance and performance improvements.

- Simplified model-registry construction, discrete-gamma fallback helpers, and
  PHYLIP/ancestral sequence writers.
- Avoided allocating per-site branch-length arrays in `_mutate_sequence` when
  substitution rates are uniform across sites, and cached ASCII translation
  tables for sequence output encoding.
- Added pre-flight validation rejecting `-a 0` before gamma quantile evaluation.

## 1.3.0

Packaging and module reorganisation.

- Added PEP 621 `pyproject.toml`, console entry point `seq-sim`, and single-source
  version in `_version.py`.
- Extracted CLI manual text into `manual.py` and rate heterogeneity setup into
  `core/rate_heterogeneity.py`, replacing module-level rate globals with
  explicit state dictionaries.
- Added golden-output SHA-256 regression tests in `tests/test_golden.py`.

## 1.2.2

- Extended `-wa` (ancestral sequence output) to support NEXUS (`-o n`) and
  relaxed PHYLIP (`-o r`) formats.
- Removed unused helper symbols and updated CLI usage examples.

## 1.2.1

- Restricted NEXUS bracket preservation to `[f=N]` and `[r=N]` tokens so
  standard comments such as `[R]` are stripped cleanly.
- Added UTF-8 BOM handling and locale encoding fallback when reading input files.
- Added parser checks rejecting duplicate branch lengths or duplicate `[f=N]`,
  `[r=N]`, and branch-rate annotations on the same node.

## 1.2.0

Bug fixes and model data updates:

- Fixed `-i` (invariable sites) when used without `-a` or `-c`.
- Fixed duplicate subtree traversal under `-wa`.
- Updated the built-in LG exchangeability matrix to match PAML / Seq-Gen 1.3.5
  and documented the DCMUT-corrected JTT and PAM matrices.
- Validated tree rootedness and positive height for `-d` before creating output
  directories.
- Reported zero relative rates for invariable sites under `-wr` when combined
  with `-i` and `-c`.
- Preserved `[f=N]` and `[r=N]` branch annotations during NEXUS preprocessing
  and allowed annotations before or after branch lengths.

## 1.1.0

- Changed amino-acid models to default to their empirical equilibrium
  frequencies (added `-f e` for uniform frequencies).
- Enabled custom `-r` and `--alt-rates` exchangeability vectors on named
  amino-acid models.

## 1.0.0

Initial release supporting nucleotide (JC69, HKY, F84, GTR) and amino-acid
(JTT, WAG, PAM, BLOSUM62, LG, GENERAL) substitution models, gamma/codon/invariable
rate variation, node-heterogeneous models, partitioned alignments, NEXUS input,
and PHYLIP/NEXUS/FASTA output.
