# Changelog

All notable changes to seq_sim_v2. Versions follow semantic versioning; the
numerical output is byte-stable within a version for a fixed seed and options.

## 1.4.1

Site and branch heterogeneity following the GFmix model, plus named internal
nodes. The baseline for this entry is a 1.3.1 checkout: the 1.4.0 and 1.4.1
releases cannot be separated from the artefacts available here, so the changes
are reported as one 1.3.1 to 1.4.1 delta.

Simulation output for a fixed seed and options is unchanged for every model and
option that already existed. The suite passes at 58 tests, of which 11 are new.

Added:

- `core/cseries.py` and `data/cseries/`: the C10-C60 empirical profile mixtures
  and their weights, transcribed from IQ-TREE 2's `model/modelmixture.cpp` and
  cross-checked against the GFmix distribution. See `data/cseries/SOURCES.md`.
- `core/profile_mixture.py`: `branch_transform`, which applies the GFmix branch
  transform to one profile frequency vector, plus the per-site profile assignment
  and the per-node model build.
- `tests/test_gfmix.py` (11 tests): an independent reimplementation of the branch
  transform is compared against the code, the sum-to-one property is checked,
  `gamma = (1, 1)` is checked to be the identity, and the per-site profile counts
  are checked against the weights.
- Four options, all documented in `--manual`: `--gfmix-cseries`,
  `--gfmix-weights`, `--gfmix-branch-hetero`, `--gfmix-log`.

Changed:

- `core/evolve.py` (+104/-17): a profile-mixture path, so the root draw and the
  branch mutation can each use a different model per site.
- `core/aa_models.py` (+8): `empirical_rates()`, which returns the 190
  exchangeabilities of a named empirical model, so a profile mixture can reuse one
  model's exchangeabilities across many profiles.
- `tree/tree_parser.py` (+16/-5): internal node labels are stored rather than
  skipped, so a branch can be addressed by name.
- `seq_sim.py` (+166/-2) and `manual.py` (+28): the new options, the branch-file
  parser, the gamma application step, and the GFmix log writer.

## 1.3.1

Code-quality pass. Simulation output is unchanged: the golden hashes in
`tests/test_golden.py` were regenerated from real output and now pass.

- Style: dropped `from __future__` imports and the type-annotation scaffolding
  on internal helpers; removed docstrings that only restated a function name.
- Structure: the four near-identical registry branches collapsed into one loop;
  `discrete_gamma`'s three nested closures became candidate functions; the
  ancestral writers all use one node walk; PHYLIP strict/relaxed share a writer.
- Performance: `_mutate_sequence` no longer allocates a per-site length array
  when rates are uniform; `set_vector` skips the array checks the caller already
  guarantees; sequence encoding goes through a cached byte translation table.
- Fixes: `-a 0` is rejected with a message instead of a ZeroDivisionError;
  `pyproject.toml` no longer points its Homepage at an unrelated project; the
  `-i`/`-d`/`--ancestor` paths keep their pre-flight validation.

## 1.3.0

Internal restructuring for maintainability and packaging. Simulation output is
byte-identical to 1.2.2 for the same seed and options.

- Packaging: added `pyproject.toml` (PEP 621, NumPy/SciPy dependencies, console
  script `seq-sim`) and an MIT `LICENSE`; the version now has a single source in
  `_version.py`.
- Structure: extracted the ~440-line manual into `manual.py`; moved the rate
  model into `core/rate_heterogeneity.py`; removed the `rate_info_global`
  module-level global in favour of explicit parameter passing.
- Typing: added `from __future__ import annotations` and type hints throughout;
  introduced a `SubstitutionModel` Protocol and a `RateInfo` TypedDict.
- Comments/docstrings: English only, concise, key notes retained.
- Tests: added a golden-output regression test (47 tests total).
- Docs: translated the design and changelog to English under `docs/`.

## 1.2.2

- Unified module usage string to `python -m seq_sim_v2`.
- `-wa` now supports NEXUS output (`-o n`); `-o r` writes full (relaxed) names.
- Removed dead code: `_GTR_RATE_LABELS`, `_inv_sqrt_pi`, `count_trees`.
- Full re-verification of all documented line references.

## 1.2.1

- NEXUS `[f=N]`/`[r=N]` detection made strict (`f=`/`r=` only), so ordinary
  comments such as `[R]` are stripped again.
- Input reading: UTF-8 with a locale fallback and warning for non-UTF-8 files.
- Parser rejects duplicate branch lengths and duplicate `f`/`r`/rate
  annotations.

## 1.2.0

Correctness fixes from the audit:

- `-i` alone now actually assigns invariable sites.
- `-wa` no longer writes each subtree twice.
- Removed the ineffective `-p/--max-partitions`.
- Replaced the LG matrix with the authoritative values (PAML/Seq-Gen).
  Documented that JTT and PAM are the DCMUT-corrected matrices.
- `-d` validates rootedness and positive tree height before creating output.
- `-wr` reports zero rates for invariable sites under `-i -c`.
- NEXUS preprocessing preserves `[f=N]`/`[r=N]`.
- Diagnostic error for missing model-registry entries; annotations accepted
  before or after branch lengths; `--ancestor` length/character validation;
  UTF-8 I/O; unified comma/space separators.

## 1.1.0

- Amino-acid models now default to their empirical frequencies (previously
  uniform). Added `-f e` for equal frequencies; named AA models honour
  `-r`/`--alt-rates`; model-aware length validation.

## 1.0.0

Initial release: JC69/HKY/F84/GTR and JTT/WAG/PAM/BLOSUM/LG/GENERAL models,
gamma/codon/invariable rate heterogeneity, node-heterogeneous models,
partitioned data, NEXUS input, and PHYLIP/NEXUS/FASTA output.
