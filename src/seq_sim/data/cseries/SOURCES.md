# C-series data provenance

`C10..C60.freq.dat` and `C10..C60.weight.dat` hold the C-series empirical profile
mixtures of:

    Si Quang, L., Gascuel, O., & Lartillot, N. (2008). Empirical profile mixture
    models for phylogenetic reconstruction. Bioinformatics, 24(20), 2317-2323.
    DOI 10.1093/bioinformatics/btn445.

**Transcription source:** IQ-TREE 2, `model/modelmixture.cpp`, the blocks
attributed to "CAT-CXX profile mixture model of Le, Gascuel & Lartillot (2008)".

**Cross-validation.** Frequencies were checked against two external sources.

1. Against IQ-TREE 2's `model/modelmixture.cpp`, the declared transcription
   source: C10-C60 agree to `max|difference| = 0.000e+00`, that is, entry by
   entry in the same amino-acid order.
2. Against the GFmix source package (`gfmixf.tar.gz`, E. Susko,
   https://www.mathstat.dal.ca/~tsusko/software.html): C10, C20, C30, C40 and C60
   agree to `max|difference|` between `4.8e-07` and `5.0e-07` in absolute
   frequency, consistent with the precision at which GFmix stores its tables.

**Archived discrepancy.** The GFmix `C50.aafreq.dat` does **not** match the
Si Quang 2008 C50: `max|difference| = 5.313e-01` in absolute frequency, or a
maximum relative difference of 1789%. It is recorded here as an archived
discrepancy in the GFmix distribution, and the IQ-TREE C50 is used throughout.

**Amino-acid order** (columns of every `.freq.dat` file):

    A R N D C Q E G H I L K M F P S T W Y V

**Rows and weights.** Each `.freq.dat` file holds one profile per row, in the
order the weights in the matching `.weight.dat` file apply. Every row sums to 1
and the weights sum to 1, both to within 1e-6; `load_cseries` enforces both.

**A note on scale.** The GFmix *estimator* used with these profiles is a maximum
likelihood implementation and does not impose the weighted-average constraint
`sum_c w_c mu_c = 1` on the branch multipliers. Some GFmix-based simulators do
impose it, which makes their multipliers differ from the absolute multipliers this
package takes by a constant scale factor. See the weighted-average constraint
note on the McCarthy et al. (2026) row in `docs/REFERENCES.md`.
