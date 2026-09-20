# References

Each entry records where the code uses it and how it was verified.

- **T1** primary source fetched, or metadata fetched from Crossref, or a
  byte-level comparison against the source file, with a locator.
- **T2** carried from this repository or from a reference implementation's own
  bibliography; not independently re-fetched.
- **UNVERIFIED** not obtained. Do not repeat these as fact before checking.

## 1. Substitution models

| Entry | Used by | Tier |
|---|---|---|
| Rambaut, A. and Grassly, N. C. (1997) Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees. *Bioinformatics*, 13(3): 235-238. https://doi.org/10.1093/bioinformatics/13.3.235 | design of the simulation engine (`core/evolve.py`, `tree/tree_parser.py`) | T1 |
| Jukes, T. H. and Cantor, C. R. (1969) Evolution of protein molecules. In Munro, H. N. (ed.) *Mammalian Protein Metabolism*, Academic Press, New York, pp. 21-132. https://doi.org/10.1016/b978-1-4832-3211-9.50009-7 | JC69 (`core/nuc_models.py`) | T1 |
| Hasegawa, M., Kishino, H. and Yano, T. (1985) Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. *Journal of Molecular Evolution*, 22: 160-174 | HKY (`core/nuc_models.py`) | T2 |
| Felsenstein, J. (1984) Distance methods for inferring phylogenies: a justification. *Evolution*, 38(1): 16-24. https://doi.org/10.1111/j.1558-5646.1984.tb00255.x | F84 (`core/nuc_models.py`) | T1 |
| Tavare, S. (1986) Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17: 57-86 | GTR (`core/nuc_models.py`) | T2, see UNVERIFIED |
| Jones, D. T., Taylor, W. R. and Thornton, J. M. (1992) The rapid generation of mutation data matrices from protein sequences. *Computer Applications in the Biosciences*, 8(3): 275-282. https://doi.org/10.1093/bioinformatics/8.3.275 | JTT (`core/aa_models.py`) | T1 |
| Whelan, S. and Goldman, N. (2001) A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Molecular Biology and Evolution*, 18(5): 691-699. https://doi.org/10.1093/oxfordjournals.molbev.a003851 | WAG (`core/aa_models.py`) | T2 |
| Dayhoff, M. O., Schwartz, R. M. and Orcutt, B. C. (1978) A model of evolutionary change in proteins. In Dayhoff, M. O. (ed.) *Atlas of Protein Sequence and Structure*, Vol. 5, Suppl. 3, National Biomedical Research Foundation, Washington DC, pp. 345-352 | PAM (`core/aa_models.py`) | T2 |
| Henikoff, S. and Henikoff, J. G. (1992) Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89: 10915-10919. https://doi.org/10.1073/pnas.89.22.10915 | BLOSUM62 (`core/aa_models.py`) | T1 |
| Kosiol, C. and Goldman, N. (2005) Different versions of the Dayhoff rate matrix. *Molecular Biology and Evolution*, 22(2): 193-199. https://doi.org/10.1093/molbev/msi005 | the DCMUT correction on the bundled JTT and PAM values | T1 |
| Le, S. Q. and Gascuel, O. (2008) An improved general amino acid replacement matrix. *Molecular Biology and Evolution*, 25(7): 1307-1320. https://doi.org/10.1093/molbev/msn067 | LG (`core/aa_models.py`) | T2 |
| Foster, P. G. (2004) Modeling compositional heterogeneity. *Systematic Biology*, 53(3): 485-495. https://doi.org/10.1080/10635150490445779 | the per-node model registry (`core/models.py`) and the `[f=N]` / `[r=N]` annotations (`tree/tree_parser.py`) | T1 |

Four notes on form.

- The 1997 Seq-Gen paper appeared in *Computer Applications in the Biosciences*
  (CABIOS), renamed *Bioinformatics* in 1998. Crossref records the journal as
  *Bioinformatics*; Seq-Gen's own manual still cites the old name. Both refer to
  the same article.
- The Kosiol & Goldman title is often quoted as "Different versions of the
  Dayhoff and JTT models of protein evolution". That is not the published title;
  Crossref records *Different Versions of the Dayhoff Rate Matrix*. Crossref
  dates the article 2004-10-13 (online), while the printed issue is volume 22(2),
  2005. This project cites the printed year.
- Yang (1994) is two different papers and both are cited here. The gamma
  algorithm is Yang, Z. (1994) Maximum likelihood phylogenetic estimation from
  DNA sequences with variable rates over sites: approximate methods. *Journal of
  Molecular Evolution*, 39: 306-314. https://doi.org/10.1007/bf00160154 (T1).
  The substitution-pattern paper is Yang, Z. (1994) Estimating the pattern of
  nucleotide substitution. *Journal of Molecular Evolution*, 39: 105-111
  (T2; Crossref has the volume and DOI https://doi.org/10.1007/bf00178256 but no
  page range).
- Foster's reference implementation of the same idea is P4,
  https://github.com/pgfoster/p4-phylogenetics, distributed under the GNU GPL v2.
  The per-node model registry here is an independent implementation written from
  the paper: it contains no P4 code and shares no identifiers with it, and the
  `(comp_idx, rmat_idx)` keying and the `[f=N]`/`[r=N]` annotation syntax are this
  project's own. P4's licence therefore does not reach this package, which is MIT.
  If you use the node-heterogeneous features, cite Foster (2004).

## 2. Site heterogeneity

- Si Quang, L., Gascuel, O. and Lartillot, N. (2008) Empirical profile mixture
  models for phylogenetic reconstruction. *Bioinformatics*, 24(20): 2317-2323.
  https://doi.org/10.1093/bioinformatics/btn445 (T2; the data themselves were
  verified numerically, see section 5). Used by `core/cseries.py` and
  `core/profile_mixture.py`.

- Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D.,
  von Haeseler, A. and Lanfear, R. (2020) IQ-TREE 2: new models and efficient
  methods for phylogenetic inference in the genomic era. *Molecular Biology and
  Evolution*, 37(5): 1530-1534. https://doi.org/10.1093/molbev/msaa015 (T1).
  This is the transcription source for the C-series files: IQ-TREE 2,
  `model/modelmixture.cpp`, blocks labelled "CAT-CXX profile mixture model of
  Le, Gascuel & Lartillot (2008)".

## 3. Branch heterogeneity

- McCarthy, C. G. P., Susko, E., Harada, R. and Roger, A. J. (2026) Modeling
  Site-and-Branch-Heterogeneity with GFmix. *Systematic Biology*, advance
  article, syag040. https://doi.org/10.1093/sysbio/syag040 (T1). Source of
  `core/profile_mixture.branch_transform`; equations 1.4 and 1.7 give the
  rescaled profile and the constraint that fixes the normalising factor.

- McCarthy, C. G. P., Susko, E. and Roger, A. J. (2025) Modeling
  site-and-branch-heterogeneity with GFmix. *bioRxiv* 2025.08.07.669136, preprint
  version 1, posted 9 August 2025. https://doi.org/10.1101/2025.08.07.669136 (T1).
  Equation numbers quoted in this project follow this version, which the
  implementation was written against. The caption differs between the two
  records and each entry reproduces its own; the published version adds a fourth
  author (Harada), who is not on the preprint record.

- The weighted-average constraint `sum_c w_c mu_c = 1` (Eq. 1.8) is not imposed by
  this code. The source states that it "is not essential and not present in the
  maximum likelihood implementations GFP and GFF", and `gfmixf`, the estimator used
  to check parameter recovery, is one of those ML implementations. A simulator
  that does impose the constraint differs from this one by a scale factor, which
  must not be read as a bias.

- Susko, E. (2026) *gfmixf*, version 1.0, 8 September 2026, Dalhousie
  University. https://www.mathstat.dal.ca/~tsusko/software.html (T2; version
  string read from the driver header, not independently re-fetched).

## 4. Numerical algorithms

| Entry | Used by | Tier |
|---|---|---|
| AS 32: Bhattacharjee, G. P. (1970) The incomplete gamma integral. *Applied Statistics*, 19: 285-287 | gamma quantiles (`core/gamma.py`) | T2 |
| AS 70: Odeh, R. E. and Evans, J. O. (1974) The percentage points of the normal distribution. *Applied Statistics*, 23: 96-97 | gamma quantiles (`core/gamma.py`) | T2 |
| AS 91: Best, D. J. and Roberts, D. E. (1975) The percentage points of the chi-squared distribution. *Applied Statistics*, 24: 385-388 | gamma quantiles (`core/gamma.py`) | T2 |

The discrete gamma follows Yang (1994), section 1. Seq-Gen 1.3.5 implements the
same conditional-mean convention and additionally exposes a median variant
(`source/gamma.c`, `DiscreteGamma` with a `median` argument). The two
implementations were not compared numerically.

## 5. Cross-validation of the bundled data

Measured for this document. Each row is reproducible by comparing the repository
files against the named source with a plain element-wise difference.

| Claim | Measured | Against |
|---|---|---|
| JTT, WAG, PAM, BLOSUM62, LG exchangeabilities and frequencies | `max\|difference\| = 0.000e+00`, same order | Seq-Gen 1.3.5, `source/aamodels.c` |
| C10-C60 profile frequencies | `max\|difference\| = 0.000e+00` | IQ-TREE 2, `model/modelmixture.cpp` |
| C10-C40 and C60 profile frequencies | `max\|difference\|` between `4.8e-07` and `5.0e-07` | GFmix distribution `C*.aafreq.dat`, nine significant digits |
| GFmix `C50.aafreq.dat` | `max\|difference\| = 5.313e-01` | archived discrepancy, see `data/cseries/SOURCES.md` |
| Rescaled profile sums to 1 | `\|sum - 1\| <= 8.882e-16`, 6 sets x 25 `(gamma_G, gamma_F)` pairs | `core/profile_mixture.branch_transform` |
| `gamma_G = gamma_F = 1` is the identity | `max\|difference\| = 1.110e-16` | same |
| `gamma_G` / `gamma_F` shared by every profile on a branch | `max\|difference\| = 0.000e+00`, 10 profiles | `core/profile_mixture.build_profile_models` |

## UNVERIFIED

| Item | Why it is open | How to close it |
|---|---|---|
| Tavare (1986), *Lectures on Mathematics in the Life Sciences* 17: 57-86 | absent from Crossref; the series is not indexed there | a library holding the series, or another implementation's bibliography that states the page range |
| Page range of Yang (1994), *Journal of Molecular Evolution*, 39 | Crossref has the volume and DOI but no pages; 105-111 comes from Seq-Gen's bibliography | the printed article |
| Volume number of the Jukes & Cantor (1969) chapter | Crossref's record has no volume field. The page range 21-132 is confirmed, and the next chapter in the same book starts at page 133 | the printed book |
| Equation numbering in the published version of McCarthy et al. (2026) | the publisher returns HTTP 403 to automated fetches and Europe PMC marks the article subscription-only | institutional access; until then quote equation numbers as "per preprint v1" |
