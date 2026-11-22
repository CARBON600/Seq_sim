# References

`Check`: **T1** primary source read, or metadata from Crossref, or a byte-level
comparison with a locator; **T2** carried from a bibliography, not checked
independently; **UNVERIFIED** not obtained, listed at the end.

## Models and algorithms

| Citation | Used by | Check |
|---|---|---|
| Rambaut, A. and Grassly, N. C. (1997) Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees. *Bioinformatics*, 13(3): 235-238. https://doi.org/10.1093/bioinformatics/13.3.235 (the journal was *CABIOS* in 1997; Crossref uses the current name) | `core/evolve.py`, `tree/tree_parser.py` | T1 |
| Jukes, T. H. and Cantor, C. R. (1969) Evolution of protein molecules. In Munro, H. N. (ed.) *Mammalian Protein Metabolism*, Academic Press, New York, pp. 21-132. https://doi.org/10.1016/b978-1-4832-3211-9.50009-7 | JC69 | T1 |
| Hasegawa, M., Kishino, H. and Yano, T. (1985) Dating of the human-ape splitting by a molecular clock of mitochondrial DNA. *Journal of Molecular Evolution*, 22: 160-174 | HKY | T2 |
| Felsenstein, J. (1984) Distance methods for inferring phylogenies: a justification. *Evolution*, 38(1): 16-24. https://doi.org/10.1111/j.1558-5646.1984.tb00255.x | F84 | T1 |
| Tavaré, S. (1986) Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17: 57-86 | GTR | T2, see open items |
| Jones, D. T., Taylor, W. R. and Thornton, J. M. (1992) The rapid generation of mutation data matrices from protein sequences. *Computer Applications in the Biosciences*, 8(3): 275-282. https://doi.org/10.1093/bioinformatics/8.3.275 | JTT | T1 |
| Whelan, S. and Goldman, N. (2001) A general empirical model of protein evolution derived from multiple protein families using a maximum-likelihood approach. *Molecular Biology and Evolution*, 18(5): 691-699. https://doi.org/10.1093/oxfordjournals.molbev.a003851 | WAG | T2 |
| Dayhoff, M. O., Schwartz, R. M. and Orcutt, B. C. (1978) A model of evolutionary change in proteins. In Dayhoff, M. O. (ed.) *Atlas of Protein Sequence and Structure*, Vol. 5, Suppl. 3, National Biomedical Research Foundation, Washington DC, pp. 345-352 | PAM | T2 |
| Henikoff, S. and Henikoff, J. G. (1992) Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89: 10915-10919. https://doi.org/10.1073/pnas.89.22.10915 | BLOSUM62 | T1 |
| Kosiol, C. and Goldman, N. (2005) Different versions of the Dayhoff rate matrix. *Molecular Biology and Evolution*, 22(2): 193-199. https://doi.org/10.1093/molbev/msi005 (online 2004-10-13; the printed issue is cited) | DCMUT correction on the bundled JTT and PAM values | T1 |
| Le, S. Q. and Gascuel, O. (2008) An improved general amino acid replacement matrix. *Molecular Biology and Evolution*, 25(7): 1307-1320. https://doi.org/10.1093/molbev/msn067 | LG | T2 |
| Foster, P. G. (2004) Modeling compositional heterogeneity. *Systematic Biology*, 53(3): 485-495. https://doi.org/10.1080/10635150490445779 | per-node model registry (`core/models.py`), `[f=N]` / `[r=N]` annotations (`tree/tree_parser.py`) | T1 |
| Yang, Z. (1994) Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. *Journal of Molecular Evolution*, 39: 306-314. https://doi.org/10.1007/bf00160154 | discrete gamma (`core/gamma.py`) | T1 |
| Yang, Z. (1994) Estimating the pattern of nucleotide substitution. *Journal of Molecular Evolution*, 39: 105-111. https://doi.org/10.1007/bf00178256 | substitution pattern | T2 |
| Si Quang, L., Gascuel, O. and Lartillot, N. (2008) Empirical profile mixture models for phylogenetic reconstruction. *Bioinformatics*, 24(20): 2317-2323. https://doi.org/10.1093/bioinformatics/btn445 | `core/cseries.py`, `core/profile_mixture.py` | T2 |
| Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A. and Lanfear, R. (2020) IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, 37(5): 1530-1534. https://doi.org/10.1093/molbev/msaa015 | transcription source for the C-series files (`model/modelmixture.cpp`) | T1 |
| McCarthy, C. G. P., Susko, E., Harada, R. and Roger, A. J. (2026) Modeling Site-and-Branch-Heterogeneity with GFmix. *Systematic Biology*, advance article, syag040. https://doi.org/10.1093/sysbio/syag040 | `core/profile_mixture.branch_transform`, equations 1.4 and 1.7; the Eq. 1.8 constraint on `sum_c w_c mu_c` is not applied | T1 |
| McCarthy, C. G. P., Susko, E. and Roger, A. J. (2025) Modeling site-and-branch-heterogeneity with GFmix. *bioRxiv* 2025.08.07.669136, preprint version 1, posted 9 August 2025. https://doi.org/10.1101/2025.08.07.669136 | equation numbers follow this version; the published version adds Harada as fourth author | T1 |
| Susko, E. (2026) *gfmixf*, version 1.0, 8 September 2026, Dalhousie University. https://www.mathstat.dal.ca/~tsusko/software.html | comparison of the C-series data | T2 |
| Bhattacharjee, G. P. (1970) The incomplete gamma integral. *Applied Statistics*, 19: 285-287 (AS 32) | not used; SciPy supplies `gammainc` and `gamma.ppf` | T2 |
| Odeh, R. E. and Evans, J. O. (1974) The percentage points of the normal distribution. *Applied Statistics*, 23: 96-97 (AS 70) | not used; SciPy supplies `gammainc` and `gamma.ppf` | T2 |
| Best, D. J. and Roberts, D. E. (1975) The percentage points of the chi-squared distribution. *Applied Statistics*, 24: 385-388 (AS 91) | not used; SciPy supplies `gammainc` and `gamma.ppf` | T2 |

## Numerical routines

| PAML 4.10.10 | seq_sim |
|---|---|
| `DiscreteGamma`, boundaries from `QuantileGamma`, bin means from `IncompleteGamma` (`src/tools.c:2601`) | `discrete_gamma` (`core/gamma.py:85`), boundaries from `scipy.stats.gamma.ppf`, means from `scipy.special.gammainc` |
| `QuantileGamma`, a macro over `QuantileChi2` (`src/chi2.c:18`) | `scipy.stats.gamma.ppf` |
| `IncompleteGamma`, AS32 series and continued fraction (`src/tools.c:2441`) | `scipy.special.gammainc` |
| `rndgamma`, Marsaglia & Tsang (2000) rejection loop (`src/tools.c:1901`) | `numpy.random.Generator.gamma` (`utils/__init__.py:33`) |

Seq-Gen 1.3.5 uses the same conditional-mean convention (`source/gamma.c`) and,
like PAML, also offers a median variant.

## Data checks

| Quantity | Measured | Against |
|---|---|---|
| JTT, WAG, PAM, BLOSUM62, LG exchangeabilities and frequencies | `max\|difference\| = 0.000e+00`, same order | Seq-Gen 1.3.5, `source/aamodels.c` |
| C10-C60 profile frequencies | `max\|difference\| = 0.000e+00` | IQ-TREE 2, `model/modelmixture.cpp` |
| C10-C40 and C60 profile frequencies | `max\|difference\|` between `4.8e-07` and `5.0e-07` | GFmix `C*.aafreq.dat`, nine significant digits |
| GFmix `C50.aafreq.dat` | `max\|difference\| = 5.313e-01` | archived, see `data/cseries/SOURCES.md` |
| Rescaled profile sums to 1 | `\|sum - 1\| <= 6.7e-16` over the 6 sets and the 25 `(gamma_G, gamma_F)` pairs of `{0.5, 1, 1.5, 2, 2.5}^2` with each stored row normalised; `4.8e-10` on the rows as shipped, which carry 10 significant digits | `core/profile_mixture.branch_transform`; `tests/test_gfmix.py::TestBranchTransform::test_sums_to_one_over_the_grid` |
| `gamma_G = gamma_F = 1` is the identity | `max\|difference\| = 1.110e-16` | same |
| `gamma_G` / `gamma_F` shared by every profile on a branch | `max\|difference\| = 0.000e+00`, 10 profiles | `core/profile_mixture.build_profile_models` |

## Open items

| Item | Why it is open |
|---|---|
| Tavaré (1986), *Lectures on Mathematics in the Life Sciences* 17: 57-86 | absent from Crossref; the series is not indexed there |
| Page range of Yang (1994), *Journal of Molecular Evolution*, 39 | Crossref has the volume and DOI but no pages; 105-111 comes from Seq-Gen's bibliography |
| Volume number of the Jukes & Cantor (1969) chapter | Crossref's record has no volume field; pages 21-132 are confirmed, and the next chapter in the same book starts at 133 |
| Equation numbering in the published version of McCarthy et al. (2026) | the typeset article was not accessible; numbering is quoted per preprint v1 |
