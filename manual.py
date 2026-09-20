"""Usage manual text for seq_sim_v2."""


def get_manual(version):
    """Return the full command-line manual stamped with *version*."""
    return r"""

================================================================================
  seq_sim v""" + version + r""" — Monte Carlo Simulator of Molecular Sequence Evolution
================================================================================

seq_sim evolves nucleotide or amino acid sequences along a phylogenetic tree
under a specified substitution model. It supports both closed-form models
(JC69, HKY, F84) and eigendecomposition-based models (GTR, JTT, WAG, PAM,
BLOSUM62, LG, GENERAL) with rate heterogeneity (Gamma, codon-position rates,
invariable sites), branch-length scaling, partitioned data, heterogeneous
models per node, and NEXUS/BEAST/MrBayes input.

================================================================================
  I. QUICK START
================================================================================

  (1) Simulate 1000 bp of DNA under HKY with ts/tv=4.0 on a rooted tree:
        $ seq_sim -m HKY -l 1000 -t 4.0 tree.nwk

  (2) Simulate 500 amino-acid sites under the WAG model:
        $ seq_sim -m WAG -l 500 tree.nwk

  (3) Simulate with Gamma rate heterogeneity (4 discrete categories, shape=0.5):
        $ seq_sim -m GTR -a 0.5 -g 4 tree.nwk

  (4) Normalise tree height to 0.5 substitutions/site (requires rooted tree):
        $ seq_sim -m HKY -d 0.5 tree.nwk

  (5) Read tree from stdin, write FASTA to stdout:
        $ cat tree.nwk | seq_sim -m JC69 -l 2000

  (6) Use an ancestral root sequence instead of random generation:
        $ seq_sim -m HKY -l 1000 --ancestor "ACGTACGT..." tree.nwk

  (7) Simulate three independent datasets (replicates) per tree:
        $ seq_sim -m HKY -l 1000 -n 3 tree.nwk


================================================================================
  II. COMMAND-LINE PARAMETER REFERENCE
================================================================================

  POSITIONAL ARGUMENTS
    treefile              Newick or NEXUS tree file. Use "-" or omit for stdin.

  REQUIRED OPTIONS
    -m, --model NAME      Substitution model (REQUIRED). Supported values:
                          NUCLEOTIDE  — JC69, HKY, F84, GTR, REV (=GTR)
                          AMINO ACID  — JTT, WAG, PAM, BLOSUM, LG, GENERAL

  MODEL PARAMETERS
    -t, --tstv FLOAT      Transition/transversion ratio (HKY, F84 only).
                          Must be > 0. Omit for equal exchangeability.

    -f, --freqs FLOAT...  Equilibrium frequencies. 4 values for nucleotides
                          (A,C,G,T order), 20 for amino acids.
                          Omit for default frequencies (equal for nucleotides,
                          model-specific empirical frequencies for amino acids).
                          Use 'e' for equal frequencies.

    -r, --rates FLOAT...  Exchangeability (relative rate) parameters.
                          6 values for GTR (AC,AG,AT,CG,CT,GT order) or
                          190 values for any AA model (overrides the built-in
                          rates; GENERAL requires them).
                          Omit for defaults (all 1.0 for nuc).

    --alt-freqs VEC       Alternative composition vector (comma-separated:
                          "0.3,0.2,0.3,0.2"). Repeat for multiple vectors.
                          Used with [f=N] annotations for per-node models.

    --alt-rates VEC       Alternative exchangeability vector (comma-
                          separated: "1,2,1,1,2,1"). Repeat for multiple.

    --alt-tstv FLOAT      Alternative ts/tv ratio. Repeat for multiple.

    --alt-freqs-file PATH File of alternative frequency vectors (one per line).

    --alt-rates-file PATH File of alternative exchangeability vectors.

  GFmix SITE-AND-BRANCH HETEROGENEITY  (amino-acid models only)
    --gfmix-cseries CXX  C-series profile mixture (C10, C20, C30, C40,
                          C50 or C60) for amino-acid SITE heterogeneity
                          (Si Quang, Gascuel & Lartillot 2008). Each site
                          draws one of the CXX empirical frequency profiles
                          with probability equal to its mixture weight.
                          Requires an AA model (-m JTT/WAG/PAM/BLOSUM/LG/
                          GENERAL). Mutually exclusive with --alt-freqs/
                          --alt-rates/--alt-tstv.

    --gfmix-weights W     Per-profile weights (CXX values summing to 1,
                          comma or space separated). Default: the C-series
                          built-in weights.

    --gfmix-branch-hetero PATH
                          File mapping tip names to branch gamma_G and
                          gamma_F (one 'name gamma_G gamma_F' per line).
                          On those branches the profile frequencies are
                          rescaled per McCarthy et al. (2026) Eq. 1.4:
                          GARP amino acids by gamma_G, FYMINK by gamma_F,
                          the remaining ten by 1, then renormalised so the
                          profile sums to 1. GARP = G,A,R,P (GC-rich
                          codons); FYMINK = F,Y,M,I,N,K (AT-rich codons).

    --gfmix-log           Write a log of the per-site profile assignment,
                          profile weights and branch gammas to
                          <prefix>.gfmix.log.

  SIMULATION SIZE
    -l, --length INT      Sequence length in sites.  Default: 1000.
    -n, --num-datasets N  Number of independent replicates per tree.  Default: 1.

  RATE HETEROGENEITY (mutually exclusive groups: -c | -a/-g | none)
    -a, --gamma SHAPE     Gamma distribution shape parameter (alpha).
                          Typically 0.2–2.0. Enables continuous gamma
                          rate variation across sites.

    -g, --gamma-cats N    Number of discrete gamma categories (2–32).
                          Must be used WITH -a.  Common choices: 4, 8, 16.

    -c, --codon-rates R1 R2 R3
                          Codon-position rates (three values; nucleotide
                          models only). Automatically normalised so mean = 1.
                          Example: -c 1 0.5 2

    -i, --invariable PROP Proportion of invariable sites, 0 <= p < 1.
                          Sites with rate = 0. Can combine with -a/-g.

  BRANCH-LENGTH SCALING (mutually exclusive: -s | -d)
    -s, --branch-scale X  Multiply all branch lengths by X.  Default: 0 (off).

    -d, --tree-scale H    Scale tree so root-to-tip height = H subs/site.
                          Requires a ROOTED tree. Typical values:
                          0.1 (weak signal), 0.5 (moderate), 1.0 (saturated).
                          Tip: use this to control overall divergence level.

    -br, --branch-rates   Interpret [value] annotations in the Newick tree as
                          per-branch rate multipliers:  node.param = value.

    --normalize-branch-rates
                          Scale branch rate multipliers so their weighted
                          mean (by branch length) = 1.0.  Only makes sense
                          with -br.  Implicitly enables -br.

  OUTPUT
    -o, --format FMT      Output format:
                            p,phylip  — strict PHYLIP (10-char padded names)
                            r,relaxed — Relaxed PHYLIP (full names)
                            n,nexus   — NEXUS DATA block
                            f,fasta   — FASTA (default)
                          Append 's' for separate files per replicate:
                            ps, rs, ns, fs  (fs is always separate).

    -y, --output-prefix DIR
                          Output file/directory prefix.  Default: "output".

    -u, --suffix EXT      Output file suffix.  Defaults: .phy, .nex, .fasta

    -wa, --write-ancestral
                          Write internal (ancestral) node sequences.
                          Outputs in the same format as tips.

    -wr, --write-rates    Write per-site relative substitution rates.

    --ancestor SEQ        Use the given sequence string as the ancestral root
                          rather than generating a random root sequence.
    -x, --insert-file PATH
                          Text file to insert after each dataset in the output.

    -z, --seed INT        Random-number seed for reproducibility. Omit for
                          time-based seed (printed to stderr).

    -q, --quiet           Suppress all informational output (progress, seed,
                          model summary, timing).

  MISCELLANEOUS
    --version             Print version and exit.
    --manual              Print this comprehensive manual and exit.
    -h, --help            Print short argument summary and exit.


================================================================================
  III. SUBSTITUTION MODELS
================================================================================

  NUCLEOTIDE MODELS (4 states: A, C, G, T)
  ─────────────────────────────────────────
  JC69    Jukes-Cantor 1969. Equal base frequencies, equal substitution rates.
          No extra parameters needed.

  HKY     Hasegawa-Kishino-Yano 1985. Closed-form transition probabilities.
          Uses equilibrium frequencies and a transition/transversion ratio.
          Parameters: -f (freqs), -t (ts/tv).  Omit -t for equal rates.

  F84     Felsenstein 1984. Similar to HKY but with a different parametrisation
          of the transition/transversion bias.
          Parameters: -f (freqs), -t (ts/tv).

  GTR/REV General Time-Reversible. Symmetric eigendecomposition of the 4x4
          instantaneous rate matrix. The most flexible nucleotide model.
          Parameters: -f (freqs), -r (6 exchangeability rates).

  AMINO ACID MODELS (20 states: ARNDCQEGHILKMFPSTWYV)
  ─────────────────────────────────────────────────────
  JTT     Jones-Taylor-Thornton 1992 (DCMUT-corrected values, Kosiol &
          Goldman 2005). Empirical model derived from large protein
          databases. Good general-purpose model.

  WAG     Whelan & Goldman 2001. Empirically estimated from globular proteins.

  PAM     Dayhoff PAM (Point Accepted Mutation) matrix (DCMUT-corrected values,
          Kosiol & Goldman 2005). Historically important; note that some
          exchangeabilities are zero due to limited data.

  BLOSUM  BLOSUM62-based model. Widely used in protein alignment and phylogeny.

  LG      Le & Gascuel 2008. Modern empirical model, recommended for most
          protein analyses.

  GENERAL User-supplied empirical matrix. Requires:
          -f (20 frequencies) and -r (190 exchangeability rates).

  For any AA model, you may:
    - Omit -f to use the model's default (empirical) frequencies.
    - Pass -f with 20 values to override frequencies.
    - Pass -f e for equal frequencies.
    - Omit -r (the model's preloaded rates are used).
    - Pass -r with 190 values to override the rates (GENERAL requires -r).


================================================================================
  IV. RATE HETEROGENEITY
================================================================================

  Rate heterogeneity allows different sites to evolve at different speeds.

  GAMMA DISTRIBUTION (+G)
    The relative substitution rate at each site is drawn from a Gamma
    distribution with shape = alpha and mean = 1.0.

    Continuous gamma (-a SHAPE):
      Each site receives an independent Gamma variate.
      More realistic but slower (one random draw per site per replicate).

    Discrete gamma (-a SHAPE -g CATEGORIES):
      The Gamma distribution is approximated by N equiprobable categories.
      The mean of the N rate values is exactly 1.0.
      Faster because P(t) matrices are precomputed per category.
      Typical values: -g 4 (standard), -g 8 (more precise).

    Shape parameter guidelines:
      alpha < 0.5  — strong rate variation (many slow, few fast sites)
      alpha = 1.0  — moderate variation
      alpha > 2.0  — weak variation (approaching equal rates)

  CODON-POSITION RATES (-c)
    Models three rate classes corresponding to codon positions 1, 2, 3.
    Example: -c 1 0.5 2 means position 1 = baseline, position 2 slower,
    position 3 twice as fast.  Values are auto-normalised to mean = 1.0.

    Mutually exclusive with Gamma (-a/-g).

  INVARIABLE SITES (+I)
    A fraction of sites never change (rate = 0). The remaining sites evolve
    at rates scaled by 1/(1-p_inv) to preserve the mean rate.
    Can be combined with Gamma models.
    Example: -a 0.5 -g 4 -i 0.2  → 20% invariable + Gamma(0.5, 4 cat).


================================================================================
  V. BRANCH-LENGTH SCALING
================================================================================

  BRANCH SCALE (-s)
    Multiply every branch length by a constant factor.
    -s 2.0  → branches twice as long (more divergence)
    -s 0.1  → branches one-tenth as long (very conserved)
    Use with unrooted OR rooted trees.

  TREE SCALE (-d)
    Normalise the whole tree so the root-to-tip distance (tree height)
    equals exactly H substitutions per site. Requires a rooted tree.
    All branches are scaled proportionally.

    Recommended values:
      -d 0.1  — weak signal, recent divergence
      -d 0.5  — moderate (good for most phylogenetics)
      -d 1.0  — strong signal
      -d 2.0  — saturated (many multiple hits)

    WARNING: if tree.total_length > 2.0 a saturation warning is printed.

  PER-BRANCH RATE MULTIPLIERS (-br)
    Newick trees with [value] annotations after branch lengths:
      (A:0.1[1.5],B:0.2[0.8])
    The [1.5] and [0.8] values become per-branch rate multipliers applied
    as:  effective_length = branch_length * scale * value.

    Use --normalize-branch-rates to globally rescale these multipliers
    so their weighted mean (by branch length) is exactly 1.0.


================================================================================
  VI. OUTPUT FORMATS
================================================================================

  PHYLIP (p)         Strict interleaved format. Taxon names truncated to 10
                      characters. Compatible with RAxML, PhyML, PAML.

  Relaxed PHYLIP (r) Like PHYLIP but full taxon names (space-separated).
                      Compatible with IQ-TREE, RAxML-NG.

  NEXUS (n)          NEXUS DATA block with Dimensions, Format, and Matrix
                      sections. Compatible with MrBayes, PAUP*, BEAST.

  FASTA (f)          Standard FASTA format, 72-char wrapped sequences.
                      Default format. Always writes separate files per
                      replicate when using multiple replicates.

  Separate files (s) Append 's' to format: -o ps, -o rs, -o ns, -o fs.
                      Creates one output file per replicate/dataset:
                        output_01.phy, output_02.phy, ...

  ANCESTRAL SEQUENCES (-wa)
    Writes sequences for ALL nodes (tips + internal) with auto-numbered
    labels. Internal nodes are numbered starting from (num_tips + 1).

  SITE RATES (-wr)
    Prints a table of per-site relative substitution rates to stderr.
    Site 1  0.85321
    Site 2  1.45210
    ...
    Useful for debugging or downstream analyses of rate heterogeneity.


================================================================================
  VII. NODE-HETEROGENEOUS MODELS (Advanced)
================================================================================

  seq_sim supports different substitution models on different branches
  via Newick annotations [f=N] and [r=N] placed after branch lengths:

    Composition index [f=N]:  selects the N-th frequency vector.
    Rate-matrix index [r=N]:  selects the N-th rate matrix / ts/tv set.

  Example — a tree where the outgroup uses different base frequencies:
    ((A:0.1[f=0],B:0.2[f=0]):0.05[f=0],C:0.3[f=1]):0.0;

  The default frequency vector (from -f) is index 0. Alternative vectors
  from --alt-freqs are indices 1, 2, ... The model registry maps each
  (comp_idx, rmat_idx) pair to a separate model instance.

  This feature generalises "mixture models" and "partitioned models" by
  allowing every branch to have its own substitution process.


================================================================================
  VIII. NEXUS INPUT (BEAST, MrBayes, IQ-TREE)
================================================================================

  seq_sim automatically detects NEXUS format and preprocesses:

    → Strips all [...] annotations (handles nesting).
    → Extracts the TRANSLATE table (numeric taxon IDs → names).
    → Extracts clean Newick trees from BEGIN TREES blocks.
    → Detects [&R] (rooted) and [&U] (unrooted) flags per tree.
    → Parses TREE and UTREE commands.

  Supported NEXUS sources:
    BEAST TreeAnnotator (MCC trees with [&height], [&posterior], etc.)
    MrBayes .con.tre (consensus trees)
    IQ-TREE .treefile, RAxML bipartitions (when wrapped in NEXUS)

  Case-insensitive keywords. Multiple TREES blocks are merged.


================================================================================
  IX. COMMON WORKFLOWS
================================================================================

  (a) Maximum-likelihood phylogeny simulation (standard):
        seq_sim -m GTR -a 0.5 -g 4 -l 1000 tree.nwk

  (b) Protein evolution with WAG+Gamma:
        seq_sim -m WAG -a 1.0 -g 8 -l 500 -o f protein_tree.nwk

  (c) NEXUS output for MrBayes analysis:
        seq_sim -m HKY -l 2000 -t 4.0 -o n tree.nwk

  (d) Generate multiple replicates for parametric bootstrap:
        seq_sim -m GTR -l 1000 -n 100 -o fs -z 42 tree.nwk

  (e) Set tree height to exactly 0.3 subs/site (rooted tree):
        seq_sim -m HKY -d 0.3 -l 500 rooted_tree.nwk

  (f) Simulate from BEAST posterior tree distribution:
        seq_sim -m GTR -l 1000 -a 1.0 -g 4 beast.trees

  (g) Use a fixed ancestral sequence (e.g., for hypothesis testing):
        seq_sim -m JC69 -l 100 --ancestor "$(cat ancestral.txt)" tree.nwk

  (h) Codon-model simulation with strong 3rd-position bias:
        seq_sim -m HKY -c 1 0.3 2.5 -l 999 tree.nwk
        (Note: length should be divisible by 3 for codon data.)

  (i) Heterogeneous model: two different GTR matrices on two clades:
        seq_sim -m GTR -r 1,2,1,1,2,1 --alt-rates "3,1,3,1,3,1" tree.nwk
        (Use with [r=0] and [r=1] annotations in the Newick string.)


================================================================================
  X. NOTES AND TROUBLESHOOTING
================================================================================

  1. ROOTED vs. UNROOTED TREES
     -d (tree-scale) requires a rooted tree. A rooted Newick tree has
     exactly two children at the root:  ((A,B),C);
     An unrooted tree has three:        ((A,B),C,D);
     Unrooted trees cannot be height-normalised. Use -s instead.

  2. SATURATION
     When tree height (total root-to-tip length) exceeds 2.0, sequences
     may be too saturated for phylogenetic signal. A warning is printed.
     Use -d to reduce the effective height.

  3. REPRODUCIBILITY
     Always set -z with a fixed seed for reproducible results:
       seq_sim -m HKY -l 1000 -z 42 tree.nwk
     Without -z, a time-based seed is used and printed to stderr.

  4. CODON RATE NORMALISATION
     Codon rates are automatically normalised such that the weighted mean
     over positions 1,2,3 equals 1.0:
       normalised_rate[i] = raw[i] * 3 / sum(raw)

  5. PARTITION MARKERS
     Trees may include [LENGTH] or [LENGTH,RATE] markers at the start.
     These define partition lengths and relative evolutionary rates for
     partitioned analyses (e.g., gene 1 = 300 bp, gene 2 = 700 bp).

  6. INVARIABLE + GAMMA
     When combining -i and -a/-g, the mean of the variable-site rates is
     scaled by 1/(1-p_inv) so the overall expected rate remains 1.0.

  7. LARGE FILES / PERFORMANCE
     For long sequences (>50kbp), large trees (>1000 tips), or many
     replicates, use -q to suppress progress output. The discrete gamma
     method (-g N) is significantly faster than continuous gamma.

  8. ANCESTRAL SEQUENCES
     The --ancestor string length must match -l. Characters must be valid
     state symbols for the chosen model (e.g., A,C,G,T for nucleotides;
     ARNDCQEGHILKMFPSTWYV for amino acids).

  9. NEXUS TREE NAMING
     Multiple TREES blocks in a single NEXUS file are supported and
     merged. TRANSLATE tables are also merged across blocks.

  10. MULTIPLE DATASETS / REPLICATES
      -n N generates N independent replicate datasets from each tree.
      Each replicate gets its own random root sequence and rate assignments.
      Use -o fs to write each replicate to a separate FASTA file.

================================================================================
  seq_sim v""" + version + """
  Monte Carlo simulator of molecular sequence evolution (Seq-Gen compatible)
================================================================================
"""
