# Examples

`tree.nwk` is a four-taxon Newick tree. `partitioned.nex` holds two partitions
(60 and 40 sites) with `[f=1]` and `[r=1]` annotations on the edges into H1 and
H2.

```bash
# Newick, nucleotide model, NEXUS output -> out/out.nex
seq-sim -m HKY -l 200 -z 12345 -o n -y out examples/tree.nwk

# C10 profile mixture, amino-acid model -> out_gfmix/out_gfmix_1.fasta
seq-sim -m WAG -l 100 -z 7 --gfmix-cseries C10 -o f -y out_gfmix examples/tree.nwk

# Two partitions summing to 100 sites -> out_nex/out_nex.phy
seq-sim -m GTR -l 100 -z 7 -r 1,2,1,1,2,1 \
    --alt-freqs 0.3,0.2,0.2,0.3 --alt-rates 2,1,1,1,1,2 \
    -o p -y out_nex examples/partitioned.nex
```

Branch lengths are used as given unless `-s` or `-d` rescales the tree. `-wa`
adds ancestral sequences, `-wr` prints the per-site rate table to stderr.
