# DSS input format

`Rscript_BsObjM_paired_All.R` reads one file per sample via `read.table(header = TRUE)`.
Each file must be whitespace-delimited with exactly these 4 columns:

| column | meaning |
|--------|---------|
| `chr`  | chromosome / contig name |
| `pos`  | genomic position (1-based) |
| `N`    | total read coverage at this position |
| `X`    | methylated read count at this position (X <= N) |

Example format

chr             pos   N   X
NW_019862748.1  1965  5   5
NW_019862748.1  2371  5   4
NW_019862748.1  5150  8   1
NW_019862748.1  5318  9   2
NW_019862748.1  5494  10  2
NW_019862748.1  5672  8   2
NW_019862748.1  5673  6   3
NW_019862748.1  5703  10  2
NW_019862748.1  5704  6   1
`
