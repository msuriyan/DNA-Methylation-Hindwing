# DNA-Methylation-Hindwing

Nanopore-based multi-omics pipeline for hindwing samples, integrating pairwise
differential gene expression (DEG), differential methylation (DMG), and
differential splicing (DSG) across sample comparisons.

## What DEG / DMG / DSG mean here

| Term | Meaning |
|------|---------|
| DEG | Pairwise differential gene expression between two samples/conditions |
| DMG | Pairwise differential methylation analysis between two samples/conditions (5mC) |
| DSG | Differentially spliced genes from a pairwise comparison |

## Repository structure

```
.
├── Differential methylation analysis/   # scripts: Nanopore mapping -> 5mC/5hmC calling -> DSS methylation testing
├── DEG/                                  # output: per-comparison differential expression results (gene-annotated)
├── DMG/                                  # output: per-comparison differential methylation results (gene-annotated, *_DMR_Annotated.txt)
├── DSG/                                  # output: per-comparison differential splicing results (gene-annotated)
├── DEG_DSG_DMG_Plot_Run.txt              # integration script: combines DEG + DMG + DSG per comparison
├── Methylation run bash_script.txt       # Nanopore mapping / methylation-calling pipeline (see note below)
└── README.md
```

`DEG/`, `DMG/`, and `DSG/` hold the output of each upstream analysis --
result files with gene annotation attached -- which are then read as input by
`DEG_DSG_DMG_Plot_Run.txt` to produce the combined Venn diagrams, scatter
plots, and gene lists per comparison.

## Notes

- `Methylation run bash_script.txt` is the Nanopore mapping + 5mC/5hmC
  calling pipeline. If this is the same pipeline as the scripts inside
  `Differential methylation analysis/`, consider consolidating to one
  location to avoid the two drifting out of sync.
- `DEG_DSG_DMG_Plot_Run.txt` and `Methylation run bash_script.txt` are
  plain-text (`.txt`) files containing R and bash code respectively
