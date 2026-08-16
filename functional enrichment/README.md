# Purpose
This R script is used to perform functional enrichment analysis (GO and KEGG) on selected lists of genes (DEGs or DSGs) from multiple comparisons.

# Usage

1. Perform functional annotation for the genome (if GO and/or KEGG annotations are not available) of your species, as described in the original publication. 

2. Create files required for the enrichment analysis, which include:

A list of all GO ID vs descriptions (GO2name.csv)

A list of all KEGG ID vs descriptions (KEGG2name.csv)

A list of all GO ID vs genes from your species, as obtained in the GO annotation of the genome (go2gene.txt, for B. anynana v1.2 genome)

A list of all KEGG ID vs genes from your species, as obtained in the KEGG annotation of the genome (kegg2gene.txt, for B. anynana v1.2 genome)

A table of genes of interest (such as DEGs or DSGs) from multiple comparisons. Here we provided shortlisted DSGs (FDR<0.05, |changes of PSI|>0.1) between seasonal forms from all stages (mylist.csv).

3. Perform functional enrichment analysis using functional_enrichment.R







