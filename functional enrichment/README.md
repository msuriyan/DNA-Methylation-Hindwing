# Purpose

This R script performs Gene Ontology (GO) and KEGG functional enrichment analyses on selected gene lists, such as the overlap between differentially expressed genes (DEGs) and differentially methylated genes (DMGs), from multiple pairwise comparisons. 
The enrichment analysis follows the workflow described by Tian and Monteiro (2022) for DEG functional enrichment.

# Usage

1. GO and KEGG information for the B. anynana 1.2 genome was obtained from Tian and Monteiro, MBE 2022. 

2. Create files required for the enrichment analysis, which include:

A list of all GO ID vs descriptions (GO2name.csv)

A list of all KEGG ID vs descriptions (KEGG2name.csv)

A list of all GO ID vs genes from your species, as obtained in the GO annotation of the genome (go2gene.txt, for B. anynana v1.2 genome)

A list of all KEGG ID vs genes from your species, as obtained in the KEGG annotation of the genome (kegg2gene.txt, for B. anynana v1.2 genome)

A table of genes of interest (DMG x DEGs) from multiple comparisons. Genelist.csv

3. Perform functional enrichment analysis using functional_enrichment.R







