#################### Functional enrichment: DMGs and DMG ∩ DEG overlap ####################
# Modified from Shen and Monteiro 2022. 
#EXAMPLE / TEMPLATE - adjust file names, column names, and comparison labels
# to match your actual DMG (differentially methylated genes) and DEG
# (differentially expressed genes) tables before running.
#
# Runs the same GO/KEGG compareCluster() enrichment as functional_enrichment_polished.R,
# but for two different gene sets:
#   1) DMGs only, per comparison
#   2) The overlap of DMGs and DEGs, per comparison (genes that are both
#      differentially methylated AND differentially expressed)

library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(stringr)

# ================================
# Input files / parameters
# ================================

go2name_file   <- "GO2name.csv"
kegg2name_file <- "KEGG2name.csv"
go2gene_file   <- "go2gene.txt"
kegg2gene_file <- "kegg2gene.txt"

# EXAMPLE input tables - replace with your real files.
# Expected format: wide, one column per comparison, gene IDs down each column
# (same layout as Genelist.csv from the DSG example), e.g.:
#   Wr60      PP50      P15       P50
#   geneA     geneC     geneE     geneG
#   geneB     geneD     geneF     NA
dmg_file <- "DMG_list.csv"   # differentially methylated genes, per comparison
deg_file <- "DEG_list.csv"   # differentially expressed genes, per comparison

comparison_names <- c("Wr60", "PP50", "P15", "P50")  # column order in DMG_list.csv / DEG_list.csv

go_pvalue_cutoff   <- 0.1
kegg_pvalue_cutoff <- 1

outdir <- "."  # change to write outputs elsewhere, e.g. "enrichment_results"

# ================================
# Import GO/KEGG annotations
# ================================

go2name_data   <- read.csv(go2name_file,   header = FALSE, sep = ",")
kegg2name_data <- read.csv(kegg2name_file, header = FALSE, sep = ",")
go2gene_data   <- read.table(go2gene_file,   header = FALSE)
kegg2gene_data <- read.table(kegg2gene_file, header = FALSE)

# ================================
# Import DMG and DEG gene lists (wide format, one column per comparison)
# ================================

dmg_all <- read.csv(dmg_file, header = TRUE, sep = ",")
names(dmg_all) <- comparison_names

deg_all <- read.csv(deg_file, header = TRUE, sep = ",")
names(deg_all) <- comparison_names

# ================================
# Build the DMG ∩ DEG overlap, per comparison
# Each column becomes the set of genes that are BOTH a DMG and a DEG
# in that comparison.
# ================================

overlap_all <- lapply(comparison_names, function(cn) {
  dmg_genes <- na.omit(dmg_all[[cn]])
  deg_genes <- na.omit(deg_all[[cn]])
  intersect(dmg_genes, deg_genes)
})
names(overlap_all) <- comparison_names

cat("DMG-DEG overlap sizes per comparison:\n")
print(sapply(overlap_all, length))

# compareCluster() accepts either a data.frame (as in dmg_all/deg_all) or a
# named list of gene vectors (as we just built for the overlap) - both work
# the same way here.

# ================================
# Helper: run compareCluster + build a labelled dotplot
# (identical to functional_enrichment_polished.R)
# ================================

run_enrichment_dotplot <- function(gene_lists, term2gene, term2name,
                                    pvalue_cutoff, show_category,
                                    label_trunc = 50) {

  result <- compareCluster(
    gene_lists,
    fun = "enricher",
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH"
  )

  clusters_present <- levels(as.factor(as.data.frame(result)$Cluster))

  p <- dotplot(result, showCategory = show_category, includeAll = TRUE, color = "pvalue") +
    scale_color_continuous(low = "orange", high = "navy") +
    ggpubr::rotate_x_text() +
    scale_y_discrete(label = function(x) stringr::str_trunc(x, label_trunc)) +
    scale_x_discrete(label = clusters_present)

  list(result = result, plot = p)
}

# ================================
# 1) DMG enrichment (GO + KEGG)
# ================================

dmg_go <- run_enrichment_dotplot(
  gene_lists = dmg_all, term2gene = go2gene_data, term2name = go2name_data,
  pvalue_cutoff = go_pvalue_cutoff, show_category = 8
)
write.csv(as.data.frame(dmg_go$result), file = file.path(outdir, "GO_DMG.csv"))
ggsave(dmg_go$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "GO_DMG.tiff"))
ggsave(dmg_go$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "GO_DMG.eps"))

dmg_kegg <- run_enrichment_dotplot(
  gene_lists = dmg_all, term2gene = kegg2gene_data, term2name = kegg2name_data,
  pvalue_cutoff = kegg_pvalue_cutoff, show_category = 100
)
write.csv(as.data.frame(dmg_kegg$result), file = file.path(outdir, "KEGG_DMG.csv"))
ggsave(dmg_kegg$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "KEGG_DMG.tiff"))
ggsave(dmg_kegg$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "KEGG_DMG.eps"))

# ================================
# 2) DMG ∩ DEG overlap enrichment (GO + KEGG)
# ================================

overlap_go <- run_enrichment_dotplot(
  gene_lists = overlap_all, term2gene = go2gene_data, term2name = go2name_data,
  pvalue_cutoff = go_pvalue_cutoff, show_category = 8
)
write.csv(as.data.frame(overlap_go$result), file = file.path(outdir, "GO_DMG_DEG_overlap.csv"))
ggsave(overlap_go$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "GO_DMG_DEG_overlap.tiff"))
ggsave(overlap_go$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "GO_DMG_DEG_overlap.eps"))

overlap_kegg <- run_enrichment_dotplot(
  gene_lists = overlap_all, term2gene = kegg2gene_data, term2name = kegg2name_data,
  pvalue_cutoff = kegg_pvalue_cutoff, show_category = 100
)
write.csv(as.data.frame(overlap_kegg$result), file = file.path(outdir, "KEGG_DMG_DEG_overlap.csv"))
ggsave(overlap_kegg$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "KEGG_DMG_DEG_overlap.tiff"))
ggsave(overlap_kegg$plot, width = 6.5, height = 4.5, filename = file.path(outdir, "KEGG_DMG_DEG_overlap.eps"))