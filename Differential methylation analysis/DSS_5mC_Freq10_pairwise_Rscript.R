# =============================================================================
# DSS_5mC_Freq10_pairwise_Rscript.R
#
# Purpose: run DSS's differential methylation testing for 5mC across 10 specific
# pairwise group comparisons (defined below), using the combined BSseq object
# built by Rscript_BsObjM_paired_All.R.
#
# For each comparison this script:
#   1. subsets the full BSseq object down to just the two groups being compared
#   2. runs DMLtest() (per-site differential methylation test, with smoothing)
#   3. calls significant DMLs (differentially methylated loci) with callDML()
#      and separately filters to a stricter fdr < 0.05 & |diff| >= 0.1 set
#   4. calls DMRs (differentially methylated regions) with callDMR(), which
#      groups nearby significant sites into regions
#
# Outputs per comparison (g1_vs_g2):
#   <g1>_vs_<g2>_5mC_DMLTest_full.csv    -- full per-site test statistics
#   <g1>_vs_<g2>5mC_DML_all.txt          -- all called DMLs (delta=0.05, p.threshold=1)
#   <g1>_vs_<g2>_5mC_DML_filtered.txt    -- DMLs filtered to fdr<0.05 & |diff|>=0.1
#   <g1>_vs_<g2>_5mC_DMR.csv             -- called DMRs
# =============================================================================

library(DSS)
library(bsseq)
library(parallel)

## 5mC
setwd("/hpctmp/snm_02/Methylome/Freq10/DSS")

BSobj_5mC <- readRDS("BSseq_5mC_Freq10_all_samples.rds")

# Define sample groups
# -- maps each condition/group label to its two replicate sample names
groups <- list(
  D15 = c("D15_3", "D15_M"),
  D50 = c("D50_1", "D50_M"),
  D60 = c("D60_3", "D60_M"),
  DP  = c("DP_1", "DP_M"),
  P15 = c("P15_1", "P15_M"),
  P50 = c("P50_1", "P50_M"),
  W60 = c("W60_2", "W60_M"),
  WP  = c("WP_3", "WP_M")
)

# Specify the comparisons to run
# -- each pair is one group-vs-group DSS test (10 comparisons total)
comparisons <- list(
  c("W60", "WP"),
  c("WP", "P15"),
  c("P15", "P50"),
  c("D60", "DP"),
  c("DP", "D15"),
  c("D15", "D50"),
  c("W60", "D60"),
  c("WP", "DP"),
  c("P15", "D15"),
  c("P50", "D50")
)



for (pair in comparisons) {
  g1 <- pair[1]
  g2 <- pair[2]

  samples_g1 <- groups[[g1]]
  samples_g2 <- groups[[g2]]
  selected_samples <- c(samples_g1, samples_g2)

  # Subset BSseq object
  # -- pulls out just the 4 samples (2 per group) relevant to this comparison
  BS_sub <- BSobj_5mC[, sampleNames(BSobj_5mC) %in% selected_samples]

  message("Running DMLtest for ", g1, " vs ", g2)

  # Run DML test
  # -- DSS's core per-site differential methylation test between the two groups,
  # with smoothing enabled to borrow information from neighboring CpG sites
  dml_test <- DMLtest(BS_sub,
                      group1 = samples_g1,
                      group2 = samples_g2,
                      smoothing = TRUE)

  # Save full DML results
  # -- every tested site, before any significance filtering
  dml_test_file <- paste0(g1, "_vs_", g2, "_5mC_DMLTest_full.csv")
  write.csv(dml_test, file = dml_test_file, row.names = FALSE, quote = FALSE)

  # Filter by FDR < 0.05
  # -- callDML() here uses delta=0.05, p.threshold=1, i.e. returns all called DMLs
  # ranked by significance (not yet hard-filtered -- the stricter filter is next)
  dml_filtered <- callDML(dml_test, delta = 0.05, p.threshold = 1)
  # Filter DML: fdr < 0.05 & |diff.Methy| ≥ 0.1
  # -- the actual significance + effect-size cutoff used to call a site a real DML
  dml_fdr <- subset(dml_filtered, fdr < 0.05 & abs(diff) >= 0.1)


  # Save raw and filtered DML
  write.table(dml_filtered, file = paste0(g1, "_vs_", g2, "5mC_DML_all.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(dml_fdr, file = paste0(g1, "_vs_", g2, "_5mC_DML_filtered.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

			  
  # Call DMRs
  # -- groups significant, nearby DML sites into regions: merges sites within
  # 250bp (dis.merge), requires >=3 CpGs per region (minCG), >=50% of CpGs in
  # the region significant (pct.sig), and a minimum region length of 20bp (minlen)
  dmr <- callDMR(dml_test,
                 p.threshold = 0.05,
                 delta = 0.05,
                 dis.merge = 250,
                 minCG = 3,
                 pct.sig = 0.5,
                 minlen = 20)

  # Save DMR results
  dmr_file <- paste0(g1, "_vs_", g2, "_5mC_DMR.csv")
  write.csv(dmr, file = dmr_file, row.names = FALSE, quote = FALSE)

  message("Finished: ", g1, " vs ", g2)
}
