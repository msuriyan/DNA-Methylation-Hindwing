# =============================================================================
# Rscript_BsObjM_paired_All.R
#
# Purpose: build a single bsseq "BSseq" object holding 5mC methylation counts
# for all 16 samples, from the per-sample DSS-format bed files produced earlier
# in the pipeline (frequency-filtered, >10%, coverage >=5). This BSseq object is
# the input DSS_5mC_Freq10_pairwise_Rscript.R reads back in to run pairwise
# DMLtest/callDML/callDMR comparisons.
#
# Each input file must be in DSS's required 4-column format (see example below).
#
# Output: BSseq_5mC_Freq10_all_samples.rds -- one R object file containing
# every sample's methylation counts, ready for downstream DSS testing.
# =============================================================================

library(DSS)
library(bsseq)
library(parallel)

## 5mC
setwd("/hpctmp/snm_02/Methylome/Freq10")

# List your sample files
# -- one DSS-format bed per sample (16 samples total: 8 groups x 2 replicates)
files_5mC <- c("D15_3_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "D15_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "D50_1_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "D50_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "D60_3_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "D60_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "DP_1_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "DP_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "P15_1_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "P15_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "P50_1_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "P50_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "W60_2_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "W60_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "WP_3_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed",
           "WP_M_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed")

# Sample names
# -- strips the filename suffix, leaving just the sample ID (e.g. "D15_3") to
# label each sample inside the BSseq object
sampleNames_5mC <- gsub("_BANY1.2_noH_pileup_Cov5_m_Freq10_DSS.bed", "", files_5mC)

# Read and build BSseq object
# -- reads every file in files_5mC into a table (each has columns chr, pos, N,
# X -- see format note below), then makeBSseqData() combines all 16 samples
# into one BSseq object, aligning sites across samples by chr/pos
BSobj_5mC <- makeBSseqData(lapply(files_5mC, read.table, header = TRUE), sampleNames_5mC)

# View summary
# -- prints the BSseq object's dimensions/structure to console as a sanity check
BSobj_5mC

# Save the BSseq object
# -- written once here so the pairwise comparison script (DSS_5mC_Freq10_pairwise_Rscript.R)
# can just readRDS() this instead of rebuilding it from the raw bed files every time
saveRDS(BSobj_5mC, file = "BSseq_5mC_Freq10_all_samples.rds")
