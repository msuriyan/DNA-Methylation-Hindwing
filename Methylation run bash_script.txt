#!/bin/bash

# =========================================================
# Step 1: Dorado basecalling (with 5mC, 5hmC, 6mA models)
# =========================================================

module load dorado/0.3.1-foss-2022a-CUDA-11.7.0

MODEL_DIR=models
INPUT_DIR=raw_data

for i in $(seq -w 1 16); do
elements=barcode${i}
dorado basecaller ${MODEL_DIR}/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 --emit-moves --modified-bases-models ${MODEL_DIR}/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mC@v2,${MODEL_DIR}/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2,${MODEL_DIR}/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_6mA@v2 ${INPUT_DIR}/${elements}/ > ${elements}_calls.bam
done


# =========================================================
# Step 2: Align to reference with minimap2
# =========================================================

module load miniconda
source activate minimap2

REF=reference/BANY_v1.2.fa

for elements in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 barcode13 barcode14 barcode15 barcode16 D15_M D50_M D60_M DP_M P15_M P50_M W60_M WP_M
do
samtools fastq -@96 -T mv,MM,ML ${elements}_calls.bam | minimap2 -ax map-ont --secondary=no ${REF} - -y -t 96 | samtools view -@96 -o ${elements}_BANY1.2_T2.bam
done


# =========================================================
# Step 3: modkit pileup (methylation calling, with 5hmC)
# =========================================================

module load modkit/0.2.8-rc1-GCC-12.3.0

REF=reference/BANY_v1.2.fa
BAMDIR=.
OUTDIR=modkit_pileupOut_withH

mkdir -p "$OUTDIR"

# note: input bam here must already be sorted + indexed (samtools sort / samtools index)

for elements in D15_3 D15_M D50_1 D50_M D60_3 D60_M DP_1 DP_M P15_1 P15_M P50_1 P50_M W60_2 W60_M WP_3 WP_M
do
    echo "Running modkit pileup with h included for ${elements}"

    modkit pileup \
        -t 24 \
        --chunk-size 108 \
        -r "$REF" \
        --filter-threshold A:0.9 \
        --filter-threshold C:0.9 \
        --prefix ${elements}_BANY1.2_withH_A0.9C0.9 \
        "$BAMDIR/${elements}_BANY1.2_sort.bam" \
        "$OUTDIR/${elements}_BANY1.2_withH_pileup0.9.bed" \
        --log-filepath "$OUTDIR/${elements}_withH_pileup_0.9.log"
done


# =========================================================
# Step 4: filter pileup by coverage (>=5)
# =========================================================

for elements in D15_3 D15_M D50_1 D50_M D60_3 D60_M DP_1 DP_M P15_1 P15_M P50_1 P50_M W60_2 W60_M WP_3 WP_M
do
awk '$5 >= 5' ${elements}_BANY1.2_withH_pileup0.9.bed > ${elements}_BANY1.2_withH_pileup0.9_Cov5.bed
done


# =========================================================
# Step 5: filter by modification type (h = 5hmC)
# =========================================================

for elements in D15_3 D15_M D50_1 D50_M D60_3 D60_M DP_1 DP_M P15_1 P15_M P50_1 P50_M W60_2 W60_M WP_3 WP_M
do
#awk '$4 == "a"' ${elements}_BANY1.2_noH_pileup_Cov5.bed > ${elements}_BANY1.2_noH_pileup_Cov5_a.bed
#awk '$4 == "m"' ${elements}_BANY1.2_noH_pileup_Cov5.bed > ${elements}_BANY1.2_noH_pileup_Cov5_m.bed
awk '$4 == "h"' ${elements}_BANY1.2_withH_pileup0.9_Cov5.bed > ${elements}_BANY1.2_withH_pileup0.9_Cov5_h.bed
done


# =========================================================
# Step 6: filter by methylation frequency (m, >10%)
# =========================================================

for elements in D15_3 D15_M D50_1 D50_M D60_3 D60_M DP_1 DP_M P15_1 P15_M P50_1 P50_M W60_2 W60_M WP_3 WP_M
do
awk '$11 > 10.00' ${elements}_BANY1.2_withH_pileup0.9_Cov5_m.bed > ${elements}_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed
done


# =========================================================
# Step 7: pairwise intersect replicate pairs, then merge into "All" bed
# (this produces All_BANY1.2_m_merged2.bed / All_BANY1.2_h_merged2.bed used in step 8)
# =========================================================

# m -- uses the Freq10-filtered files from step 6
bedtools intersect -wa -a D15_3_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b D15_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > D15_m.bed
bedtools intersect -wa -a D50_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b D50_1_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > D50_m.bed
bedtools intersect -wa -a D60_3_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b D60_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > D60_m.bed
bedtools intersect -wa -a DP_1_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b DP_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > DP_m.bed
bedtools intersect -wa -a P15_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b P15_1_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > P15_m.bed
bedtools intersect -wa -a P50_1_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b P50_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > P50_m.bed
bedtools intersect -wa -a W60_2_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b W60_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > W60_m.bed
bedtools intersect -wa -a WP_M_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b WP_3_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed > WP_m.bed

cat D15_m.bed D50_m.bed D60_m.bed DP_m.bed P15_m.bed P50_m.bed W60_m.bed WP_m.bed | sort -k1,1 -k2,2n > All_BANY1.2_m_merged2.bed

# h -- uses the plain Cov5 files (no frequency filter for h)
bedtools intersect -wa -a D15_3_BANY1.2_withH_pileup0.9_Cov5_h.bed -b D15_M_BANY1.2_withH_pileup0.9_Cov5_h.bed > D15_h.bed
bedtools intersect -wa -a D50_M_BANY1.2_withH_pileup0.9_Cov5_h.bed -b D50_1_BANY1.2_withH_pileup0.9_Cov5_h.bed > D50_h.bed
bedtools intersect -wa -a DP_1_BANY1.2_withH_pileup0.9_Cov5_h.bed -b DP_M_BANY1.2_withH_pileup0.9_Cov5_h.bed > DP_h.bed
bedtools intersect -wa -a P15_M_BANY1.2_withH_pileup0.9_Cov5_h.bed -b P15_1_BANY1.2_withH_pileup0.9_Cov5_h.bed > P15_h.bed
bedtools intersect -wa -a W60_2_BANY1.2_withH_pileup0.9_Cov5_h.bed -b W60_M_BANY1.2_withH_pileup0.9_Cov5_h.bed > W60_h.bed
bedtools intersect -wa -a WP_M_BANY1.2_withH_pileup0.9_Cov5_h.bed -b WP_3_BANY1.2_withH_pileup0.9_Cov5_h.bed > WP_h.bed
bedtools intersect -wa -a D60_3_BANY1.2_withH_pileup0.9_Cov5_h.bed -b D60_M_BANY1.2_withH_pileup0.9_Cov5_h.bed > D60_h.bed
bedtools intersect -wa -a P50_1_BANY1.2_withH_pileup0.9_Cov5_h.bed -b P50_M_BANY1.2_withH_pileup0.9_Cov5_h.bed > P50_h.bed

cat D15_h.bed D50_h.bed D60_h.bed DP_h.bed P15_h.bed P50_h.bed W60_h.bed WP_h.bed | sort -k1,1 -k2,2n > All_BANY1.2_h_merged2.bed


# =========================================================
# Step 8: intersect each sample with the merged "All" bed, per mod type
# =========================================================

for elements in D15_3 D15_M D50_1 D50_M D60_3 D60_M DP_1 DP_M P15_1 P15_M P50_1 P50_M W60_2 W60_M WP_3 WP_M
do
    # m -- uses the Freq10-filtered file from step 6
    bedtools intersect -wa -a ${elements}_BANY1.2_withH_pileup0.9_Cov5_Freq10_m.bed -b All_BANY1.2_m_merged2.bed > ${elements}_BANY1.2_withH_pileup0.9_Cov5_Freq10_m_All.bed

    # h -- unchanged, no frequency filter applied; since h (5hmC) is usually low-abundance genome-wide, we capture all h modifications rather than restricting by frequency
	
    bedtools intersect -wa -a ${elements}_BANY1.2_withH_pileup0.9_Cov5_h.bed -b All_BANY1.2_h_merged2.bed > ${elements}_BANY1.2_withH_pileup0.9_Cov5_h_All.bed
done


