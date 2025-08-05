#!/bin/bash
#SBATCH --job-name=generate_heatmap
#SBATCH --cpus-per-task=8
#SBATCH --time=16:00:00
#SBATCH --mem=64G

##################################################################
##                     Variables to Modify                      ##
##################################################################

# ============================================================== #
# 1) Set your BigWig files for each sample name:                 #
#    Syntax: bigwigs[SampleName]="path/to/file.bw"               #
# -------------------------------------------------------------- #
# Example:                                                       #
# bigwigs[Sample1]="path/to/Sample1_consensus_peaks.bw"          #
# bigwigs[Sample2]="path/to/Sample2_consensus_peaks.bw"          #
# bigwigs[Sample3]="path/to/Sample3_consensus_peaks.bw"          #
#                                                                #
# ============================================================== #
# 2) Set your BED file with regions of interest:                 #
# -------------------------------------------------------------- #
# Example:                                                       #
# BED_FILE="results/overlaps/MidpointOverlaps.bed"               #
#                                                                #
# ============================================================== #
# 3) Set output file names:                                      #
# -------------------------------------------------------------- #
# MATRIX_OUT="results/heatmap/Beds_Specific_Samples.gz"          #
# HEATMAP_OUT="results/heatmap/HeatPlots_Specific_Samples.png"   #
#                                                                #
# ============================================================== #
# 4) Set your sample order (1-n: order of samples as presented): #
# -------------------------------------------------------------- #
# SAMPLE_ORDER="1 2 3"                                           #
#                                                                #
# ============================================================== #
# USAGE:                                                         #
# After modifying the variables above,                           #
# submit this script to SLURM with:                              #
#                                                                #
#   sbatch make_heatplots.sh                                     #
# ============================================================== #

##################################################################
##                     User-defined variables                   ##
##################################################################

# ---- Paths to BigWig files ----
# Use associative arrays for clarity
declare -A bigwigs

# Example: sample name => path
bigwigs[H3K27ac_C]="results/consensusPeaks/H3K27ac_C_consensus_peaks.bw"
bigwigs[H3K27ac_T]="results/consensusPeaks/H3K27ac_T_consensus_peaks.bw"
bigwigs[H3K27me3_C]="results/consensusPeaks/H3K27me3_C_consensus_peaks.bw"
bigwigs[H3K27me3_T]="results/consensusPeaks/H3K27me3_T_consensus_peaks.bw"

# ---- Bed file with regions ----
BED_FILE="results/overlaps/MidpointOverlaps.bed"

# ---- Output files ----
MATRIX_OUT="results/heatmap/Beds_Specific_Samples.gz"
HEATMAP_OUT="results/heatmap/HeatPlots_Specific_Samples.png"

# ---- Heatmap parameters ----
BEFORE_REGION=3000
AFTER_REGION=3000
PROCESSORS=8
REFPOINT_NAME="Center"
REGIONS_NAME="All Midpoints"

# ---- Define sample order (space-separated numbers, order matches bigwig input order) ----
SAMPLE_ORDER="1 2 3 4"

##################################################################
##                      Generate input lists                    ##
##################################################################

# Build the list of BigWig files in sample order
BIGWIG_LIST=""
SAMPLES_LABEL=""

for sample in ${SAMPLE_ORDER}; do
    BIGWIG_LIST+="${bigwigs[$sample]} "
    SAMPLES_LABEL+="${sample} "
done

##################################################################
##                        Run computeMatrix                     ##
##################################################################

echo "Running computeMatrix..."
computeMatrix reference-point \
    -S ${BIGWIG_LIST} \
    -R ${BED_FILE} \
    --outFileName ${MATRIX_OUT} \
    -a ${AFTER_REGION} -b ${BEFORE_REGION} \
    --numberOfProcessors ${PROCESSORS}

##################################################################
##                         Run plotHeatmap                      ##
##################################################################

echo "Running plotHeatmap..."
plotHeatmap \
    -m ${MATRIX_OUT} \
    -out ${HEATMAP_OUT} \
    --dpi 1000 \
    --sortUsing sum \
    --sortUsingSamples ${SAMPLE_ORDER} \
    --refPointLabel "${REFPOINT_NAME}" \
    --regionsLabel "${REGIONS_NAME}" \
    --samplesLabel ${SAMPLES_LABEL}

echo "Heatmap generation complete!"
