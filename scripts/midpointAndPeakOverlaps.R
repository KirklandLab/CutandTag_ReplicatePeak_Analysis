# midpointAndPeakOverlaps.R
# Author: Kevin Boyd
# Purpose: Generate midpoint and full overlap BED files from a list of reproducible peaks,
#          and additionally compute unique sample midpoints (unique to only one sample).
# Date Modified: 7/21/2026

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# All but the last argument are input BED files; the last argument is the output directory for overlaps.
input_beds <- args[-length(args)]
output_dir <- args[length(args)]

# Ensure the output directory ends with a slash
if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
}

# -----------------------------
# Original Overlap and Midpoint Calculation
# -----------------------------

# Import BED files using BED-aware coordinate conversion.
Beds.gr <- lapply(input_beds, function(path) {
    import(path, format = "BED")
})

# Merge genuinely overlapping ranges, but do not merge ranges that
# are merely adjacent under BED interval semantics.
AllBeds.gr <- reduce(
    do.call("c", Beds.gr),
    min.gapwidth = 0L,
    ignore.strand = TRUE
)

# Generate overlap metadata:
# For each input GRanges, check which ranges in AllBeds.gr overlap it.
# Combine the resulting logical vectors into a data frame and attach as metadata.
peakOverlaps <- lapply(Beds.gr, function(gr) { AllBeds.gr %over% gr }) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  # Use the file basenames (with "_consensus_peaks.bed" removed) as column names
  set_names(gsub("_consensus_peaks.bed", "", basename(input_beds))) %>%
  { mcols(AllBeds.gr) <- .; AllBeds.gr }

# Helper before calculating midpoints
make_one_base_midpoints <- function(gr) {
    if (length(gr) == 0L) {
        return(gr)
    }

    # Convert the GRanges start back to a BED start, calculate the
    # conventional BED midpoint, then construct a width-one GRanges.
    bed_start <- start(gr) - 1L
    bed_midpoint <- floor((bed_start + end(gr)) / 2)

    ranges(gr) <- IRanges(
        start = bed_midpoint + 1L,
        width = rep.int(1L, length(gr))
    )

    gr
}

# Calculate midpoints for the consensus overlaps.
midpointOverlaps <- make_one_base_midpoints(peakOverlaps)

# Export the overall midpoint and full overlaps files to the specified output directory.
export.bed(midpointOverlaps, paste0(output_dir, "MidpointOverlaps.bed"))
export.bed(peakOverlaps, paste0(output_dir, "PeakOverlaps.bed"))

# -----------------------------
# Additional: Calculate Unique Sample Midpoints
# -----------------------------
# For each sample (i.e. each column in the overlap metadata), extract the consensus peaks
# that are contributed by only that sample (where the row sum is 1 and that sample is TRUE).
# Then recalculate the midpoint for those peaks and export them as a BED file in the consensus peaks folder

# Assume the consensus peaks folder is the same as the folder of the input BED files.
consensus_dir <- dirname(input_beds[[1]])

# Convert the metadata (overlap matrix) to a data frame for easy subsetting.
overlap_df <- as.data.frame(mcols(peakOverlaps))

# Loop over each sample column.
for(sample in colnames(overlap_df)) {
    # Find indices where this sample is TRUE and the total number of TRUE values in that row is 1.
    unique_idx <- which(overlap_df[[sample]] & (rowSums(overlap_df) == 1))
    
    # Subset the consensus GRanges to obtain peaks unique to this sample.
    unique_peaks <- make_one_base_midpoints(peakOverlaps[unique_idx])

    # Define the output filename. The file will be placed in the consensus peaks folder.
    out_file <- file.path(consensus_dir, paste0(sample, "_unique_MP.bed"))
    
    # Export the unique midpoints BED file.
    export.bed(unique_peaks, out_file)
}
