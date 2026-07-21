#!/usr/bin/env Rscript

# Author: Kevin Boyd
# Date: March 22, 2025
# Purpose: Filter merged peak set to include only those peaks that overlap with a minimum
# number of replicate peaks, and write the result as a BED file.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(stringr)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript filterConsensusFromMergedPeakSet.R merged.narrowPeak rep1.narrowPeak rep2.narrowPeak ... min_overlap out.bed")
}

merged_peak_file <- args[1]
replicate_files  <- args[2:(length(args) - 2)]
min_overlap      <- as.numeric(args[length(args) - 1])
output_bed       <- args[length(args)]

# ---- Helper function to read narrowPeak files as GRanges ----

narrowpeak_colnames <- c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
narrowpeak_colclasses <- c("character", "integer", "integer", "character", "integer", "character",
                           "numeric", "numeric", "numeric", "integer")

read_narrowpeak_as_granges <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t", col.names = narrowpeak_colnames, colClasses = narrowpeak_colclasses)
  # Clean up strand: replace '.' with '*'
  dt$strand[!dt$strand %in% c("+", "-", "*")] <- "*"
  
  # narrowPeak uses 0-based BED starts; GRanges uses 1-based starts.
  gr <- makeGRangesFromDataFrame(
    as.data.frame(dt),
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    starts.in.df.are.0based = TRUE,
    keep.extra.columns = FALSE
  )
  gr
}

# ---- Load merged and replicate peak sets ----

merged_gr <- read_narrowpeak_as_granges(merged_peak_file)
replicate_gr_list <- lapply(replicate_files, read_narrowpeak_as_granges)

# ---- For each merged peak, count how many replicates it overlaps ----

overlap_counts <- sapply(replicate_gr_list, function(rep_gr) {
  overlaps <- countOverlaps(merged_gr, rep_gr)
  as.numeric(overlaps > 0)
})

# Sum across all replicates to get total overlap support
support_counts <- rowSums(overlap_counts)

# Filter peaks supported by at least min_overlap replicates
consensus_gr <- merged_gr[support_counts >= min_overlap]

# ---- Write filtered peaks to output BED file ----

# Create data.frame for writing as BED (chrom, start, end)
# Convert the internal 1-based GRanges coordinates back to BED.
consensus_bed <- data.frame(
  chrom = as.character(seqnames(consensus_gr)),
  start = start(consensus_gr) - 1L,
  end   = end(consensus_gr)
)

# Write without quotes and no row/column names
fwrite(consensus_bed, file = output_bed, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

message("Consensus peaks written to: ", output_bed)
