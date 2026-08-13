###############################################
# Script name: makeEulerPlotOfBedOverlaps.R
# Author: Kevin Boyd
# Updated: August 13, 2026
# Purpose: Generate an Euler plot of overlapping BED files.
###############################################

library(data.table)
library(GenomicRanges)
library(eulerr)
library(ggplotify)
library(magrittr)
library(stringr)

# Parse command-line args
args <- commandArgs(trailingOnly = TRUE)
bed_files <- unlist(strsplit(args[1], ","))  
set_names <- unlist(strsplit(args[2], ","))  
output_rds <- args[3]
output_pdf <- args[4]
font_size <- as.numeric(args[5])
colors <- strsplit(args[6], ",")[[1]]   
pdf_width <- as.numeric(args[7])
pdf_height <- as.numeric(args[8])

# Function to read a BED/narrowPeak file as GRanges (handles . strand)
read_bed_gr <- function(bed_file) {
  dt <- fread(bed_file, header = FALSE)
  dt[[2]] <- as.integer(round(as.numeric(dt[[2]])))
  dt[[3]] <- as.integer(round(as.numeric(dt[[3]])))

  # Sanitize strand column (convert '.' → '*')
  strand_col <- if (ncol(dt) >= 6) {
    strand_vals <- as.character(dt[[6]])
    strand_vals[!strand_vals %in% c("+", "-", "*")] <- "*"
    strand_vals
  } else {
    "*"
  }

  gr <- GRanges(
    seqnames = dt[[1]],
    ranges = IRanges(start = dt[[2]], end = dt[[3]]),
    strand = strand_col
  )
  # BED intervals are 0-based, half-open; GRanges intervals are
  # 1-based, closed. Add 1 only to the BED start coordinate.
  gr <- GRanges(
    seqnames = dt[[1]],
    ranges = IRanges(
      start = dt[[2]] + 1L,
      end = dt[[3]]),
    strand = strand_col)

  if (ncol(dt) > 3) {
    mcols(gr) <- dt[, -(1:3)]
  }

  return(gr)
}

# Read all BEDs into GRanges
gr_list <- lapply(bed_files, read_bed_gr)
names(gr_list) <- set_names

# Create a unified GRanges set
all_combined <- do.call("c", unname(gr_list))
all_union <- reduce(all_combined)

# Create presence/absence matrix
mcols(all_union) <- do.call(
  cbind,
  lapply(seq_along(gr_list), function(i) {
    overlaps <- findOverlaps(all_union, gr_list[[i]])
    in_set <- rep(0, length(all_union))
    in_set[queryHits(overlaps)] <- 1
    in_set
  })
)
colnames(mcols(all_union)) <- set_names

# Build and plot Euler
presence_absence_mat <- as.matrix(mcols(all_union))
eulerr_options(
  labels = list(fontsize = font_size),
  quantities = list(fontsize = font_size - 2, padding = grid::unit(100, "mm")),
  legend = list(fontsize = font_size, vgap = 0.01)
)

EulerPlot <- presence_absence_mat %>%
  euler(shape = "ellipse") %>%
  plot(quantities = TRUE, legend = TRUE, adjust_labels = TRUE, fills = colors) %>%
  as.ggplot()

# Save both outputs
saveRDS(EulerPlot, file = output_rds)
pdf(output_pdf, width = pdf_width, height = pdf_height)
print(EulerPlot)
dev.off()
