#!/usr/bin/Rscript
### This script is made to take data from a bam file resulted from a read mapping approach against a reference genome
### and will calculte the coverage and will find windows of 1 Mb in the reference genome that shows significant high 
### number of mapped reads. It will compare the coverage of each 1 Mb windows against each other and will output a
### a plot highlighting only these windows that showed an enrichment of mapped reads.
# Load necessary packages
if (!require(install.packages(c("data.table", "GenomicRanges", "dplyr", "BiocManager", quietly=TRUE))))
if (!require(BiocManager::install(c("karyoploteR", "GenomicAlignments", quietly=TRUE))))

library(Rsamtools)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomeInfoDb)
library(dplyr)
library(regioneR)
library(karyoploteR)

# Function to read BAM file and calculate coverage
read_and_calculate_coverage <- function(bam_file) {
  bam <- readGAlignments(bam_file)
  coverage <- coverage(bam)
  return(coverage)
}
# Function to calculate coverage in 1 Mb windows
calculate_coverage <- function(bam, window_size = 1e6) {
  reads <- readGAlignments(bam)
  seq_lengths <- seqlengths(bam)
  windows <- tileGenome(seq_lengths, tilewidth = window_size, cut.last.tile.in.chrom = TRUE)
  coverage <- coverage(reads)
  window_coverage <- binnedAverage(windows, coverage, "coverage")
  return(window_coverage)
}
# Function to identify significant windows based on coverage
find_significant_windows <- function(window_coverage, threshold = 0.01) {
  coverage_values <- as.integer(window_coverage$coverage)
  coverage_values[is.na(coverage_values) | coverage_values < 0] <- 0
  
  p_values <- sapply(coverage_values, function(x) poisson.test(x)$p.value)
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  significant_windows <- window_coverage[adjusted_p_values <= threshold]
  
  return(significant_windows)
}
# Function to create a custom plot using karyoploteR
plot_coverage <- function(window_coverage, significant_windows, bam_file, output_file) {
  # Set up the karyoploteR plot
  kp <- plotKaryotype(genome = "your_genome_of_interest") ### Replace 'your_genome_of_interest' by the GRanges object of your genome
  
  # Plot the BAM file density
  kp <- kpPlotBamDensity(kp, data = bam_file, r0 = 0, r1 = 0.5, col = "blue")
  
  # Add significant regions
  significant_gr <- as(significant_windows, "GRanges")
  kpPlotRegions(kp, data = significant_gr, col = "red", r0 = 0.6, r1 = 0.7)
  
  # Add axis for better visualization
  kpAxis(kp, r0 = 0, r1 = 0.5)
  
  # Save the plot as an EPS file
  dev.copy2eps(file = output_file, width = 861 / 72, height = 537 / 72)
  dev.off()
}
# Main function to orchestrate the steps
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Usage: plot_coverage.R <bam_file> <output_file>")
  }
  
  bam_file <- args[1]
  output_file <- args[2]
  
  coverage <- read_and_calculate_coverage(bam_file)
  window_coverage <- calculate_coverage(coverage)
  significant_windows <- find_significant_windows(window_coverage)
  
  plot_coverage(window_coverage, significant_windows, bam_file, output_file)
}

# Run the main function
main()
