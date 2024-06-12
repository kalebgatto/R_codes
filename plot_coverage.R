#!/usr/bin/Rscript
### This script is made to take data from a bam file resulted from a read mapping approach against a reference genome
### and will calculte the coverage and will find windows of 1 Mb in the reference genome that shows significant high 
### number of mapped reads. It will compare the coverage of each 1 Mb windows against each other and will output a
### a plot highlighting only these windows that showed an enrichment of mapped reads.
# Load necessary packages
if (!require(install.packages(c("data.table", "dplyr", "BiocManager", quietly=TRUE))))
if (!require(BiocManager::install(c("karyoploteR", "GenomicRanges", "GenomicAlignments", "GenomeInfoDb", "Rsamtools", "regioneR", quietly=TRUE))))

library(Rsamtools)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomeInfoDb)
library(dplyr)
library(regioneR)
library(karyoploteR)

# Specify the BAM file path
bam_file <- "~/Dados_R/Pxen_plots/Pxen_mapping_X_laevis_genome_sorted.bam"

# Read the BAM file
bam <- BamFile(bam_file)

# Define function to calculate coverage
calculate_coverage <- function(bam, window_size = 1e6) {
  # Get reads from BAM file
  reads <- readGAlignments(bam)
  
  # Get chromosome lengths
  seq_lengths <- seqlengths(bam)
  
  # Create GRanges object with 1 Mb windows
  windows <- tileGenome(seq_lengths, tilewidth = window_size, cut.last.tile.in.chrom = TRUE)
  
  # Calculate coverage
  coverage <- coverage(reads)
  
  # Summarize coverage over windows
  window_coverage <- binnedAverage(windows, coverage, "coverage")
  
  return(window_coverage)
}

# Calculate coverage with 1 Mb windows
window_coverage <- calculate_coverage(bam)

# Function to identify significant windows based on coverage
find_significant_windows <- function(window_coverage, threshold = 0.01) {
  # Ensure coverage values are nonnegative and finite, and convert to integers
  coverage_values <- as.integer(window_coverage$coverage)
  coverage_values[is.na(coverage_values) | coverage_values < 0] <- 0
  
  # Perform statistical test (e.g., Poisson test)
  p_values <- sapply(coverage_values, function(x) poisson.test(x)$p.value)
  
  # Adjust p-values for multiple testing (e.g., Benjamini-Hochberg method)
  adjusted_p_values <- p.adjust(p_values, method = "BH")
  
  # Find significant windows
  significant_windows <- window_coverage[adjusted_p_values <= threshold]
  
  return(significant_windows)
}

# Identify significant windows
significant_windows <- find_significant_windows(window_coverage)

# Function to create a custom plot using karyoploteR
plot_coverage <- function(window_coverage, significant_windows, bam_file) {
  kp <- plotKaryotype(genome="hg19")  # Adjust genome as needed
  density <- kpPlotBAMDensity(kp, data = bam_file, col="blue", window.size = 1e6, r0 = 0, r1 = 0.5)
  kpAxis(density, ymax = density$latest.plot$computed.values$max.density, cex = 0.66)
  kpAddBaseNumbers(density, tick.dist = 20000000)
  
  # Add significant regions
  significant_gr <- as(significant_windows, "GRanges")
  kpPlotRegions(kp, data = significant_gr, col = "red", r0 = 0.6, r1 = 0.7)
}

# Plot the coverage and significant regions
plot_coverage(window_coverage, significant_windows, bam_file)
