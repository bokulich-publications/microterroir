#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(vegan))

# Parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript mantel_analysis.R <micro_matrix> <spatial_matrix> <matrix_2>")
}

micro_matrix_file <- args[1]
spatial_matrix_file <- args[2]
matrix_2_file <- args[3]

# Read the distance matrices
matrix_micro <- read.table(micro_matrix_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
matrix_spatial <- read.table(spatial_matrix_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
matrix_2 <- read.table(matrix_2_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# Find common IDs
common_ids <- Reduce(intersect, list(rownames(matrix_micro), rownames(matrix_spatial), rownames(matrix_2)))

if (length(common_ids) == 0) {
  stop("Error: No common IDs found among matrices.")
}

# Subset matrices to common IDs
matrix_micro <- matrix_micro[common_ids, common_ids, drop=FALSE]
matrix_spatial <- matrix_spatial[common_ids, common_ids, drop=FALSE]
matrix_2 <- matrix_2[common_ids, common_ids, drop=FALSE]

# convert to dist objects
matrix_micro_dist <- as.dist(as.matrix(matrix_micro))
matrix_spatial_dist <- as.dist(as.matrix(matrix_spatial))
matrix_2_dist <- as.dist(as.matrix(matrix_2))


# Perform Mantel tests and extract relevant values
mantel_micro_spatial_Pearson <- mantel(matrix_micro_dist, matrix_spatial_dist, method="pearson", permutations=999)
mantel_micro_matrix2_Pearson <- mantel(matrix_micro_dist, matrix_2_dist, method="pearson", permutations=999)
mantel_micro_spatial_spearman <- mantel(matrix_micro_dist, matrix_spatial_dist, method="spearman", permutations=999)
mantel_micro_matrix2_spearman <- mantel(matrix_micro_dist, matrix_2_dist, method="spearman", permutations=999)

mantel_partial <- mantel.partial(matrix_micro_dist, matrix_2_dist, matrix_spatial_dist, method="pearson", permutations=999)
mantel_partial2 <- mantel.partial(matrix_micro_dist, matrix_spatial_dist, matrix_2_dist, method="pearson", permutations=999)


# Extract Mantel statistic r and significance from each result
mantel_results <- data.frame(
  Test = c("Microbial vs Spatial (P)", "Microbial vs Spatial (S)",  "Microbial vs Matrix 2 (P)", "Microbial vs Matrix 2 (S)", "Micro. vs Matrix 2 (ctrl. for Spatial)", "Micro. vs Spatial (ctrl. for Matrix 2)" ),
  Mantel_statistic_r = c(mantel_micro_spatial_Pearson$statistic, mantel_micro_spatial_spearman$statistic, mantel_micro_matrix2_Pearson$statistic, mantel_micro_matrix2_spearman$statistic, mantel_partial$statistic, mantel_partial2$statistic),
  Significance = c(mantel_micro_spatial_Pearson$signif, mantel_micro_spatial_spearman$signif,  mantel_micro_matrix2_Pearson$signif, mantel_micro_matrix2_spearman$signif, mantel_partial$signif, mantel_partial2$signif)
)

print(mantel_results)
