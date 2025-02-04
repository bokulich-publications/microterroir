#!/usr/bin/env Rscript

# This is specifically for calculating PERMANOVA for LC-MS 

# Load required libraries
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: permanova.R <distance_matrix.tsv> <metadata.csv> <formula>")
}

# Parse arguments
dist_file <- args[1]
meta_file <- args[2]
formula_string <- args[3]  # Formula passed as a string

# Load distance matrix
dist_matrix <- as.dist(read.delim(dist_file, row.names = 1, check.names = FALSE))

# Load metadata
metadata <- read_csv(meta_file, col_types = cols())

# Ensure row names match in both data
common_samples <- intersect(rownames(as.matrix(dist_matrix)), metadata$sample_id)
dist_matrix <- as.dist(as.matrix(dist_matrix)[common_samples, common_samples])
metadata <- metadata[metadata$sample_id %in% common_samples, ]

# Construct formula properly
formula <- as.formula(paste("~", formula_string))  # Keep only metadata terms

# Run PERMANOVA with by-term breakdown
permanova_result <- adonis2(dist_matrix ~ ., data = metadata[, all.vars(formula)], 
                            permutations = 999, by = "margin")

# Extract and print results as a table
result_df <- as.data.frame(permanova_result)[, c("R2", "Pr(>F)")]
result_df <- result_df %>%
  mutate(Term = rownames(permanova_result)) %>%
  select(Term, R2, `Pr(>F)`)  # Ensure correct order

# Print output in a readable format
cat("\nPERMANOVA Results:\n")
print(result_df, row.names = FALSE)
