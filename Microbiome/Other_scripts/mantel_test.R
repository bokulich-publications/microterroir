#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(vegan))

# Parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript mantel_analysis.R <folder_path_1> <folder_path_2> <geodesic_matrix>")
}

folder_path_1 <- args[1]
folder_path_2 <- args[2]
geodesic_matrix <- args[3]

# Read the geodesic distance matrix
matrix_spatial <- read.table(geodesic_matrix, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# Function to perform Mantel test
perform_mantel <- function(micro_matrix_file, spatial_matrix) {
  # Read the microbial distance matrix
  matrix_micro <- read.table(micro_matrix_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

  # Find common IDs
  common_ids <- Reduce(intersect, list(rownames(matrix_micro), rownames(spatial_matrix)))

  if (length(common_ids) == 0) {
    stop("Error: No common IDs found among matrices.")
  }

  # Subset matrices to common IDs
  matrix_micro <- matrix_micro[common_ids, common_ids, drop=FALSE]
  spatial_matrix <- spatial_matrix[common_ids, common_ids, drop=FALSE]

  # Convert to dist objects
  matrix_micro_dist <- as.dist(as.matrix(matrix_micro))
  spatial_matrix_dist <- as.dist(as.matrix(spatial_matrix))

  # Perform Mantel test for Spearman correlation
  mantel_result <- mantel(matrix_micro_dist, spatial_matrix_dist, method="spearman", permutations=999)

  return(list(Mantel_statistic_r = mantel_result$statistic, Significance = mantel_result$signif))
}

# Function to detect available distance matrices in a given folder
find_distance_matrices <- function(folder_path) {
  possible_types <- c(
    "bray_curtis_distance_matrix",
    "jaccard_distance_matrix",
    "unweighted_unifrac_distance_matrix",
    "weighted_unifrac_distance_matrix"
  )
  
  available_types <- possible_types[file.exists(file.path(folder_path, possible_types, "distance-matrix.tsv"))]
  return(available_types)
}

# Initialize result data frame
results <- data.frame(
  Test = character(),
  Mantel_statistic_r = numeric(),
  Significance = numeric(),
  row.names = NULL
)

# Perform Mantel tests for both folder paths
for (folder_path in c(folder_path_1, folder_path_2)) {
  available_distance_types <- find_distance_matrices(folder_path)
  
  for (distance_type in available_distance_types) {
    matrix_path <- file.path(folder_path, distance_type, "distance-matrix.tsv")

    mantel_result <- perform_mantel(matrix_path, matrix_spatial)

    # Append results to the table
    temp_results <- data.frame(
      Test = paste(folder_path, distance_type, sep = "-"),
      Mantel_statistic_r = mantel_result$Mantel_statistic_r,
      Significance = mantel_result$Significance
    )
    results <- rbind(results, temp_results)

    # Write results to a file in the current folder path
    result_file_path <- file.path(folder_path, "mantel_results.tsv")
    write.table(results, file=result_file_path, sep="\t", quote=FALSE, row.names=FALSE)
  }
}

# Print the results table to the console
print(results)
