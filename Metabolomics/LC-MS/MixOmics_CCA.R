# MixOmics 

#BiocManager::install("mixOmics")

library(mixOmics) 

setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/Metabolomics/Multiomics_preparedData/PostMLF21_SensoryMetabolitesMicrobiome/')


### new try to read in 

# Read the CSV without using the first row as column names
df_raw <- read.csv("MFA_table.csv", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
# Extract the group names from the first row
groups <- as.character(df_raw[1, ])
# Extract the variable names from the second row
variables <- as.character(df_raw[2, ])
# Remove the first two rows (which are now stored separately)
df_clean <- df_raw[-c(1, 2), ]
# Assign proper column names
colnames(df_clean) <- variables
# Convert back to numeric
df_clean[] <- lapply(df_clean, as.numeric)
# Assign the group information back
group_info <- data.frame(Variable = variables, Group = groups)
# Split the data into separate data frames by group
data_list <- split.default(df_clean, groups)

# Split the data into separate data frames for each group
groups <- c("LCMS_Pos", "LCMS_Neg", "Wine Chemistry", "Plots", "Climate", "Fungi", "Bacteria", "GCMS", "Sensory")

# Combine metabolomics data
metabolomics <- list(
  LCMS_Pos = data_list[["LCMS_Pos"]],
  LCMS_Neg = data_list[["LCMS_Neg"]],
  GCMS = data_list[["GCMS"]])

# Combine microbiome data
microbiome <- list(
  Fungi = data_list[["Fungi"]],
  Bacteria = data_list[["Bacteria"]])

# Combine metadata
metadata <- list(
  Wine_Chemistry = data_list[["Wine Chemistry"]],
  Plots = data_list[["Plots"]],
  Climate = data_list[["Climate"]],
  Sensory = data_list[["Sensory"]])




