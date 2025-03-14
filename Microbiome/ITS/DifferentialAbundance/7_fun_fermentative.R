# ABUNDANCE FERMENTATIVE ORGANISMS ON BERRIES 

set.seed(42) # for reproducible stochastic processes
library(phyloseq)
library(ANCOMBC)
library(ComplexHeatmap)
library(pheatmap)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) # for combining multiple plots
library(microViz)
library(tidyverse)
library(stringr)
library(microbiome)
library(RColorBrewer)
library(viridis)
library(legendary)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});","}")))


setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/FermentativeYeasts')


### 1. LOAD DATA AND TURN TO PHYLOSEQ OBJECT 
# load csv (created on the cloud in 6_fun_investigateTaxa.ipynb)
otu_data <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/berries_harvest_ASVs_rarefied_labled.tsv", sep = '\t', row.names = 1, check.names = FALSE)

# Create the taxonomy table by separating taxonomy strings
tax_table <- data.frame(Taxon = colnames(otu_data)) %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies"), 
           sep = ";", fill = "right")

# Ensure otu_table is a matrix
otu_table <- as.matrix(otu_data)
# Make the row names of tax_table unique
rownames(tax_table) <- make.unique(colnames(otu_table))
# Reassign row names of tax_table to the columns of otu_table
colnames(otu_table) <- rownames(tax_table)
# Create the phyoseq object
otu_table_phy <- otu_table(otu_table, taxa_are_rows = FALSE)
tax_table_phy <- tax_table(as.matrix(tax_table))

# add metadata 
metadata_climate <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/Metadata/ITS_Lavaux_Climate.tsv", row.names = 2, sep = "\t", check.names = FALSE)
metadata_chemistry <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/Metadata/ITS_Lavaux_BerryChemistry.tsv", row.names = 2, sep = "\t", check.names = FALSE)

metadata_climate <- metadata_climate %>%
  select(Plot_ID, Year, Altitude, Average_slope, Average_radiation, 
         #Soil_thickness, Soil_type, Soil_depth, Hydromorphie, Geology,
          median_rh, median_temperature
         #maximum_rh, minimum_rh, cv_rh, GDD,maximum_temperature, minimum_temperature, cv_temperature
        )

metadata_chemistry <- metadata_chemistry %>%
  select(Tartrate_gL, Malate_gL, Glucose_gL, Fructose_gL)
         
metadata <- merge(metadata_climate, metadata_chemistry, by = "row.names", all = TRUE)
rownames(metadata) <- metadata$Row.names  # Set row names to the sample IDs
metadata <- metadata[, -1]  # Remove the 'Row.names' column

sample_data_phy <- sample_data(metadata)

# Merge into a Phyloseq Object
physeq <- phyloseq(otu_table_phy, tax_table_phy, sample_data_phy)

# drop empty cols 
physeq <- physeq %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)

# some unidentified Genera etc.
physeq <- physeq %>% tax_fix(unknowns = c("s__unidentified", "g__unidentified", "f__unidentified"))


############################################################     
# Filter to fermentative organisms 
# Define the genera of interest with "g__" prefix
genera_of_interest <- c("Saccharomyces", "Saccharomycodes", "Schizosaccharomyces",
                        "Candida", "Torulaspora", "Debaryomyces", "Issatchenkia",
                        "Pichia", "Kluyveromyces", "Metschnikowia", "Hanseniaspora",
                        "Kloeckera", "Rhodotorula", "Brettanomyces", "Dekkera",
                        "Zygosaccharomyces")
genera_of_interest <- paste0("g__", genera_of_interest)
# Filter phyloseq object to keep only selected genera
physeq_yeasts <- subset_taxa(physeq, Genus %in% genera_of_interest)


############################################################     

### 2. HEATMAP 
cols <- setNames(viridis::viridis_pal()(length(unique(samdat_tbl(physeq)$Year))), 
                 unique(samdat_tbl(physeq)$Year))

# log2
physeq_yeasts %>%
  tax_transform("log2", rank = "Genus", zero_replace = 1) %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_widtsh = 0.3, size = grid::unit(1, "cm"), ylim = 0:1)
    ),
    colors = heat_palette(palette='Purples', rev=TRUE),
    sample_anno = sampleAnnotation(
      Year = anno_sample("Year"),
      col = list(Year = cols), border = TRUE
    )
  )

# presence / absence 
physeq_yeasts %>%
  tax_transform("binary", rank = "Genus") %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_widtsh = 0.3, size = grid::unit(1, "cm"), ylim = 0:1)
    ),
    colors = heat_palette(palette='Purples', rev=TRUE),
    sample_anno = sampleAnnotation(
      Year = anno_sample("Year"),
      col = list(Year = cols), border = TRUE
    )
  )

# Calculate prevalence at the Genus level
physeq_genus <- tax_glom(physeq_yeasts, taxrank = "Genus")  # Aggregate at Genus level
otu_genus <- abundances(physeq_genus)  # Extract OTU table
# Convert to binary presence/absence
otu_binary <- (otu_genus > 0) * 1  
# Calculae prevalence (fraction of samples where each genus is present)
prev_df <- data.frame(
  Genus = tax_table(physeq_genus)[, "Genus"],
  Prevalence = rowSums(otu_binary) / ncol(otu_binary))
print(prev_df)
write.csv(prev_df, "FermYeast_Genus_Prevalence.csv", row.names = FALSE)

# Genus Prevalence 
# g__Saccharomyces 0.340694006
# g__Hanseniaspora 0.735015773
# g__Rhodotorula 0.003154574
# g__Torulaspora 0.009463722
# g__Candida 0.003154574


############################################################     
 
### 2. ANCOM-BC MUTLIVARIATE 

# remove samples with missing data in any metadata column:
physeq_yeast_clean <- subset_samples(physeq_yeasts, complete.cases(sample_data(physeq)))
dim(physeq_yeasts@otu_table)
# 340  32
dim(physeq_yeast_clean@otu_table)
# 314  32
 
output <- ancombc2(
  data = physeq_yeast_clean, 
  tax_level = 'Species',
  fix_formula = "Altitude + Average_slope + Average_radiation + median_rh + median_temperature + Tartrate_gL + Malate_gL + Glucose_gL + Fructose_gL", 
  rand_formula = NULL,  
  pseudo_sens = TRUE,  
  p_adj_method = "holm", 
  prv_cut = 0.1, # note, setting this to 0 keeps all yeasts but then we cannot control for year anymore  
  s0_perc = 0.05,
  struc_zero = FALSE, 
  neg_lb = FALSE,
  alpha = 0.05, 
  n_cl = 2, 
  verbose = TRUE,
  global = FALSE,  # interested at specific taxa not at global comparison 
  pairwise = FALSE,  # Only needed for categorical variables
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
)


# resuls dataframe
res_df <- output$res

# Select lfc columns, excluding lfc_(Intercept)
lfc_cols <- grep("^lfc_", colnames(res_df), value = TRUE)
lfc_cols <- setdiff(lfc_cols, "lfc_(Intercept)")  # Remove lfc_(Intercept)

lfc_matrix <- as.matrix(res_df[, lfc_cols])
rownames(lfc_matrix) <- res_df$taxon

# Identify corresponding q_ values
q_cols <- gsub("^lfc_", "q_", lfc_cols)
q_matrix <- as.matrix(res_df[, q_cols])

# Identify differential abundance locations
diff_cols <- gsub("^lfc_", "diff_", lfc_cols)
diff_matrix <- as.matrix(res_df[, diff_cols])

# Format q-values with 3 decimals and add a star for differentially abundant values
annotation_matrix <- matrix("", nrow = nrow(lfc_matrix), ncol = ncol(lfc_matrix))
for (i in seq_len(nrow(q_matrix))) {
  for (j in seq_len(ncol(q_matrix))) {
    formatted_q <- sprintf("%.3f", as.numeric(q_matrix[i, j]))  # Format to 3 decimals
    if (diff_matrix[i, j] == "TRUE") {
      formatted_q <- paste0(formatted_q, "*")  # Append star if differentially abundant
    }
    annotation_matrix[i, j] <- formatted_q
  }}

# Function to lighten a color
lighten_color <- function(color, factor = 0.5) {
  col <- col2rgb(color) / 255  
  col <- col + (1 - col) * factor  
  rgb(col[1], col[2], col[3], maxColorValue = 1)}
# Define colors and lighten them
low_color <- lighten_color("#39558CFF", 0.5)  # Lighter blue
high_color <- lighten_color("#B8DE29FF", 0.5) # Lighter yellow

# Create a color gradient with white at zero
# Define a symmetric range around zero
limit <- max(abs(lfc_matrix))  # Get the largest absolute value
# Create breaks symmetrically around zero
breaks <- seq(-limit, limit, length.out = 100)
# Generate colors with white at 0
color_palette <- colorRampPalette(c(low_color, "white", high_color))(length(breaks) - 1)

# Remove 'lfc_' from column names
colnames(lfc_matrix) <- gsub("^lfc_", "", colnames(lfc_matrix))

# Transpose the heatmap data
lfc_matrix <- t(lfc_matrix)
annotation_matrix <- t(annotation_matrix)

# Adjust legend size
legend_breaks <- seq(-limit, limit, length.out = 5)  # Fewer legend ticks
legend_labels <- round(legend_breaks, 2)  # Keep readable numbers

# Plot heatmap
pheatmap(
  lfc_matrix, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  display_numbers = annotation_matrix,
  fontsize_number = 10, 
  color = color_palette, 
  breaks = breaks,  # Ensure white is centered at 0
  legend_breaks = legend_breaks,  
  legend_labels = legend_labels,
  legend_title = "LFC"  # Add legend title
)


##############################

### 2. ANCOM-BC SEPERATE MODELS  

# Variables to iterate over
variables <- c("Altitude", "Average_slope" , "Average_radiation","median_rh", "median_temperature", 
               "Tartrate_gL", "Malate_gL", "Glucose_gL", "Fructose_gL")

# Initialize lists to store results
lfc_list <- list()
q_list <- list()
diff_list <- list()

tax_names <- NULL  # Placeholder for taxon names

for (var in variables) {
  output <- ancombc2(
    data = physeq_yeast_clean, 
    tax_level = 'Species',
    fix_formula = var, 
    rand_formula = NULL,  
    pseudo_sens = TRUE,  
    p_adj_method = "holm", 
    prv_cut = 0.1,  
    s0_perc = 0.05,
    struc_zero = FALSE, 
    neg_lb = FALSE,
    alpha = 0.05, 
    n_cl = 2, 
    verbose = TRUE,
    global = FALSE,  
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 100),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
  )
  
  res_df <- output$res
  
  if (is.null(tax_names)) tax_names <- res_df$taxon  # Store taxon names
  
  lfc_list[[var]] <- res_df[[paste0("lfc_", var)]]
  q_list[[var]] <- res_df[[paste0("q_", var)]]
  diff_list[[var]] <- res_df[[paste0("diff_", var)]]
}

# Combine results into matrices
lfc_matrix <- do.call(cbind, lfc_list)
q_matrix <- do.call(cbind, q_list)
diff_matrix <- do.call(cbind, diff_list)
rownames(lfc_matrix) <- tax_names
rownames(q_matrix) <- tax_names
rownames(diff_matrix) <- tax_names

# Format q-values with stars for differentially abundant values
annotation_matrix <- matrix("", nrow = nrow(q_matrix), ncol = ncol(q_matrix))
for (i in seq_len(nrow(q_matrix))) {
  for (j in seq_len(ncol(q_matrix))) {
    formatted_q <- sprintf("%.3f", as.numeric(q_matrix[i, j]))
    if (diff_matrix[i, j] == TRUE) {
      formatted_q <- paste0(formatted_q, "*")
    }
    annotation_matrix[i, j] <- formatted_q
  }
}

# Define a color palette for the heatmap
lighten_color <- function(color, factor = 0.5) {
  col <- col2rgb(color) / 255  
  col <- col + (1 - col) * factor  
  rgb(col[1], col[2], col[3], maxColorValue = 1)
}

low_color <- lighten_color("#39558CFF", 0.5)
high_color <- lighten_color("#B8DE29FF", 0.5)

limit <- max(abs(lfc_matrix), na.rm = TRUE)
breaks <- seq(-limit, limit, length.out = 100)
color_palette <- colorRampPalette(c(low_color, "white", high_color))(length(breaks) - 1)

# Transpose matrices for heatmap
lfc_matrix <- t(lfc_matrix)
annotation_matrix <- t(annotation_matrix)

# Plot heatmap
pheatmap(
  lfc_matrix, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  display_numbers = annotation_matrix,
  fontsize_number = 10, 
  color = color_palette, 
  breaks = breaks,
  legend_breaks = seq(-limit, limit, length.out = 5),
  legend_labels = round(seq(-limit, limit, length.out = 5), 2),
  legend_title = "LFC"
)
