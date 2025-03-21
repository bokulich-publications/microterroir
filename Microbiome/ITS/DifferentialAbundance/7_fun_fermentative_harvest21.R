# ABUNDANCE FERMENTATIVE ORGANISMS ON BERRIES (Harvest 2021 only)

set.seed(42) 
library(phyloseq)
library(ANCOMBC)
library(ComplexHeatmap)
library(pheatmap)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
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
otu_data <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/berries_harvest21_cOTUs_rarefied_labled.tsv", sep = '\t', row.names = 1, check.names = FALSE,encoding = "UTF-8")
dim(otu_data)


# Create the taxonomy table by separating taxonomy strings
tax_table <- data.frame(Taxon = colnames(otu_data)) %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
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
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

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
    colors = heat_palette(palette='Purples', rev=TRUE)
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
# 218   6
dim(physeq_yeast_clean@otu_table)
# 218   6

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
  neg_lb = TRUE,
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


plot_lfc_heatmap <- function(res_df) {
  # Select lfc columns, excluding lfc_(Intercept)
  lfc_cols <- grep("^lfc_", colnames(res_df), value = TRUE)
  lfc_cols <- setdiff(lfc_cols, "lfc_(Intercept)")
  
  lfc_matrix <- as.matrix(res_df[, lfc_cols])
  rownames(lfc_matrix) <- res_df$taxon
  
  # Identify corresponding q_ and diff_ values
  q_matrix <- as.matrix(res_df[, gsub("^lfc_", "q_", lfc_cols)])
  diff_matrix <- as.matrix(res_df[, gsub("^lfc_", "diff_", lfc_cols)])
  
  # Format q-values and mark differentially abundant values
  annotation_matrix <- matrix(mapply(
    function(q, diff) ifelse(diff == "TRUE", paste0(sprintf("%.3f", as.numeric(q)), "*"), sprintf("%.3f", as.numeric(q))),
    q_matrix, diff_matrix
  ), nrow = nrow(lfc_matrix), ncol = ncol(lfc_matrix))
  # Function to lighten a color
  lighten_color <- function(color, factor = 0.5) {
    col <- col2rgb(color) / 255  
    col <- col + (1 - col) * factor  
    rgb(col[1], col[2], col[3], maxColorValue = 1)  }
  # Define colors
  low_color <- lighten_color("#39558CFF", 0.5)  # Lighter blue
  high_color <- lighten_color("#B8DE29FF", 0.5) # Lighter yellow
  
  # Define color scale
  limit <- max(abs(lfc_matrix))  
  breaks <- seq(-limit, limit, length.out = 100)
  color_palette <- colorRampPalette(c(low_color, "white", high_color))(length(breaks) - 1)
  
  # Clean column names
  colnames(lfc_matrix) <- gsub("^lfc_", "", colnames(lfc_matrix))
  
  # Transpose matrices
  lfc_matrix <- t(lfc_matrix)
  annotation_matrix <- t(annotation_matrix)
  
  # Define legend
  legend_breaks <- seq(-limit, limit, length.out = 5)
  legend_labels <- round(legend_breaks, 2)
  
  # Generate heatmap
  pheatmap(
    lfc_matrix, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    display_numbers = annotation_matrix,
    fontsize_number = 10, 
    color = color_palette, 
    breaks = breaks,
    legend_breaks = legend_breaks,  
    legend_labels = legend_labels,
    legend_title = "LFC"
  )
}

# PLOT
plot_lfc_heatmap(output$res)




