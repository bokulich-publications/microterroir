# ABUNDANCE FERMENTATIVE ORGANISMS ON BERRIES 

set.seed(42) # for reproducible stochastic processes
library(phyloseq)
library(ANCOMBC)
library(ComplexHeatmap)
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
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))


setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance')


### 1. LOAD DATA AND TURN TO PHYLOSEQ OBJECT 
# load csv (created on the cloud in 6_fun_investigateTaxa.ipynb)
#otu_data <- read.csv("FermentativeOrganisms_crOTU90_annotated.csv", row.names = 1, check.names = FALSE)
otu_data <- read.csv("berries_harvest_ASVs_rarefied_labled.tsv", sep = '\t', row.names = 1, check.names = FALSE)


# to manually make a phyloseq object we create the taxa table
tax_table <- data.frame(Taxon = colnames(otu_data)) %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right")
rownames(tax_table) <- colnames(otu_data)  # Set OTU names as row names

# convert to phyloseq objects 
otu_table_phy <- otu_table(as.matrix(otu_data), taxa_are_rows = FALSE)
tax_table_phy <- tax_table(as.matrix(tax_table))

# add metadata 
metadata <- read.csv("ITS_lavaux.tsv", row.names = 2, sep = "\t", check.names = FALSE)
sample_data_phy <- sample_data(metadata)

# Merge into a Phyloseq Object
physeq <- phyloseq(otu_table_phy, tax_table_phy, sample_data_phy)

# drop empty cols 
physeq <- physeq %>%
  filter_taxa(function(x) sum(x) > 0, prune = TRUE)


############################################################     

### 2. HEATMAP 
cols <- setNames(viridis::viridis_pal()(length(unique(samdat_tbl(physeq)$Plot_ID))), 
                 unique(samdat_tbl(physeq)$Plot_ID))

# log2
physeq %>%
  tax_transform("log2", rank = "Genus", zero_replace = 1) %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_widtsh = 0.3, size = grid::unit(1, "cm"), ylim = 0:1)
    ),
    colors = heat_palette(palette='Purples', rev=TRUE),
    sample_anno = sampleAnnotation(
      Plot_ID = anno_sample("Plot_ID"),
      col = list(Plot_ID = cols), border = TRUE
    )
  )

# presence / absence 
physeq %>%
  tax_transform("binary", rank = "Genus") %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_widtsh = 0.3, size = grid::unit(1, "cm"), ylim = 0:1)
    ),
    colors = heat_palette(palette='Purples', rev=TRUE),
    sample_anno = sampleAnnotation(
      Plot_ID = anno_sample("Plot_ID"),
      col = list(Plot_ID = cols), border = TRUE
    )
  )

# Calculate prevalence at the Genus level
physeq_genus <- tax_glom(physeq, taxrank = "Genus")  # Aggregate at Genus level
otu_genus <- abundances(physeq_genus)  # Extract OTU table
# Convert to binary presence/absence
otu_binary <- (otu_genus > 0) * 1  
# Calculae prevalence (fraction of samples where each genus is present)
prev_df <- data.frame(
  Genus = tax_table(physeq_genus)[, "Genus"],
  Prevalence = rowSums(otu_binary) / ncol(otu_binary))
print(prev_df)

# Genus Prevalence 
# g__Saccharomyces 0.340694006
# g__Hanseniaspora 0.735015773
# g__Rhodotorula 0.003154574
# g__Torulaspora 0.009463722
# g__Candida 0.003154574


############################################################     

### 2. ANCOM-BC

# Define the factors to test
factors <- c("Altitude", "Average_slope", "Exposition", "Average_radiation")

# Initialize list to store results
results_list <- list()


res <- ancombc2(
  data = physeq, tax_level = 'Genus',
  fix_formula = "Altitude + Average_slope", rand_formula = "(1 | Year)",
  pseudo_sens = TRUE,  p_adj_method = "holm", 
  prv_cut = 0.10, s0_perc = 0.05,
  group = "Average_slope", struc_zero = FALSE, neg_lb = FALSE,
  alpha = 0.05, n_cl = 2, verbose = TRUE,
  global = FALSE, pairwise = TRUE, 
  dunnet = FALSE, trend = FALSE,
  iter_control = list(tol = 1e-5, max_iter = 20, 
                      verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
  trend_control = NULL)
)



# Combine results from all factors into a single dataframe
beta_results <- bind_rows(results_list, .id = "Factor") %>%
  pivot_wider(names_from = Factor, values_from = Beta)




Heatmap(
  beta_df %>% spread(Factor, Beta) %>% column_to_rownames("Taxa"),
  name = "Beta Coefficient",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
  show_row_names = TRUE,
  show_column_names = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method = "complete"
)

