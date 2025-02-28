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
                    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});","}")))


setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance')


### 1. LOAD DATA AND TURN TO PHYLOSEQ OBJECT 
# load csv (created on the cloud in 6_fun_investigateTaxa.ipynb)
# all plant pathogens
#otu_data <- read.csv("FunPathogens_crOTU90_annotated.csv", row.names = 1, check.names = FALSE) --> don't trust OTU
otu_data <- read.csv("berries_harvest_ASVs_rarefied_labled.tsv", sep = '\t', row.names = 1, check.names = FALSE)

# Create the taxonomy table by separating taxonomy strings
tax_table <- data.frame(Taxon = colnames(otu_data)) %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies"), 
           sep = ";", fill = "right")
# Ensure row names are unique before setting them
otu_names <- colnames(otu_data)
otu_names <- make.unique(otu_names)  # Make OTU names unique
rownames(tax_table) <- otu_names  # Assign unique names as row names

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

# some unidentified Genera 
physeq <- physeq %>% tax_fix(unknowns = c("s__unidentified", "g__unidentified", "f__unidentified"))


############################################################     
# Filter to grapevine specific genera

# Define the genera of interest with "g__" prefix
genera_of_interest <- c("Botrytis", "Phaeomoniella", "Phaeoacremonium", "Neofusicoccum", 
                        "Ilyonectria", "Lasiodiplodia", "Fomitiporia", "Diplodia", 
                        "Diaporthe", "Neoscytalidium", "Biscogniauxia", 
                        "Plasmopara", "Erysiphaceae", "Neonectria")
# Add "g__" prefix
genera_of_interest <- paste0("g__", genera_of_interest)
# Filter phyloseq object to keep only selected genera
physeq_GVpathogens <- subset_taxa(physeq, Genus %in% genera_of_interest)
physeq_GVpathogens


############################################################     

### 2. HEATMAP 
cols <- setNames(viridis::viridis_pal()(length(unique(samdat_tbl(physeq)$Year))), 
                 unique(samdat_tbl(physeq)$Year))

# log2
physeq_GVpathogens %>%
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
physeq_GVpathogens %>%
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
physeq_genus <- tax_glom(physeq, taxrank = "Genus")  # Aggregate at Genus level
otu_genus <- abundances(physeq_genus)  # Extract OTU table
# Convert to binary presence/absence
otu_binary <- (otu_genus > 0) * 1  
# Calculae prevalence (fraction of samples where each genus is present)
prev_df <- data.frame(
  Genus = tax_table(physeq_genus)[, "Genus"],
  Prevalence = rowSums(otu_binary) / ncol(otu_binary))
# Save to CSV
#write.csv(prev_df, "FunPathogens_Genus_Prevalence.csv", row.names = FALSE)
write.csv(prev_df, "FunGrapevinePathogens_Genus_Prevalence.csv", row.names = FALSE)

print(prev_df)
