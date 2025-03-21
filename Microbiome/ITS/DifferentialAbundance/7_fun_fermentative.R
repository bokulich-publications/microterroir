# ABUNDANCE FERMENTATIVE ORGANISMS ON BERRIES 

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
library(gridExtra)
library(grid)
library(tidyr)
library(purrr)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});","}")))


setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/FermentativeYeasts')


### 1. LOAD DATA AND TURN TO PHYLOSEQ OBJECT 
# load csv (created on the cloud in 6_fun_investigateTaxa.ipynb)
otu_data <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/berries_harvest_cOTUs_rarefied_labled.tsv", sep = '\t', row.names = 1, check.names = FALSE,encoding = "UTF-8")
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
p <- physeq_yeasts %>%
  tax_transform("log2", rank = "Genus", zero_replace = 1) %>%
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_widtsh = 1, size = grid::unit(0.8, "cm"), ylim = 0:1)
    ),
    colors = heat_palette(palette='Purples', rev=TRUE),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 8),
    sample_anno = sampleAnnotation(
      Year = anno_sample("Year"),
      col = list(Year = cols), border = TRUE)
    )
png("/Users/lfloerl/Desktop/MICROTERROIR/Figures/FermYeast_Heatmap_Log2.png", width = 1600, height = 600, res = 300)
draw(p)  # Explicitly draw the heatmap
dev.off()

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



# TEST statistical singnificant presence of yeast between years! 
# Extract metadata
year_factor <- physeq_yeasts@sam_data[["Year"]]
# Create presence/absence matrix with years
otu_year <- cbind(t(otu_binary), Year = year_factor)
otu_year <- data.frame(otu_year)

# 1. Drop genera with a total prevalence of < 10%
otu_filtered <- otu_year %>%
  select(-Year) %>% # Temporarily remove 'Year' for calculations
  mutate_all(~ ifelse(. > 0, 1, 0)) %>% # Convert to presence/absence (1/0)
  select(where(~ mean(.) >= 0.1)) %>% # Select genera with prevalence >= 10%
  cbind(otu_year$Year) %>% # Re-add 'Year'
  rename(Year = `otu_year$Year`)

# 2. Print the occurrence of each genus per year
occurrence_per_year <- otu_filtered %>%
  pivot_longer(cols = -Year, names_to = "Genus", values_to = "Presence") %>%
  group_by(Year, Genus) %>%
  summarise(Occurrence = sum(Presence), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = Occurrence, values_fill = 0)

print(occurrence_per_year)

# 3. Calculate statistical significance using Fisher's exact test
years <- unique(otu_filtered$Year)
genera <- names(otu_filtered)[names(otu_filtered) != "Year"]

fisher_results <- map_dfr(genera, function(genus) {
  map_dfr(combn(years, 2, simplify = FALSE), function(year_pair) {
    year1 <- year_pair[1]
    year2 <- year_pair[2]
    
    data_subset <- otu_filtered %>%
      filter(Year %in% c(year1, year2)) %>%
      select(Year, !!genus)
    
    contingency_table <- table(data_subset)
    
    fisher_test <- fisher.test(contingency_table)
    
    data.frame(
      Genus = genus,
      Year1 = year1,
      Year2 = year2,
      p_value = fisher_test$p.value,
      OR = fisher_test$estimate
    )
  })
})

print(fisher_results)

#Optional: adjust p values for multiple comparisons using bonferroni or FDR
fisher_results_adj <- fisher_results %>%
  group_by(Genus) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni")) %>%
  ungroup()
 
print(fisher_results_adj)

#Genus                                                                                                                             Year1 Year2 p_value    OR  p_adj
#  1 k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycodaceae.g__Hanseniaspora.s__unidentified              2021  2022 0.00396 0.371 0.0119
#2 k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycodaceae.g__Hanseniaspora.s__unidentified              2021  2023 0.137   0.581 0.206 
#3 k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycodaceae.g__Hanseniaspora.s__unidentified              2022  2023 0.306   1.56  0.306 
#4 k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetaceae.g__Saccharomyces.s__Saccharomyces_cerevisiae  2021  2022 1       0.982 1     
#5 k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetaceae.g__Saccharomyces.s__Saccharomyces_cerevisiae  2021  2023 0.738   0.867 1     
#6 k__Fungi.p__Ascomycota.c__Saccharomycetes.o__Saccharomycetales.f__Saccharomycetaceae.g__Saccharomyces.s__Saccharomyces_cerevisiae  2022  2023 0.831   0.883 1 



# Assuming physeq_yeasts is your phyloseq object
p_log2 <- physeq_yeasts %>%
  tax_transform("log2", rank = "Genus", zero_replace = 1) %>%
  otu_table() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Genus", values_to = "log2_abundance") %>%
  left_join(sample_data(physeq_yeasts) %>% data.frame() %>% rownames_to_column("Sample"), by = "Sample") %>%
  select(Sample, Genus, log2_abundance, Year) # Assuming "Year" is in your sample data


# Wilcoxon Rank-Sum (Pairwise Comparisons)
wilcoxon_results <- p_log2 %>%
  group_by(Genus) %>%
  do({
    year_pairs <- combn(unique(.$Year), 2, simplify = FALSE)
    map_dfr(year_pairs, function(pair) {
      year1 <- pair[1]
      year2 <- pair[2]
      data_pair <- filter(., Year %in% c(year1, year2))
      test_result <- wilcox.test(log2_abundance ~ Year, data = data_pair)
      data.frame(Year1 = year1, Year2 = year2, p.value = test_result$p.value, statistic = test_result$statistic)
    })
  })

print(wilcoxon_results)

# Multiple comparisons adjustment
wilcoxon_results_adj <- wilcoxon_results %>%
  group_by(Genus) %>%
  mutate(p_adj = p.adjust(p.value, method = "bonferroni")) %>%
  ungroup()

print(wilcoxon_results_adj)

#Genus            Year1 Year2   p.value statistic    p_adj
#  1 g__Hanseniaspora  2021  2022   0.00389     6847    0.0117
#2 g__Hanseniaspora  2021  2023   0.852       5321    0.852 
#3 g__Hanseniaspora  2022  2023   0.0359       918    0.0538
#4 g__Rhodotorula    2021  2022 NaN           5450  NaN     
#5 g__Rhodotorula    2021  2023 NaN           5232  NaN     
#6 g__Rhodotorula    2022  2023 NaN           1200  NaN     
#7 g__Saccharomyces  2021  2022   0.546       5216.   0.819 
#8 g__Saccharomyces  2021  2023   0.382       4902.   0.819 
#9 g__Saccharomyces  2022  2023   0.867       1180.   0.867 
#10 g__Torulaspora    2021  2022   0.502       5500    0.502 
#11 g__Torulaspora    2021  2023   0.0937      5061    0.226 
#12 g__Torulaspora    2022  2023   0.151       1150    0.226 




## LOG FOLD CHANGE PER GENUS AND YEAR 

# Assuming physeq_yeasts is your phyloseq object
p_log2 <- physeq_yeasts %>%
  tax_transform("log2", rank = "Genus", zero_replace = 1) %>%
  otu_table() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Genus", values_to = "log2_abundance") %>%
  left_join(sample_data(physeq_yeasts) %>% data.frame() %>% rownames_to_column("Sample"), by = "Sample") %>%
  select(Sample, Genus, log2_abundance, Year)

# Create a function to calculate log2 changes between year pairs
calculate_log2_changes <- function(data, genus_name, year1, year2) {
  data_year1 <- data %>%
    filter(Genus == genus_name, Year == year1) %>%
    pull(log2_abundance)
  
  data_year2 <- data %>%
    filter(Genus == genus_name, Year == year2) %>%
    pull(log2_abundance)
  
  # Ensure both years have data for the genus
  if (length(data_year1) > 0 && length(data_year2) > 0) {
    log2_changes <- outer(data_year2, data_year1, "-") %>% as.vector()
    return(data.frame(Genus = genus_name, Year_pair = paste(year1, "-", year2), Log2_change = log2_changes))
  } else {
    return(data.frame()) # Return an empty data frame if data is missing
  }
}

# Get unique genera and years
genera <- unique(p_log2$Genus)
years <- unique(p_log2$Year)

# Calculate log2 changes for all year pairs and genera
log2_changes_df <- map_dfr(genera, function(gen) {
  combn(years, 2, simplify = FALSE) %>%
    map_dfr(function(year_pair) {
      calculate_log2_changes(p_log2, gen, year_pair[1], year_pair[2])
    })
})

# Create histograms
ggplot(log2_changes_df, aes(x = Log2_change)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_grid(Year_pair ~ Genus) +
  labs(title = "Log2 Changes in Genus Abundance Between Years",
       x = "Log2 Change",
       y = "Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 10))



############################################################     
 
### 2. ANCOM-BC MUTLIVARIATE 

# remove samples with missing data in any metadata column:
physeq_yeast_clean <- subset_samples(physeq_yeasts, complete.cases(sample_data(physeq)))
dim(physeq_yeasts@otu_table)
# 316   7
dim(physeq_yeast_clean@otu_table)
# 301   7


output <- ancombc2(
  data = physeq_yeast_clean,
  tax_level = 'Genus',
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


plot_lfc_heatmap <- function(res_df, file_name = NULL) {
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
  lighten_color <- function(color, factor = 0.3) {
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
  
  # Save or display heatmap
  if (!is.null(file_name)) {
    png(file_name, width = 1800, height = 3000, res = 600)
  }
  
  pheatmap(
    lfc_matrix, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    display_numbers = annotation_matrix,
    fontsize_number = 8,    # Smaller annotation text
    fontsize = 10,          # Smaller general font
    cellwidth = 30,         # Adjust cell width
    cellheight = 18,        # Adjust cell height
    color = color_palette, 
    breaks = breaks,
    legend_breaks = legend_breaks,  
    legend_labels = legend_labels,
    legend_title = "LFC")
  
  if (!is.null(file_name)) {
    dev.off() 
  }
}

plot_lfc_heatmap(output$res)  

plot_lfc_heatmap(output$res, file_name = "/Users/lfloerl/Desktop/MICROTERROIR/Figures/FermYeast_ANCOM_multivariat.png") 
dev.off()




##############

# Year as random effect 

# note the year needs to be a factor 
sample_data(physeq_yeast_clean)$Year <- as.factor(sample_data(physeq_yeast_clean)$Year)

output_year <- ancombc2(
  data = physeq_yeast_clean, 
  tax_level = 'Genus',
  fix_formula = "Altitude + Average_slope + Average_radiation + median_rh + median_temperature + Tartrate_gL + Malate_gL + Glucose_gL + Fructose_gL", 
  rand_formula = "(1 | Year)",  # Year as a random effect
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

plot_lfc_heatmap(output_year$res)  

plot_lfc_heatmap(output_year$res, file_name = "/Users/lfloerl/Desktop/MICROTERROIR/Figures/FermYeast_ANCOM_multivariat_RandomYear.png") 
dev.off()



