# Heat Tree of fungal berry communities 

library(metacoder)
library(dplyr)
library(ggplot2)
library(factoextra)  # For clustering visualization
library(cluster)     # For clustering functions
library(tidyverse)
library(taxa)

set.seed(42) 

setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/HeatTree')
otu_data_raw <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/berries_harvest21_cOTUs_rarefied_labled.tsv", sep = '\t', row.names = 1, check.names = FALSE,encoding = "UTF-8")

# 1. Assign unique OTU IDs with "OTU_" prefix
taxa_names <- colnames(otu_data_raw)
tax_data <- tibble(Taxa = taxa_names, OTU_ID = paste0("OTU_", 1:length(taxa_names)))

# 2. Rename columns in otu_data_raw with OTU_IDs
otu_data_renamed <- otu_data_raw
colnames(otu_data_renamed) <- tax_data$OTU_ID

# 3. Transpose the data
otu_data <- otu_data_renamed %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -Sample, names_to = "OTU_ID", values_to = "Count") %>%
  pivot_wider(id_cols = OTU_ID, names_from = Sample, values_from = Count) %>%
  select(OTU_ID, everything())

# 4. Create the tax_data with original Taxa names and OTU_ID, and taxonomy column
tax_data <- tax_data %>%
  mutate(OTU_ID = as.character(OTU_ID)) %>%
  mutate(Taxa = str_replace_all(Taxa, ";;", ";Unassigned;")) %>%
  mutate(Taxa = str_replace_all(Taxa, ";$", ";Unassigned")) %>%
  mutate(Taxa = str_replace_all(Taxa, "^$", "Unassigned")) %>%
  separate(Taxa,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Confidence"),
           sep = ";",
           fill = "right") %>%
  mutate(Confidence = 1) %>%
  select(OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)
tax_data$taxonomy <- paste(tax_data$Kingdom, 
                           tax_data$Phylum, 
                           tax_data$Class, 
                           tax_data$Order, 
                           tax_data$Family, 
                           tax_data$Genus, 
                           tax_data$Species, 
                           sep = ";")


# 2. Sample Metadata 
metadata_climate <- read.csv("/Users/lfloerl/Desktop/MICROTERROIR/Data/DifferentialAbundance/Metadata/ITS_Lavaux_Climate.tsv", row.names = 2, sep = "\t", check.names = FALSE)
hmp_samples_df <- metadata_climate %>%
  filter(sample_type == "must", Year == 2021)

# to get from continous metadata to categorical we cluster!
# Filter the data to remove rows where either median_rh or median_temperature is NA
hmp_samples_filtered <- hmp_samples_df %>%
  filter(!is.na(median_rh) & !is.na(median_temperature))

# Determine the optimal number of clusters for median_rh using the elbow method
median_rh_df <- as.data.frame(hmp_samples_filtered$median_rh)
colnames(median_rh_df) <- "median_rh"
fviz_nbclust(median_rh_df, kmeans, method = "wss")

# Apply k-means clustering for median_rh with optimal k (set manually after checking the plot)
optimal_k_rh <- 2  # Change this based on the elbow method result
kmeans_result_rh <- kmeans(median_rh_df, centers = optimal_k_rh)

# Assign cluster labels ensuring the higher RH cluster is "high"
cluster_means <- aggregate(median_rh_df$median_rh, list(cluster = kmeans_result_rh$cluster), mean)
high_cluster <- cluster_means$cluster[which.max(cluster_means$x)]

hmp_samples_filtered <- hmp_samples_filtered %>%
  mutate(rh_category = ifelse(kmeans_result_rh$cluster == high_cluster, "high", "low"))

# Convert to factor for plotting
hmp_samples_filtered$rh_category <- factor(hmp_samples_filtered$rh_category, levels = c("low", "high"))

# Scatterplot for median_rh with labeled clusters
ggplot(hmp_samples_filtered, aes(x = rownames(hmp_samples_filtered), y = median_rh, 
                                 color = rh_category, shape = as.factor(Year))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("low" = "blue", "high" = "red")) +
  labs(x = "Sample ID", y = "Median RH", color = "Category", shape = "Year") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



# Determine the optimal number of clusters for TEMPERATURE using the elbow method
median_temp_df <- as.data.frame(hmp_samples_filtered$median_temperature)
colnames(median_temp_df) <- "median_temperature"
fviz_nbclust(median_temp_df, kmeans, method = "wss")

# Apply k-means clustering for median_temperature with k = 3
optimal_k_temp <- 3  
kmeans_result_temp <- kmeans(hmp_samples_filtered$median_temperature, centers = optimal_k_temp)

# Calculate the mean value of each cluster
cluster_means_temp <- tapply(hmp_samples_filtered$median_temperature, kmeans_result_temp$cluster, mean)

# Sort clusters by their mean values
sorted_clusters_temp <- order(cluster_means_temp)

# Assign category labels based on the sorted clusters
hmp_samples_filtered <- hmp_samples_filtered %>%
  mutate(temp_category = case_when(
    kmeans_result_temp$cluster == sorted_clusters_temp[1] ~ "low",
    kmeans_result_temp$cluster == sorted_clusters_temp[2] ~ "mid",
    kmeans_result_temp$cluster == sorted_clusters_temp[3] ~ "high"
  ))

# Convert to factor for plotting with all three levels
hmp_samples_filtered$temp_category <- factor(hmp_samples_filtered$temp_category, levels = c("low", "mid", "high"))

# Scatterplot for median_temperature with labeled clusters
ggplot(hmp_samples_filtered, aes(x = rownames(hmp_samples_filtered), y = median_temperature, 
                                 color = temp_category, shape = as.factor(Year))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("low" = "blue", "mid" = "green", "high" = "red")) +
  labs(x = "Sample ID", y = "Median Temp", color = "Category", shape = "Year") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


# make to a tibble 
hmp_samples_filtered <- hmp_samples_filtered %>%
  rownames_to_column(var = "SampleID") %>%
  rename(
    Name = Plot_ID,
    Type = sample_type,
    Year = Year,
    Plot = Plot,
    Altitude = Altitude,
    Average_slope = Average_slope,
    Exposition = Exposition,
    Average_radiation = Average_radiation,
    Soil_thickness = Soil_thickness,
    Soil_type = Soil_type,
    Soil_depth = Soil_depth,
    Hydromorphie = Hydromorphie,
    Hydromorphie_code = Hydromorphie_code,
    Geology = Geology,
    Quadrant = Quadrant,
    Cluster = Cluster,
    Plot_PCA_kMeans_Cluster = Plot_PCA_kMeans_Cluster,
    Unnamed = `Unnamed: 0`,
    average_rh = average_rh,
    median_rh = median_rh,
    maximum_rh = maximum_rh,
    minimum_rh = minimum_rh,
    cv_rh = cv_rh,
    GDD = GDD,
    average_temperature = average_temperature,
    median_temperature = median_temperature,
    maximum_temperature = maximum_temperature,
    minimum_temperature = minimum_temperature,
    accumulated_temperature = accumulated_temperature,
    cv_temperature = cv_temperature,
    rh_category = rh_category,
    temp_category = temp_category
  ) %>%
  select(SampleID, Name, Year, Plot, Altitude, Average_slope, Exposition, Average_radiation, average_rh, median_rh, median_temperature, id, temp_category, rh_category)


# Identify the common SampleID between hmp_samples_filtered and the column names of otu_data
common_samples <- intersect(hmp_samples_filtered$SampleID, colnames(otu_data))
# Subset hmp_samples_filtered to include only the rows with common SampleID
hmp_samples <- hmp_samples_filtered %>%
  filter(SampleID %in% common_samples)

# Subset otu_data to include only the columns corresponding to the common SampleID, but keep the OTU_ID column
otu_data_subset <- otu_data %>%
  select(OTU_ID, all_of(common_samples))


# 5. Join the data frames
otu_data_joined <- left_join(otu_data_subset, tax_data, by = "OTU_ID")

# make object 
obj <- parse_tax_data(otu_data_joined,
                      class_cols = "taxonomy",
                      class_sep = ";",
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))


########################################################################
# Plotting! 

# BASIC PLOT 
# no to: 
#     layout = "davidson-harel", initial_layout = "davidson-harel" 
#     layout = "automatic"
# circular:  layout = "reingold-tilford", initial_layout = "reingold-tilford"


obj %>% 
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", initial_layout = "reingold-tilford")



########################################################################
# COMPARE RH 

# comparing heat trees between RH 
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols =  hmp_samples$SampleID,
                                       groups = hmp_samples$rh_category)

set.seed(2)
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = low, 
            initial_layout = "re", layout = "da")
set.seed(2)
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE) %>% # subset to the class rank
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = high, 
            initial_layout = "re", layout = "da")

# Differential heat trees
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = hmp_samples$SampleID)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols =  hmp_samples$SampleID,
                                      groups = hmp_samples$rh_category)
# wilcox test (corrected)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
# ONLY show sign. enriched features
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0


# genus level! 
set.seed(1)
obj %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  metacoder::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree(node_label = cleaned_names,
            node_size = n_obs, # number of OTUs
            node_color = log2_median_ratio, # difference between groups
            node_color_interval = c(-10, 10), # symmetric interval
            node_color_range = c("#39558CFF", "gray", "#B8DE29FF"), 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median counts",
            layout = "da", initial_layout = "re", # 
            output_file = "RH_21_differential_heat_tree.png")


########################################################################
# COMPARE TEMP 

# Differential heat trees
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
                                       cols = hmp_samples$SampleID)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols =  hmp_samples$SampleID,
                                      groups = hmp_samples$temp_category)
# wilcox test (corrected)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
# ONLY show sign. enriched features
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0

# genus level! 
set.seed(1)
obj %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  metacoder::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree(node_label = cleaned_names,
            node_size = n_obs, # number of OTUs
            node_color = log2_median_ratio, # difference between groups
            node_color_interval = c(-10, 10), # symmetric interval
            node_color_range = c("#39558CFF", "gray", "#B8DE29FF"), 
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median counts",
            layout = "da", initial_layout = "re", 
            title = "high vs low Temp",
            output_file = "Temp_differential_heat_tree.png")


set.seed(1)
obj %>%
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree_matrix(obj,
                   data = "diff_table",
                   node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                   node_label = taxon_names,
                   node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                   node_color_range = diverging_palette(), # The built-in palette for diverging data
                   node_color_trans = "linear", # The default is scaled by circle area
                   node_size_axis_label = "Number of OTUs",
                   node_color_axis_label = "Log2 ratio median proportions",
                   layout = "re", initial_layout = "re") 



