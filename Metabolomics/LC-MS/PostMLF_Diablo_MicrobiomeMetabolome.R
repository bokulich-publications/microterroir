# MixOmics Diablo with all PostMLF Microbiome and LC-MS, GC-MS data 


### 0. SETUP 

#BiocManager::install("mixOmics")

library(mixOmics)
library(caret)
library(BiocParallel)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)
library(reshape2)
library(tibble)


setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/Metabolomics/Multiomics_preparedData/PostMLF_MicrobiomeMetabolome/')

# Define custom colors 
custom_colors_3 <- c("#440154FF", "#238A8DFF", "#FDE725FF")  

block_colors <- c(
  'LCMS_Pos' = '#20A386FF',
  'LCMS_Neg' = '#238A8DFF',
  'GCMS' = '#32648EFF',
  'Fungi' = '#453781FF',
  'Bacteria' = '#440154FF'  )

 
### 1. LOAD DATA 

# Read the CSV without using the first row as column names
df_raw <- read.csv("MFA_table.csv", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
df_raw[2, ][df_raw[2, ] == "NP-001596.1"] <- "Hexadecanedioic_acid"

# Extract the group names from the first row
groups <- as.character(df_raw[1, ])
# Extract the variable names from the second row
variables <- as.character(df_raw[2, ])
# Remove the first two rows (which are now stored separately)
df_clean <- df_raw[-c(1, 2), ]
# Assign proper column names
colnames(df_clean) <- variables
# Assign the group information back
group_info <- data.frame(Variable = variables, Group = groups)
# Split the data into separate data frames by group
data_list <- split.default(df_clean, groups)

# Split the data into separate data frames for each group
groups <- c("LCMS_Pos", "LCMS_Neg","GCMS", "Wine Chemistry", "Plots", "Climate", "Fungi", "Bacteria", "Metadata")

# Combine metabolomics data
metabolomics <- list(
  LCMS_Pos = data_list[["LCMS_Pos"]],
  LCMS_Neg = data_list[["LCMS_Neg"]], 
  GCMS = data_list[["GCMS"]] )

# Combine microbiome data
microbiome <- list(
  Fungi = data_list[["Fungi"]],
  Bacteria = data_list[["Bacteria"]])

# Combine metadata
metadata <- list(
  Wine_Chemistry = data_list[["Wine Chemistry"]],
  Plots = data_list[["Plots"]],
  Climate = data_list[["Climate"]],
  Metadata = data_list[['Metadata']])


# Function to convert a data frame to numeric matrix
convert_to_numeric <- function(df) {
  # Convert all columns to numeric
  numeric_df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
  # Convert to matrix
  return(as.matrix(numeric_df))}

# Convert metabolomics data
metabolomics <- list(
  LCMS_Pos = convert_to_numeric(data_list[["LCMS_Pos"]]),
  LCMS_Neg = convert_to_numeric(data_list[["LCMS_Neg"]]),
  GCMS = convert_to_numeric(data_list[["GCMS"]]))

# Convert microbiome data
microbiome <- list( Fungi = convert_to_numeric(data_list[["Fungi"]]),
  Bacteria = convert_to_numeric(data_list[["Bacteria"]]))

# Combine all blocks into X
X <- c(metabolomics, microbiome)
sapply(X, ncol)
# LCMS_Pos LCMS_Neg     GCMS    Fungi Bacteria 
# 191      170      151       88       92 
# identify and remove features with zero or near-zero variance
# note: the LC-MS dataframes have very little variance so we just keep them 
X <- lapply(X, function(block) {
  nzv <- nearZeroVar(block)
  if(length(nzv) > 0) {return(block[, -nzv, drop = FALSE])
  } else {return(block) }})
# also remove columns with high colinearity 
X <- lapply(X, function(block) {
  correlations <- cor(block)
  highCorr <- findCorrelation(correlations, cutoff = 0.95) # Relaxed cutoff
  block[, -highCorr]
})
# how many left?
sapply(X, ncol)
# LCMS_Pos LCMS_Neg     GCMS    Fungi Bacteria 
# 152      164       92       19       19 
# must be a matrix 
X <- lapply(X, as.matrix)

#########################################################

### 2. SET UP DIABLO MODEL for YEAR
# note, for plot this does not work because the groups are not balanced and very small (1-3 per)

# Define outcome variable (Y) : let's try 'year' 
Y <- factor(metadata$Metadata$year)

# Design matrix (1 indicates connection between datasets)
design <- matrix(1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

#########################################################

### 3. TEST NUMBER OF FEATURES TO USE
# Test what keepX values we should use! 
# note, Complex data structure!! there are some linear dependencies, plus the sample size is too small relative to the number of variables! 
set.seed(42)
# Define range of components to test
ncomp_range <- 1:5
# Adjust test.keepX ( #nb. reduce this otherwise it won't run)
test.keepX <- lapply(X, function(block) {
  c(1, min(5, ncol(block)), min(10, ncol(block)), min(20, ncol(block)))})
# Run tune.block.splsda with adjusted parameters 
tune_results <- lapply(ncomp_range, function(ncomp) {
  tune.block.splsda(X, Y, 
                    ncomp = ncomp,
                    test.keepX = test.keepX,
                    design = design,
                    validation = "Mfold", 
                    folds = 5, 
                    nrepeat = 5,  # Reduce repeats for testing
                    near.zero.var = TRUE,
                    dist = "centroids.dist", 
                    measure = "BER", progressBar=TRUE)
})

# Extract optimal values for each ncomp
optimal_keepX <- lapply(tune_results, function(res) res$choice.keepX)
ber_values <- sapply(tune_results, function(res) min(res$error.rate)) # Extract BER
optimal_ncomp <- ncomp_range[which.min(ber_values)] # Select best ncomp

# Print suggested number of components and keepX
list(optimal_ncomp = optimal_ncomp, optimal_keepX = optimal_keepX[[which.min(ber_values)]])

# Optimal Number of Components: 2

# Block	      Component 1	  Component 2
# LCMS_Pos	  20	          1
# LCMS_Neg	  1	            5
# GCMS	      20	          1
# Fungi	      10	          1
# Bacteria	   5          	1

#########################################################

### 4. RUN DIABLO
# tuned model
diablo_tuned <- block.splsda(X = X, Y = Y, design = design, near.zero.var = TRUE, 
                             ncomp = 2, 
                             keepX = list(LCMS_Pos = c(20,1),
                                          LCMS_Neg = c(1,5),
                                          GCMS = c(20,1),
                                          Fungi = c(10,1),
                                          Bacteria = c(5,1)))


# Classification accuracy
# Perform cross-validation to assess classification accuracy
set.seed(123)  # Ensure reproducibility
perf_diablo <- perf(diablo_tuned, validation = "Mfold", folds = 9, nrepeat = 50, progressBar = TRUE)

# Print performance metrics
print(perf_diablo$error.rate)

# Compute accuracy for each data block
accuracy_list <- lapply(perf_diablo$error.rate, function(error_matrix) {
  1 - error_matrix  # Convert error rate to accuracy
})
print(accuracy_list)

# Average accuracy across blocks and components
mean_accuracy <- sapply(accuracy_list, function(acc_matrix) {
  mean(acc_matrix, na.rm = TRUE)
})

# Print mean accuracy per block
print(mean_accuracy)



# Extract proportion of explained variance
explained_var <- diablo_tuned$prop_expl_var
print(explained_var)



### 5. VISUALIZATION 
# Plot sample plot
png("/Users/lfloerl/Desktop/MICROTERROIR/Figures/diablo_PLSscores_plot.png", width = 6, height = 7.5, units = "in", res = 1000)  # Adjust resolution
plotIndiv(diablo_tuned, group = Y, ind.names = FALSE, ellipse = TRUE,
          legend = TRUE,  col.per.group = custom_colors_3)
dev.off()


# Set up the plotting area with extra space for the legend
par(mar = c(5, 4, 4, 8), xpd = TRUE)
# Create the plot
plotIndiv(diablo_tuned, 
          comp = c(1, 2),
          ind.names = FALSE, 
          legend = FALSE, 
          col.per.group = custom_colors_3,
          title = '',  ellipse = TRUE,
          style = "graphics",
          pch = 16)
# Add a single legend to the side
legend("topright", 
       legend = levels(Y),  col = custom_colors_3, 
       pch = 16, title = "Year",
       inset = c(-0.8, 0), bty = "n")
# Reset the plotting parameters
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)


### GROUP REPRESENTATION 
# Extract sample scores (latent variables) for each block
block_scores <- lapply(names(X), function(block) {
  data.frame(diablo_tuned$variates[[block]], Block = block)  # Add block name
})
# Combine all blocks into one dataframe
df_scores <- bind_rows(block_scores)
# Compute group means (centroids) per block
group_means <- df_scores %>%
  group_by(Block) %>%
  summarise(across(starts_with("comp"), mean))

# Plot using ggplot2 with custom colors
png("/Users/lfloerl/Desktop/MICROTERROIR/Figures/diablo_group_representation_plot.png", width = 4, height = 4, units = "in", res = 1000)  # Adjust resolution
ggplot(group_means, aes(x = comp1, y = comp2, color = Block)) +
  geom_point(size = 5) +   # Plot group centroids
  geom_text_repel(aes(label = Block), box.padding = 0.5, max.overlaps = Inf) +  # Prevent overlap
  theme_minimal() +
  labs(title = "",
       x = "Component 1",
       y = "Component 2") +
  scale_x_continuous(labels = scientific_format(digits = 2)) +  # Keep scientific notation
  scale_y_continuous(labels = scientific_format(digits = 2)) +  # Keep scientific notation
  scale_color_manual(values = block_colors) +  # Apply custom colors
  theme(legend.position = "none")  # Remove legend
dev.off()


# Plot variable selection
plotVar(diablo_tuned, var.names = TRUE,
        pch = list(rep(4, 152), rep(4, 164), rep(4, 92), rep(4, 19), rep(4, 19)),  
        cex = list(rep(2, 152), rep(2, 164), rep(2, 92), rep(2, 19), rep(2, 19)),  
        col = list(rep(block_colors["LCMS_Pos"], 152),  
                   rep(block_colors["LCMS_Neg"], 164),  
                   rep(block_colors["GCMS"], 92),  
                   rep(block_colors["Fungi"], 19),  
                   rep(block_colors["Bacteria"], 19)),  
        title = 'DIABLO Variable Plot')


# make a custom Variable Plot (from Jorge!)
plotVar(diablo_tuned, comp = c(1,2), title = "PLS-DA Variable Loadings", style = "ggplot2", 
            var.names = TRUE, 
            pch = list(rep(4, 152), rep(4, 164), rep(4, 92), rep(4, 19), rep(4, 19)),  
            cex = list(rep(2, 152), rep(2, 164), rep(2, 92), rep(2, 19), rep(2, 19))) 

# Export the data 
plsda_loadings <- plotVar(diablo_tuned, comp = c(1,2), title = "PLS-DA Variable Loadings", style = "ggplot2", 
                                  var.names = FALSE, 
                                  pch = list(rep(4, 152), rep(4, 164), rep(4, 92), rep(4, 19), rep(4, 19)),  
                                  cex = list(rep(2, 152), rep(2, 164), rep(2, 92), rep(2, 19), rep(2, 19))) 


# get the top 10 loadings
top10_vars <- plsda_loadings %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(Contribution = abs(x) + abs(y)) %>%
  arrange(desc(Contribution)) %>%
  slice_head(n = 10)
print(top10_vars$Variable)

# rename them 
top10_vars_named <- top10_vars %>%
  mutate(NewName = case_when(
    Variable == "Hanseniaspora" ~ "Hanseniaspora sp.",
    Variable == "Farnesol...2Z.6Z..." ~ "Farnesol, (2Z,6Z)-",
    Variable == "X5.Hydroxytryptophol.1" ~ "Hydroxytryptophol",
    Variable == "Cer.NS.d34.3" ~ "C16-0(Palmitoyl)ceramide",
    Variable == "Lavandulyl..tetrahydro.." ~ "Tetrahydro lavandulyl acetate",
    Variable == "Norepinephrine" ~ "Norepinephrine",
    Variable == "Azelaic.acid" ~ "Azelaic acid",
    Variable == "Saccharomyces_cerevisiae" ~ "Saccharomyces cerevisiae",
    Variable == "Phe.Tyr" ~ "Phenylalanyltyrosine",
    Variable == "L.Saccharopine.1" ~ "L-Saccharopine",
    TRUE ~ Variable
  ))

# Join renamed labels back to full loading data
loadings_df <- plsda_loadings %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(Contribution = abs(x) + abs(y)) %>%
  left_join(top10_vars_named %>% select(Variable, NewName), by = "Variable") %>%
  mutate(Label = if_else(!is.na(NewName), NewName, ""),
         Highlight = if_else(!is.na(NewName), "Top 10", "Others"))


# Manual Plot
png("/Users/lfloerl/Desktop/MICROTERROIR/Figures/diablo_variable_plot_manual.png", width = 7, height = 7, units = "in", res = 1000)  # Adjust resolution
ggplot(loadings_df, aes(x = x, y = y, color = Highlight, fill = Highlight)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  annotate("path", x = 0.5 * cos(seq(0, 2*pi, length.out = 100)),
           y = 0.5 * sin(seq(0, 2*pi, length.out = 100)), color = "darkgray", linetype = "solid") +
  annotate("path", x = 1.0 * cos(seq(0, 2*pi, length.out = 100)),
           y = 1.0 * sin(seq(0, 2*pi, length.out = 100)), color = "black", linetype = "solid") +
  
  geom_point(aes(size = Highlight), shape = 21, stroke = 1) +
  geom_text_repel(aes(label = Label), max.overlaps = 10, size = 4, box.padding = 0.4) +
  
  theme_bw(base_size = 12) +
  theme(legend.position = "none") + 
  xlab("Component 1") +
  ylab("Component 2") +
  ggtitle("") +
  
  scale_fill_manual(values = c("Top 10" = "#FFA500", "Others" = "gray60")) +
  scale_colour_manual(values = c("Top 10" = "#FFA500", "Others" = "gray60")) +
  scale_size_manual(values = c("Top 10" = 3, "Others" = 1.5)) +
  coord_fixed()
dev.off()



# Perform cross-validation
# Get selected variables
selected_vars <- selectVar(diablo_tuned, block = 1:5)
selected_vars
diablo_perf <- perf(diablo_tuned, validation = "Mfold", folds = 3, nrepeat = 5)
# Plot performance
plot(diablo_perf)



# CIRCOS PLOT 

# Extract all variable names across all blocks
all_diablo_vars <- lapply(diablo_tuned$loadings, rownames) %>%
  unlist(use.names = FALSE) %>%
  unique()
print(all_diablo_vars)

# Define your renaming dictionary as a named vector
rename_dict <- c(
  "Hanseniaspora" = "Hanseniaspora sp.",
  "Farnesol...2Z.6Z..." = "Farnesol, (2Z,6Z)-",
  "X5.Hydroxytryptophol.1" = "5-Hydroxytryptophol",
  "Cer.NS.d34.3" = "C16-0(Palmitoyl)ceramide",
  "Lavandulyl..tetrahydro.." = "Tetrahydro lavandulyl acetate",
  "Norepinephrine" = "Norepinephrine",
  "Azelaic.acid" = "Azelaic acid",
  "Saccharomyces_cerevisiae" = "Saccharomyces cerevisiae",
  "Phe.Tyr" = "Phenylalanyltyrosine",
  "L.Saccharopine.1" = "L-Saccharopine",
  "X.g__S085.1" = "Dehalococcoidia sp.",
  "X.s__Acetobacter_pasteurianus" = "Acetobacter pasteurianus",
  "X.g__Aquabacterium" = "Aquabacterium sp.",
  "X.g__KD4.96.1" = "Chloroflexi sp.", 
  "X.s__Azospirillum_brasilense" = "Azospirillum brasilense",
  "X.s__Oenococcus_oeni" = "Oenococcus oeni",
  "Endoconidioma_populi" = "Endoconidioma populi",
  "Ascomycota.1" = "Ascomycota sp.",
  "Epicoccum_nigrum" = "Epicoccum nigrum", 
  "Cladosporium_austrohemisphaericum" = "Cladosporium austrohemisphaericum", 
  "Alternaria_subcucurbitae" = "Alternaria subcucurbitae",
  "Alternaria.1" = "Alternaria sp.",
  "Stemphylium_solani" = "Stemphylium solani",
  "X4.Methyl.2..hydroxy.bicyclo.4.3.0.non.3.en.6.methanol" = "3-Me-BCN",
  "X.Nonanoic.acid." = "Nonanoic acid",
  "X1.Octene..6.methyl." = "1-Octene-6-methyl",
  "X2.4.Nonadienol" = "2,4-Nonadienol",
  "D.Hexanal" = "Hexanal",
  "Benzaldehyde..2.4.dimethyl.." = "Benzaldehyde, 2,4-dimethyl-",
  "Ethyl.linoleate" = "Ethyl linoleate",
  "Hexanoic.acid" = "Hexanoic acid",
  "Diethyl.succinate" = "Diethyl succinate", 
  "X1.Octanol" = "1-Octanol",
  "X1.Hexanol..3.5.5.trimethyl." = "3,5,5-Trimethyl-1-hexanol",
  "X1.Hexanol..2.ethyl." = "2-Ethyl-1-hexanol",                                                                                                                                                                            
  "X1.Hexene..3.5.dimethyl." = "3,5-Dimethyl-1-hexene",
  "X2.Hexenoic.acid..ethyl.ester" = "Ethyl hexanoate", 
  "X3.Hexen.1.ol...Z.." = "3-Hexen-1-ol, (Z)-",
  "X4..benzoyloxy..2H.pyran.3.one" = "BzO", 
  "Ethyl.3.hydroxytridecanoate" = "Ethyl 3-hydroxytridecanoate",
  "X1.Heptanol" = "1-Heptanol", 
  "Hexadecanoic.acid..ethyl.ester" = "Hexadecanoic acid, ethyl ester",
  "UDP.N.acetylglucosamine" = "UDP-GlcNAc",
  "Digalacturonic.acid" = "Digalacturonic acid",
  "Suberic.acid" = "Suberic acid",
  "trans.Aconitic.acid.1" = "trans-Aconitic acid",
  "Hexadecanedioic_acid" = "Hexadecanedioic acid", 
  "X2.2.difluoro.N..2.methyl.6..trifluoromethyl.pyridin.3.yl..2..phenylthio.acetamide" = "TFPMA",
  "N..5.acetamidopentyl.acetamide" = "N-(5-acetamidopentyl)acetamide",
  "Argininosuccinic.acid.1" = "Argininosuccinic acid",  
  "Argininosuccinic.acid" = "Argininosuccinic acid",
  "L.....Methionine.1" = "Methionine",
  "X3.phenyl.5..1.2.3.thiadiazol.4.yl..1.2.4.oxadiazole" = "3-Phenyl-1,2,4-oxadiazole-5-thiol",
  "α.Aspartylphenylalanine.1" = "α.Aspartylphenylalanine",
  "X6.bromo.2..3.pyridylmethyl..2.3.dihydro.1H.benzo.de.isoquinoline.1.3.dione" = "6-Br-2-Et-BIQ",
  "N..2..1.5.dimethyl.4.nitro.1H.pyrazol.3.yl.vinyl..N.N.dimethylamine" = "1,5-dimethyl-4-nitro-1H-pyrazole",
  "X4.Acetamidobutanoic.acid" = "4-Acetamidobutanoic acid",
  "Proline.Hydroxyproline" = "Hydroxyprolyl-Proline",
  "N.Methyldioctylamine" = "N-Methyldioctylamine",
  "Mycosphaerella_tassiana"  = "Mycosphaerella tassiana"  
)

# Function to apply variable renaming across the DIABLO model
rename_variables_in_diablo <- function(diablo_model, rename_dict) {
  # Helper function to make names unique
  make_names_unique <- function(names_vec) {
    make.unique(names_vec, sep = ".")
  }
  
  # Rename and make rownames of loadings unique
  diablo_model$loadings <- lapply(diablo_model$loadings, function(block) {
    new_names <- ifelse(rownames(block) %in% names(rename_dict),
                        rename_dict[rownames(block)],
                        rownames(block))
    rownames(block) <- make_names_unique(new_names)
    block
  })
  
  # Rename and make column names of X unique
  diablo_model$X <- lapply(diablo_model$X, function(block) {
    new_names <- ifelse(colnames(block) %in% names(rename_dict),
                        rename_dict[colnames(block)],
                        colnames(block))
    colnames(block) <- make_names_unique(new_names)
    block
  })
  
  # Rename and make names$colnames unique
  diablo_model$names$colnames <- lapply(diablo_model$names$colnames, function(names_vec) {
    new_names <- ifelse(names_vec %in% names(rename_dict),
                        rename_dict[names_vec],
                        names_vec)
    make_names_unique(new_names)
  })
  
  return(diablo_model)
}


# Apply the renaming function
diablo_tuned_renamed <- rename_variables_in_diablo(diablo_tuned, rename_dict)

# Plot the Circos plot
png("/Users/lfloerl/Desktop/MICROTERROIR/Figures/diablo_circos_plot_corr90.png", width = 6.6, height = 6.5, units = "in", res = 1000)  
circosPlot(
  diablo_tuned_renamed, 
  cutoff = 0.90,
  line = TRUE,
  size.variables = 0.45,
  size.labels = 0.8,
  color.blocks = block_colors,
  color.cor = c("red", "blue"),
  var.adj = -0.5,
  title = "DIABLO Circos Plot"
)
dev.off()



# Create the Circos plot
#png("/Users/lfloerl/Desktop/MICROTERROIR/Figures/diablo_circos_plot_corr90.png", width = 6.5, height = 6.5, units = "in", res = 1000)  
circosPlot(diablo_tuned, 
           cutoff = 0.90,  # Correlation cutoff for displaying links
           line = TRUE,   # Show lines between correlated variables
           size.variables = 0.45,  # Size of variable names
           size.labels = 0.8,     # Size of block labels
           color.blocks = block_colors,  # Colors for each block
           color.cor = c("red", "blue"),  # Colors for positive and negative correlations
           var.adj = -0.5,  # Adjust variable label positions
           title = "DIABLO Circos Plot")
dev.off()

png("diablo_circos_plot_corr85.png", width = 7, height = 7, units = "in", res = 1000)  # Adjust resolution
circosPlot(diablo_tuned, 
           cutoff = 0.85,  # Correlation cutoff for displaying links
           line = TRUE,   # Show lines between correlated variables
           size.variables = 0.45,  # Size of variable names
           size.labels = 0.8,     # Size of block labels
           color.blocks = block_colors,  # Colors for each block
           color.cor = c("red", "blue"),  # Colors for positive and negative correlations
           var.adj = -0.5,  # Adjust variable label positions
           title = "DIABLO Circos Plot")
dev.off()



# Extract the similarity matrix
circos_plot <- circosPlot(diablo_tuned_renamed, 
                          cutoff = 0.90,
                          line = FALSE,
                          size.variables = 0.0001,
                          size.labels = 0.0001,
                          color.blocks = block_colors,
                          color.cor = c("red", "blue"),
                          var.adj = -0.5,
                          legend = FALSE,
                          var.names = NULL,
                          title = "")
similarity_matrix <- circos_plot

similarity_matrix_df <- melt(similarity_matrix, varnames = c("Row", "Col"), value.name = "Similarity")

# Ensure unique pairs (A-B == B-A)
similarity_matrix_df$Pair <- apply(similarity_matrix_df[, c("Row", "Col")], 1, function(x) paste(sort(x), collapse = "_"))

# Remove self-comparisons and filter for absolute similarity > 0.9
similarity_matrix_df_filtered <- similarity_matrix_df[similarity_matrix_df$Row != similarity_matrix_df$Col & 
                                                        (abs(similarity_matrix_df$Similarity) > 0.9), ]

# Remove duplicate pairs
similarity_matrix_df_filtered <- similarity_matrix_df_filtered[!duplicated(similarity_matrix_df_filtered$Pair), ]

# View the result
print(similarity_matrix_df_filtered)
# Extracting unique values from the 'Row' and 'Col' columns
unique_rows_cols <- unique(c(as.character(similarity_matrix_df_filtered$Row), as.character(similarity_matrix_df_filtered$Col)))
writeLines(unique_rows_cols, "Top_correlations_circos_unique.txt")


write.csv(similarity_matrix_df_filtered, "Top_correlations_circos.csv", row.names = FALSE)
