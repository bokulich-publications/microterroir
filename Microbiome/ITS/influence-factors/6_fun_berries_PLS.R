# MixOmics PLS with all Berry samples (different years)
## MAYBE JUST DELETE?? 

### 0. SETUP 
library(mixOmics)
library(caret)
library(BiocParallel)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)

setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/DIABLO-data/Berries')

### 1. LOAD DATA 
# Read the CSV without using the first row as column names
df_raw <- read.csv("MFA_table.csv", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
df_raw <- df_raw[-3, ]  # Remove the 3rd row
df_raw <- df_raw[, -1]  # Remove the 1st column

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
groups <- c('Chemistry', 'Climate', 'Fungi', 'Metadata', 'Plots')
# Combine metadata data
data <- list(
  Chemistry = data_list[["Chemistry"]],
  Climate = data_list[["Climate"]], 
  Plots = data_list[["Plots"]] )
# Combine microbiome data
microbiome <- list(
  Fungi = data_list[["Fungi"]])
# Year / Plot 
metadata <- list(
  Metadata = data_list[['Metadata']])


# Function to convert a data frame to numeric matrix
convert_to_numeric <- function(df) {
  # Convert all columns to numeric
  numeric_df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
  # Convert to matrix
  return(as.matrix(numeric_df))}
# Convert metabolomics data
data <- list(
  Chemistry = convert_to_numeric(data_list[["Chemistry"]]),
  Climate = convert_to_numeric(data_list[["Climate"]]),
  Plots = convert_to_numeric(data_list[["Plots"]]))
# Convert microbiome data
microbiome <- list( Fungi = convert_to_numeric(data_list[["Fungi"]]))


# Combine the predictor variables (Chemistry, Climate, Plots) into one data frame
predictors <- cbind(data$Chemistry, data$Climate, data$Plots)
# response 
response <- microbiome$Fungi  
dim(response)

# remove samples with low variance (othwise it's an issue with scaling)
response_rmLowVar <- lapply(list(response), function(block) {
  nzv <- nearZeroVar(block)
  if (length(nzv) > 0) {
    return(block[, -nzv, drop = FALSE])
  } else {
    return(block)
  }
})[[1]]
dim(response_rmLowVar)

# Standardize data (important for PLS)
predictors <- scale(predictors)
response <- scale(response_rmLowVar)
sum(is.na(response))

predictors <- as.matrix(predictors)
response <- as.matrix(response)

#########################################################
# Perform PLS and cross-validation to find the optimal number of components
pls_cv <- perf(pls(X = predictors, Y = response), validation = "Mfold", folds = 5)





#########################################################

### 2. SET UP DIABLO MODEL for PLOT 

# Define outcome variable (Y)
Y <- factor(metadata$Metadata$Year)

# Design matrix (1 indicates connection between datasets)
design <- matrix(1, ncol = length(X), nrow = length(X), 
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

#########################################################

### 3. TEST NUMBER OF FEATURES TO USE
# Test what keepX values we should use! 
# note, Complex data structure!! there are some linear dependencies, plus the sample size is too small relative to the number of variables! 
set.seed(42)
# Adjust test.keepX ( #nb. reduce this otherwise it won't run)
test.keepX <- lapply(X, function(block) {
  c(1, min(5, ncol(block)), min(10, ncol(block)))})
# Run tune.block.splsda with adjusted parameters 
folds <- min(2, min(table(Y)))  # Ensures folds don’t exceed class sizes
tune.block.splsda(X, Y, 
                    ncomp = 2,
                    test.keepX = test.keepX,
                    design = design,
                    validation = "Mfold", 
                    folds = folds, 
                    nrepeat = 5,
                    near.zero.var = TRUE,
                    dist = "centroids.dist", 
                    measure = "BER", progressBar=TRUE)


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



### 5. VISUALIZATION 
# Plot sample plot
plotIndiv(diablo_tuned, group = Y, ind.names = FALSE, 
          legend = TRUE, title = 'DIABLO: samples for plot',  col.per.group = custom_colors_3)

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
png("diablo_group_representation_plot.png", width = 6, height = 5, units = "in", res = 1000)  # Adjust resolution
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
  theme(legend.position = "right")  # Keep legend for clarity
dev.off()

# Plot variable selection
png("diablo_variable_plot.png", width = 6, height = 6, units = "in", res = 1000)  # Adjust resolution
plotVar(diablo_tuned, var.names = TRUE,
        pch = list(rep(4, 152), rep(4, 164), rep(4, 92), rep(4, 19), rep(4, 19)),  
        cex = list(rep(2, 152), rep(2, 164), rep(2, 92), rep(2, 19), rep(2, 19)),  
        col = list(rep(block_colors["LCMS_Pos"], 152),  
                   rep(block_colors["LCMS_Neg"], 164),  
                   rep(block_colors["GCMS"], 92),  
                   rep(block_colors["Fungi"], 19),  
                   rep(block_colors["Bacteria"], 19)),  
        title = 'DIABLO Variable Plot')
dev.off()



# Perform cross-validation
# Get selected variables
selected_vars <- selectVar(diablo_tuned, block = 1:5)
selected_vars
diablo_perf <- perf(diablo_tuned, validation = "Mfold", folds = 3, nrepeat = 5)
# Plot performance
plot(diablo_perf)



# Create the Circos plot
png("diablo_circos_plot_corr90.png", width = 7, height = 7, units = "in", res = 1000)  
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

