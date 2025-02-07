# MixOmics Diablo with all PostMLF Microbiome and LC-MS data 


### 0. SETUP 

#BiocManager::install("mixOmics")

library(mixOmics)

setwd('/Users/lfloerl/Desktop/MICROTERROIR/Data/Metabolomics/Multiomics_preparedData/PostMLF_MicrobiomeMetabolome/')

# Define custom colors 
custom_colors_3 <- c("#440154FF", "#238A8DFF", "#FDE725FF")  


### 1. LOAD DATA 

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
# Assign the group information back
group_info <- data.frame(Variable = variables, Group = groups)
# Split the data into separate data frames by group
data_list <- split.default(df_clean, groups)

# Split the data into separate data frames for each group
groups <- c("Metabolites_Pos", "Metabolites_Neg", "Wine Chemistry", "Plots", "Climate", "Fungi", "Bacteria", "Metadata")

# Combine metabolomics data
metabolomics <- list(
  LCMS_Pos = data_list[["Metabolites_Pos"]],
  LCMS_Neg = data_list[["Metabolites_Neg"]])

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
  LCMS_Pos = convert_to_numeric(data_list[["Metabolites_Pos"]]),
  LCMS_Neg = convert_to_numeric(data_list[["Metabolites_Neg"]]))

# Convert microbiome data
microbiome <- list( Fungi = convert_to_numeric(data_list[["Fungi"]]),
  Bacteria = convert_to_numeric(data_list[["Bacteria"]]))

# Combine all blocks into X
X <- c(metabolomics, microbiome)
sapply(X, ncol)
# LCMS_Pos LCMS_Neg    Fungi Bacteria 
# 191      170       88       91 
# identify and remove features with zero or near-zero variance
X <- lapply(X, function(block) {
  nzv <- nearZeroVar(block)
  if(length(nzv) > 0) {return(block[, -nzv, drop = FALSE])
  } else {return(block) }})
# how many left? 
sapply(X, ncol)
# LCMS_Pos LCMS_Neg    Fungi Bacteria 
# 191      170       33       21 
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


### 3. TEST NUMBER OF FEATURES TO USE
# Test what keepX values we should use! 
# note, Complex data structure!! there are some linear dependencies, plus the sample size is too small relative to the number of variables! 
set.seed(42)
# Adjust test.keepX
test.keepX <- lapply(X, function(block) {
  n <- ncol(block)
  c(5, 10, 20, 30, 50, min(100, n))})
# Run tune.block.splsda with adjusted parameters
tune_result <- tune.block.splsda(X, Y, 
                                 ncomp = 3,  # Reduced number of components
                                 test.keepX = test.keepX,
                                 design = design,
                                 validation = "Mfold", 
                                 folds = 5, 
                                 nrepeat = 10,  # Increased number of repeats
                                 dist = "max.dist", 
                                 measure = "BER",
                                 near.zero.var = TRUE,
                                 max.iter = 200)  # Increased max iterations
# how much does it suggest we use? 
print(tune_result$choice.keepX)
# LCMS_Pos: 3  3 10
# LCMS_Neg: 15  3  5
# Fungi: 3 3 3
# Bacteria: 3 3 3


### 4. RUN DIABLO
# tuned model
diablo_tuned <- block.splsda(X = X, Y = Y, design = design,  near.zero.var = TRUE, 
                             ncomp = 3,  # Using 3 components as in the tuning
                             keepX = list(LCMS_Pos = c(3, 3, 10),
                                          LCMS_Neg = c(15, 3, 5),
                                          Fungi = c(3, 3, 3),
                                          Bacteria = c(3, 3, 3)))

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
          title = '',
          ellipse = TRUE,
          style = "graphics",
          pch = 16)
# Add a single legend to the side
legend("topright", 
       legend = levels(Y), 
       col = custom_colors_3, 
       pch = 16, 
       title = "Plots",
       inset = c(-0.8, 0),
       bty = "n")
# Reset the plotting parameters
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)


# Plot variable selection
plotVar(diablo_tuned, var.names = TRUE, 
        pch = list(rep(4, 191), rep(4, 170), rep(4, 33), rep(4, 21)),
        cex = list(rep(2, 191), rep(2, 170), rep(2, 33), rep(2, 21)),
        title = 'DIABLO Variable Plot')


# Define the color palette
color_palette <- c("#440154FF", "#238A8DFF", "#B8DE29FF", "#39558CFF")
# Set up the plotting area
par(mar = c(5, 4, 4, 8), xpd = TRUE)
# Create the plot
plotVar(diablo_tuned, var.names = TRUE, 
        pch = list(rep(16, 191), rep(16, 170), rep(16, 33), rep(16, 21)),
        cex = list(rep(1.5, 191), rep(1.5, 170), rep(1.5, 33), rep(1.5, 21)),
        col = list(rep(color_palette[1], 191), 
                   rep(color_palette[2], 170), 
                   rep(color_palette[3], 33), 
                   rep(color_palette[4], 21)),
        title = 'DIABLO Variable Plot',
        legend = FALSE)
# Add a custom legend
legend("topright", 
       legend = c("LCMS_Pos", "LCMS_Neg", "Fungi", "Bacteria"),
       col = color_palette,
       pch = 16,
       pt.cex = 1.5,
       title = "Data Types",
       inset = c(-0.2, 0),
       bty = "n")

# Reset the plotting parameters
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE)


# Get selected variables
selected_vars <- selectVar(diablo_tuned, block = 1:5)
selected_vars

# Perform cross-validation
diablo_perf <- perf(diablo_tuned, validation = "Mfold", folds = 3, nrepeat = 5)
# Plot performance
plot(diablo_perf)


# Create the Circos plot
circosPlot(diablo_tuned, 
           cutoff = 0.7,  # Correlation cutoff for displaying links
           line = TRUE,   # Show lines between correlated variables
           size.variables = 0.5,  # Size of variable names
           size.labels = 0.7,     # Size of block labels
           color.blocks = c("#440154FF", "#238A8DFF", "#B8DE29FF", "#39558CFF"),  # Colors for each block
           color.cor = c("red", "blue"),  # Colors for positive and negative correlations
           var.adj = 1.2,  # Adjust variable label positions
           title = "DIABLO Circos Plot")




#########################################################
# OLD TESTING 
## larger model 
keepX <- list(
  LCMS_Pos = c(20, 30, 40),
  LCMS_Neg = c(20, 30, 40),
  Fungi = c(10, 15, 20),
  Bacteria = c(5, 10, 15))
# Ensure we're not selecting more features than available
min_features <- sapply(X, ncol)
keepX <- lapply(names(keepX), function(name) {
  pmin(keepX[[name]], min_features[name])})
names(keepX) <- names(X)
# downstream there are issues with the variance in the bacteria block 
nzv_bacteria <- caret::nearZeroVar(X$Bacteria, saveMetrics = TRUE)
zero_var_features <- rownames(nzv_bacteria)[nzv_bacteria$zeroVar]
X$Bacteria <- X$Bacteria[, !colnames(X$Bacteria) %in% zero_var_features]
keepX$Bacteria <- pmin(keepX$Bacteria, ncol(X$Bacteria))

diablo_tuned_large <- block.splsda(X = X, Y = Y, design = design, 
                                   ncomp = 3, near.zero.var = TRUE,
                                   keepX = keepX)


## compare the model performance - small sample size!! 
perf_small <- perf(diablo_tuned_small, validation = "Mfold", folds = 3, nrepeat = 5)
perf_large <- perf(diablo_tuned_large, validation = "Mfold", folds = 3, nrepeat = 5)
# Plot performance:
plot(perf_small)  
plot(perf_large)  

# different validation method
perf_large <- perf(diablo_tuned_large, validation = "loo")
perf_large$MER