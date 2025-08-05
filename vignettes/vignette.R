## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(collapse = TRUE, comment = "#>")
library(dominatR)
library(airway)
library(SummarizedExperiment)
library(ggplot2)
data("airway")
airway_se <- airway  

## ----eval=FALSE---------------------------------------------------------------
# # From GitHub (development version)
# if (!requireNamespace("remotes", quietly = TRUE))
#     install.packages("remotes")
# remotes::install_github("EthanCHEN6/dominatR_testing", force = TRUE)

## ----load-packages, message=FALSE, warning=FALSE------------------------------
# Load required packages
library(dominatR)
library(SummarizedExperiment)

# Load airway dataset
library(airway)
data(airway)
airway_se <- airway

## -----------------------------------------------------------------------------
airway_se
assayNames(airway_se)
dim(assay(airway_se))
head(assay(airway_se))

## -----------------------------------------------------------------------------
# Prepare input matrix
count_mat <- assay(airway)

# Apply min-max normalization
airway_minmax <- minmax_normalization(count_mat, new_min = 0, new_max = 1)

# Inspect structure
dim(airway_minmax)
summary(as.vector(airway_minmax))
head(airway_minmax[, 1:5])

## -----------------------------------------------------------------------------
df_scaled <- minmax_normalization(count_mat, new_min = 10, new_max = 20)
head(df_scaled)  # All columns now range from 10 to 20

## -----------------------------------------------------------------------------
se <- airway

# Option A: Overwrite the default assay
se1 <- minmax_normalization(se)
head(assay(se1))

# Option B: Write to a new assay slot
se2 <- minmax_normalization(se, new_assay_name = "minmax_counts")

## -----------------------------------------------------------------------------
## Apply quantile normalization
airway_quantile <- quantile_normalization(airway_se)

## Check result
dim(assay(airway_quantile))
summary(as.vector(assay(airway_quantile)))
head(assay(airway_quantile)[1:5, 1:5])

## -----------------------------------------------------------------------------
df <- assay(airway)
# Normalize without log2-transform
df_cpm <- cpm_normalization(df, log_trans = FALSE)
head(df_cpm[, 1:5])

# Normalize with log2-transform
df_cpm_log <- cpm_normalization(df, log_trans = TRUE)
head(df_cpm_log[, 1:5])

## -----------------------------------------------------------------------------
library(SummarizedExperiment)

# Apply in-place normalization (overwrite assay)
se1 <- cpm_normalization(airway, log_trans = FALSE)
head(assay(se1))

# Save to a new assay slot
se2 <- cpm_normalization(airway, log_trans = TRUE, new_assay_name = 
                            "cpm_logged")
head(assay(se2, "cpm_logged"))


## -----------------------------------------------------------------------------
new_counts <- matrix(sample(1:100000, nrow(airway) * ncol(airway), TRUE),
                    nrow = nrow(airway))
rownames(new_counts) <- rownames(airway)
colnames(new_counts) <- colnames(airway)

assay(airway, "new_raw") <- new_counts

se3 <- cpm_normalization(airway, assay_name = "new_raw", new_assay_name = 
                            "cpm_new_raw")
head(assay(se3, "cpm_new_raw"))

## -----------------------------------------------------------------------------
## Calculate gene length
rowData(airway_se)$gene_length <- rowData(airway_se)$gene_seq_end - 
    rowData(airway_se)$gene_seq_start

## Apply RPKM normalization
airway_se_rpkm <- rpkm_normalization(airway_se, gene_length, log_trans = TRUE)

## Check the result
dim(assay(airway_se_rpkm))  # Check the dimensions
summary(as.vector(assay(airway_se_rpkm)))  # Summary statistics for all values
head(assay(airway_se_rpkm)[1:5, 1:5])

## -----------------------------------------------------------------------------
## Calculate gene_length if not provided
rowData(airway_se)$gene_length <- rowData(airway_se)$gene_seq_end - 
    rowData(airway_se)$gene_seq_start

## Apply TPM normalization
airway_tpm <- tpm_normalization(airway_se, log_trans = TRUE)

## Check result
dim(assay(airway_tpm))
summary(as.vector(assay(airway_tpm)))
head(assay(airway_tpm)[1:5, 1:5])

## -----------------------------------------------------------------------------
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = "tpm_norm")
se <- se[1:1000, ]
#' Rename columns for consistency
colnames(se) <- paste('Column_', 1:8, sep ='')

## -----------------------------------------------------------------------------
plot_circle(
    x = se,
    n = 8,
    entropyrange     = c(0, 3),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
plot_circle(
    x = se,
    n = 8,
    entropyrange     = c(0, 1.5),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
plot_circle(
    x = se,
    n = 8,
    entropyrange     = c(2, 3),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
plot_circle(
    x = se,
    n = 8,
    column_variable_factor = 'gene_biotype',
    entropyrange     = c(2,3),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
plot_circle(
    x = se,
    n = 8,
    column_variable_factor = 'gene_biotype',
    point_size = 3,
    entropyrange     = c(0,1.5),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    assay_name = 'tpm_norm',
)

## -----------------------------------------------------------------------------
# Emphasize miRNA genes in orange
plot_circle(
    x = se,
    n = 8,
    column_variable_factor = 'gene_biotype',
    point_size = 3,
    entropyrange     = c(0,1.5),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    assay_name = 'tpm_norm',
    point_fill_colors = c('protein_coding' = 'orange'),
    point_line_colors = c('protein_coding' = 'orange')
)

## -----------------------------------------------------------------------------

se_result <- plot_circle(
    x = se,
    n = 8,
    column_variable_factor = 'gene_biotype',
    point_size = 3,
    entropyrange     = c(0,1.5),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = TRUE,
    assay_name = 'tpm_norm',
    point_fill_colors = c('protein_coding' = 'orange'),
    point_line_colors = c('protein_coding' = 'orange')
)

## -----------------------------------------------------------------------------
se_result[[1]]
head(se_result[[2]])

## -----------------------------------------------------------------------------
#' First we extract the normalized data as a data.frame:

df <- assay(se, 'tpm_norm') |> as.data.frame()
colnames(df) <- paste('Column_', 1:8, sep ='')

## -----------------------------------------------------------------------------
plot_circle(
    x = df,
    n = 8,
    entropyrange     = c(0, 3),
    magnituderange   = c(0, Inf),
    label  = 'legend', 
    output_table = FALSE
)

## -----------------------------------------------------------------------------
# Genes with entropy between 2-3 (more balanced expression)
plot_circle(
    x = df,
    n = 8,
    entropyrange     = c(2, 3),
    magnituderange   = c(0, Inf),
    label  = 'legend', 
    output_table = FALSE
)



## -----------------------------------------------------------------------------
#' Genes with entropy between 0-2 (more specialized expression)
plot_circle(
    x = df,
    n = 8,
    entropyrange     = c(0, 2),
    magnituderange   = c(0, Inf),
    label  = 'legend', 
    output_table = FALSE
)

## -----------------------------------------------------------------------------
plot_circle(
    x = df,
    n = 8,
    entropyrange     = c(0, 2),
    magnituderange   = c(0, Inf),
    label  = 'curve',
    output_table = FALSE
)

## -----------------------------------------------------------------------------
#' Emphasize expression dominance in Columns 1, 3, and 5
plot_circle(
    x = df,
    n = 8,
    entropyrange     = c(0, 2),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    background_alpha_polygon = 0.2,
    background_na_polygon = 'transparent',
    background_polygon = c('Column_1'  = 'indianred',
                        'Column_3' = 'lightblue',
                        'Column_5' = 'lightgreen'),
    point_fill_colors = c('Column_1'  = 'darkred',
                        'Column_3' = 'darkblue',
                        'Column_5' = 'darkgreen'),
    point_line_colors = c('Column_1'  = 'black',
                        'Column_3' = 'black',
                        'Column_5' = 'black')
)

# 1.2 Using factor variables

#' Add a factor column for grouping
set.seed(123)  # For reproducibility
df$factor <- sample(c('A', 'B', 'C', 'D'), size = nrow(df), replace = TRUE)

## -----------------------------------------------------------------------------
plot_circle(
    x = df,
    n = 8,
    column_variable_factor = 'factor',
    entropyrange     = c(0, 2),
    magnituderange   = c(0, Inf),
    label  = 'legend',
    output_table = FALSE,
    background_alpha_polygon = 0.2,
    background_na_polygon = 'transparent',
    background_polygon = c('Column_1'  = 'indianred',
                        'Column_3' = 'lightblue',
                        'Column_5' = 'lightgreen')
)

## -----------------------------------------------------------------------------
plot_circle(
    x = df,
    n = 8,
    column_variable_factor = 'factor',
    entropyrange     = c(0, 2),
    magnituderange   = c(0, Inf),
    label  = 'curve',
    output_table = FALSE,
    background_alpha_polygon = 0.02,
    background_na_polygon = 'transparent',
    point_fill_colors = c('A' = 'black',
                        'B' = 'gray',
                        'C' = 'white',
                        'D' = 'orange'),
    point_line_colors = c('A' = 'black',
                        'B' = 'gray',
                        'C' = 'white',
                        'D' = 'orange')
)

## -----------------------------------------------------------------------------
#' When `output_table = TRUE`, returns a list containing:
#' 1. ggplot object
#' 2. Data frame with entropy, magnitude, and dominance information
plot_result <- plot_circle(
    x = df,
    n = 8,
    point_size =  2,
    column_variable_factor = 'factor',
    entropyrange     = c(0, 2),
    magnituderange   = c(0, Inf),
    label  = 'curve',
    output_table = TRUE,
    background_alpha_polygon = 0.02,
    background_na_polygon = 'transparent',
    point_fill_colors = c('A' = 'black',
                        'B' = 'gray',
                        'C' = 'white',
                        'D' = 'orange'),
    point_line_colors = c('A' = 'black',
                        'B' = 'gray',
                        'C' = 'white',
                        'D' = 'orange')
)

# View plot
plot_result[[1]]

# View data
head(plot_result[[2]])

## -----------------------------------------------------------------------------
# Data preprocessing
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
se <- se[1:1000, ]


# Creating the circle plot data

# First we create the circle plot with output_table = TRUE to get 
# the data needed for the frequency plot. We'll use gene biotype as our 
# factor variable.

circle_data <- plot_circle(
    x = se,
    n = 8,
    column_variable_factor = 'gene_biotype',
    entropyrange = c(0, Inf),
    magnituderange = c(0, Inf),
    label = 'legend',
    output_table = TRUE,
    assay_name = 'tpm_norm'
)


## -----------------------------------------------------------------------------
freq_plot_default <- plot_circle_frequency(
    n = 8,
    circle = circle_data,
    single = TRUE,
    legend = TRUE,
    numb_columns = 1,
    filter_class = NULL,
    point_size = 2
)

# Display the plot
freq_plot_default[[1]]

# View aggregated data
head(freq_plot_default[[2]])

## -----------------------------------------------------------------------------
# Visualize each factor level in separate panels

plot_circle_frequency(
    n = 8,
    circle = circle_data,
    single = FALSE,
    legend = TRUE,
    numb_columns = 3,  # Arrange in 3 columns
    filter_class = NULL,
    point_size = 2
)

## -----------------------------------------------------------------------------
# Focus on specific gene biotypes

plot_circle_frequency(
    n = 8,
    circle = circle_data,
    single = FALSE,
    legend = TRUE,
    numb_columns = 1,  # Single column layout
    filter_class = c('protein_coding', 'snoRNA', 'miRNA'),
    point_size = 3  # Larger points for emphasis
)

## -----------------------------------------------------------------------------
# Create a combined plot showing only selected classes

plot_circle_frequency(
    n = 8,
    circle = circle_data,
    single = TRUE,
    legend = TRUE,
    numb_columns = 1,
    filter_class = c('protein_coding', 'miRNA', 'lincRNA'),
    point_size = 3
)

## -----------------------------------------------------------------------------
# Create data.frame version
df <- assay(se, 'tpm_norm') |> as.data.frame()
colnames(df) <- paste('Sample', 1:8, sep = '_')
df$gene_biotype <- rowData(se)$gene_biotype

# Create circle plot data
circle_df <- plot_circle(
    x = df,
    n = 8,
    column_variable_factor = 'gene_biotype',
    entropyrange = c(0, Inf),
    magnituderange = c(0, Inf),
    label = 'legend',
    output_table = TRUE
)

## -----------------------------------------------------------------------------
plot_circle_frequency(
    n = 8,
    circle = circle_df,
    single = FALSE,
    legend = TRUE,
    numb_columns = 2,
    filter_class = NULL,
    point_size = 1.5
)

## ----abacus-df-setup, message=FALSE, warning=FALSE----------------------------

se <- airway[1:1000, ]  
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')

# Prepare data frame
df_abacus <- as.data.frame(assay(se, "tpm_norm"))
df_abacus$gene_id <- rownames(df_abacus)
df_abacus <- df_abacus[, c("gene_id", setdiff(colnames(df_abacus), "gene_id"))]
head(df_abacus[, 1:5])

## -----------------------------------------------------------------------------
# Generate plot with minimal parameters
abacus_res <- plot_abacus(
    data              = df_abacus,               
    n                 = ncol(df_abacus) - 1,     
    x_variable        = "gene_id",               
    y_variables       = colnames(df_abacus)[-1], 
    percentiles       = 4,                       
    title             = "Gene Expression Dominance", 
    point_size        = 2,                      
    single            = TRUE                  
)

abacus_res[[1]]


## -----------------------------------------------------------------------------
## Data preparation
se <- airway[1:1000, ]  # Subset for faster computation
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
df <- as.data.frame(assay(se, 'tpm_norm'))
sample1 <- "SRR1039508"
sample2 <- "SRR1039516"

## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = se,
    column_name = c(sample1, sample2),
    col = c('lightgreen', 'indianred'),
    entropyrange = c(0, 0.1),
    maxvaluerange = c(4, 8),
    title = "SE Input: Low Entropy + Medium Expression"
)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = se,
    column_name = c(sample1, sample2),
    col = c('lightgreen', 'indianred'),
    entropyrange = c(0.1, 0.8),
    maxvaluerange = c(4, 8),
    title = "SE Input: Medium Entropy + Medium Expression"
)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = se,
    column_name = c(sample1, sample2),
    col = c('lightgreen', 'indianred'),
    entropyrange = c(0.8, 1),
    maxvaluerange = c(4, 8),
    title = "SE Input: High Entropy + Medium Expression"
)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = se,
    column_name = c(sample1, sample2),
    output_table = TRUE,
    col = c('lightgreen', 'indianred'),
    entropyrange = c(0.8, 1),
    maxvaluerange = c(4, 8)
)

str(res_rope)
head(res_rope)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    title = "Default Rope Plot"
)



## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    col = c('darkgreen', 'darkred'),
    title = "Custom Colors"
)


## -----------------------------------------------------------------------------
head(res_rope)

## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    col = c('darkgreen', 'darkred'),
    entropyrange = c(0, 0.1),
    title = "Low Entropy Genes (0-0.1)"
)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    col = c('darkgreen', 'darkred'),
    entropyrange = c(0, 0.1),
    maxvaluerange = c(2, Inf),
    title = "Low Entropy + High Expression"
)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    col = c('darkgreen', 'darkred'),
    entropyrange = c(0, 0.1),
    maxvaluerange = c(4, 8),
    title = "Low Entropy + Medium Expression"
)
head(res_rope[[2]])

## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    col = c('darkgreen', 'darkred'),
    entropyrange = c(0.1, 0.8),
    maxvaluerange = c(4, 8),
    title = "Medium Entropy + Medium Expression"
)


## -----------------------------------------------------------------------------
res_rope = plot_rope(
    x = df,
    column_name = c(sample1, sample2),
    col = c('darkgreen', 'darkred'),
    entropyrange = c(0.8, 1),
    maxvaluerange = c(4, 8),
    title = "High Entropy + Medium Expression"
)

## -----------------------------------------------------------------------------
## Minimal data preparation
se <- airway[1:1000, ]  # Subset for faster computation
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
se <- tpm_normalization(se, log_trans = TRUE, new_assay_name = 'tpm_norm')
df <- as.data.frame(assay(se, 'tpm_norm'))
samples <- c("SRR1039508", "SRR1039512", "SRR1039516")

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples
)


## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue')
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(0, 0.4)
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(0.4, 1.3)
)


## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(1.3, Inf)
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(1.2, Inf),
    maxvaluerange = c(2, Inf)
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(1.2, Inf),
    maxvaluerange = c(5, Inf)
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(1.2, Inf),
    maxvaluerange = c(10, Inf)
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = df,
    column_name = samples,
    col = c('indianred', 'lightgreen', 'lightblue'),
    entropyrange = c(1.2, Inf),
    maxvaluerange = c(2, Inf),
    plotAll = FALSE
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = se,
    column_name = samples,
    col = c('darkred', 'darkgreen', 'darkblue'),
    entropyrange = c(0, 0.4),
    maxvaluerange = c(0.1, Inf),
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = se,
    column_name = samples,
    col = c('darkred', 'darkgreen', 'darkblue'),
    entropyrange = c(0.4, 1.3),
    maxvaluerange = c(0.1, Inf),
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
res_rope = plot_triangle(
    x = se,
    column_name = samples,
    col = c('darkred', 'darkgreen', 'darkblue'),
    entropyrange = c(1.3, Inf),
    maxvaluerange = c(0.1, Inf),
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
triangle_data <- plot_triangle(
    x = se,
    column_name = samples,
    output_table = TRUE,
    entropyrange = c(1.3, Inf),
    maxvaluerange = c(0.1, Inf),
    assay_name = 'tpm_norm'
)

## -----------------------------------------------------------------------------
# View first 6 rows of the output data
head(triangle_data)

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

