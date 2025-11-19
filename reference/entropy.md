# Compute Shannon Entropy on row-normalized data

Compute Shannon Entropy on row-normalized data

## Usage

``` r
entropy(x, assay_name = NULL, new_assay_name = "Entropy")
```

## Arguments

- x:

  A data.frame (with numeric columns) or a SummarizedExperiment (with an
  assay of numeric data).

- assay_name:

  (SummarizedExperiment only) The name of the assay to transform and
  compute Entropy on. If NULL, uses the first assay.

- new_assay_name:

  If you prefer to store Q-values in a \*new\* assay, provide a name. By
  default 'Entropy'

## Value

- If `x` is a data.frame: returns the same data.frame in which numeric
  columns have been replaced by their row-wise proportions, and an
  `Entropy` column is appended.

- If `x` is a SummarizedExperiment: returns the same
  SummarizedExperiment in with a new assay (Default name is `Entropy`)
  and `rowData(x)$Entropy` is added.

## Examples

``` r
library(SummarizedExperiment)
library(airway)
data('airway')

se = airway

# Only use a random subset of 1000 rows
set.seed(123)
idx <- sample(seq_len(nrow(se)), size = min(1000, nrow(se)))
se <- se[idx, ]

# -------------------------------
# 1) Using a data.frame
# -------------------------------
df = assay(se) |> as.data.frame()
df = entropy(df)

## The function adds a new column called Entropy and transform all
## the counts accordingly
head(df)
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166  0.0000000  0.0000000  0.0000000 0.00000000  0.0000000
#> ENSG00000266931  0.0000000  0.0000000  0.0000000 0.00000000  0.0000000
#> ENSG00000104774  0.1271726  0.1250082  0.1504558 0.07903194  0.1227127
#> ENSG00000267583  0.0000000  1.0000000  0.0000000 0.00000000  0.0000000
#> ENSG00000227581  0.3333333  0.0000000  0.0000000 0.00000000  0.0000000
#> ENSG00000227317  0.0000000  0.0000000  0.0000000 0.00000000  0.0000000
#>                 SRR1039517 SRR1039520 SRR1039521   Entropy
#> ENSG00000260166  0.0000000  1.0000000  0.0000000 0.0000000
#> ENSG00000266931  0.0000000  0.0000000  0.0000000 0.0000000
#> ENSG00000104774  0.1482915  0.1190398  0.1282875 2.9791667
#> ENSG00000267583  0.0000000  0.0000000  0.0000000 0.0000000
#> ENSG00000227581  0.6666667  0.0000000  0.0000000 0.9182958
#> ENSG00000227317  0.0000000  0.0000000  0.0000000 0.0000000

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

## The function adds a new assay called 'Entropy' with the transformed
## counts.
## This name can be modified with the 'new_assay_name' parameter
## In the rowData dataframe a new column called Entropy is added.
se2 <- entropy(se, new_assay_name = 'Entropy')
se2
#> class: RangedSummarizedExperiment 
#> dim: 1000 8 
#> metadata(1): ''
#> assays(2): counts Entropy
#> rownames(1000): ENSG00000260166 ENSG00000266931 ... ENSG00000160886
#>   ENSG00000142871
#> rowData names(11): gene_id gene_name ... symbol Entropy
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample

## In case the experiment has multiple assays, the function allows you to
## choose which assay to use.
new_matrix =  matrix(data = sample(x = seq(1, 100000),
                                   size = nrow(se) * ncol(se),
                                   replace = TRUE),
                     nrow = nrow(se),
                     ncol = ncol(se))
rownames(new_matrix) = rownames(se)
colnames(new_matrix) = colnames(se)

## Creating a new assay called new counts
assay(se, 'new_counts') = new_matrix


## Saving the entropy values as Entropy_newmatrix using the assay 'new
## counts'
se2 = entropy(se,
              new_assay_name = 'Entropy_newmatrix',
              assay_name = 'new_counts')

se2
#> class: RangedSummarizedExperiment 
#> dim: 1000 8 
#> metadata(1): ''
#> assays(3): counts new_counts Entropy_newmatrix
#> rownames(1000): ENSG00000260166 ENSG00000266931 ... ENSG00000160886
#>   ENSG00000142871
#> rowData names(11): gene_id gene_name ... symbol Entropy
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
```
