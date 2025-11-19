# Min-Max Normalization

Scales each column of a matrix (or SummarizedExperiment assay) so that
the minimum value in that column is mapped to `new_min` and the maximum
value is mapped to `new_max`

## Usage

``` r
minmax_normalization(
  x,
  new_min = 0,
  new_max = 1,
  assay_name = NULL,
  new_assay_name = NULL
)
```

## Arguments

- x:

  A numeric `matrix`, `data.frame`, or `SummarizedExperiment`.

- new_min:

  The lower bound of the new range (default 0).

- new_max:

  The upper bound of the new range (default 1).

- assay_name:

  If `x` is a SummarizedExperiment, name of the assay to normalize.
  Defaults to the first assay if none is specified.

- new_assay_name:

  If `x` is a SummarizedExperiment, name of a new assay to store the
  normalized data. If `NULL`, overwrites the assay specified by
  `assay_name`.

## Value

- If `x` is a data.frame or matrix, returns a matrix of column-wise
  scaled values (same dimensions as `x`).

- If `x` is a SummarizedExperiment, returns the same
  SummarizedExperiment object with the chosen or new assay replaced by
  the scaled values.

## Examples

``` r
library(SummarizedExperiment)
library(airway)
data('airway')

se <- airway

# Only use a random subset of 1000 rows
set.seed(123)
idx <- sample(seq_len(nrow(se)), size = min(1000, nrow(se)))
se <- se[idx, ]

# -------------------------------
# 1) Using a data.frame
# -------------------------------

df <- assay(se)

df1 <- minmax_normalization(df)

apply(df1, 2, range)
#>      SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
#> [1,]          0          0          0          0          0          0
#> [2,]          1          1          1          1          1          1
#>      SRR1039520 SRR1039521
#> [1,]          0          0
#> [2,]          1          1

## Using a new range
df1 <- minmax_normalization(df, new_min = 5, new_max = 10)

apply(df1, 2, range)
#>      SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
#> [1,]          5          5          5          5          5          5
#> [2,]         10         10         10         10         10         10
#>      SRR1039520 SRR1039521
#> [1,]          5          5
#> [2,]         10         10

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

# If now new_assay_name is provided, then overwrites existing assay
se2 <- minmax_normalization(se)

apply(assay(se2), 2, range)
#>      SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
#> [1,]          0          0          0          0          0          0
#> [2,]          1          1          1          1          1          1
#>      SRR1039520 SRR1039521
#> [1,]          0          0
#> [2,]          1          1


# If new new_assay_name, normalization stored in a new object
se2 <- minmax_normalization(se, new_assay_name = 'minmax_counts')

apply(assay(se2, 'minmax_counts'), 2, range)
#>      SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
#> [1,]          0          0          0          0          0          0
#> [2,]          1          1          1          1          1          1
#>      SRR1039520 SRR1039521
#> [1,]          0          0
#> [2,]          1          1

# A specific assay can also be selected
new_matrix <-  matrix(data = sample(x = seq(1, 100000),
                                   size = nrow(se) * ncol(se),
                                   replace = TRUE),
                     nrow = nrow(se),
                     ncol = ncol(se))
rownames(new_matrix) <- rownames(se)
colnames(new_matrix) <- colnames(se)

## Creating a new assay called new counts
assay(se, 'new_counts') <- new_matrix

se2 <- minmax_normalization(se,
                           new_assay_name = 'minmax_counts_new',
                           assay_name = 'new_counts')

apply(assay(se2, 'minmax_counts_new'), 2, range)
#>      SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
#> [1,]          0          0          0          0          0          0
#> [2,]          1          1          1          1          1          1
#>      SRR1039520 SRR1039521
#> [1,]          0          0
#> [2,]          1          1

## Using a different range
se2 <- minmax_normalization(se,
                           new_assay_name = 'minmax_counts_new',
                           assay_name = 'new_counts',
                           new_min = 10,
                           new_max = 20)

apply(assay(se2, 'minmax_counts_new'), 2, range)
#>      SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
#> [1,]         10         10         10         10         10         10
#> [2,]         20         20         20         20         20         20
#>      SRR1039520 SRR1039521
#> [1,]         10         10
#> [2,]         20         20
```
