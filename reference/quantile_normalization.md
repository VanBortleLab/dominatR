# Quantile Normalization

Normalizes read counts by the quantile normalization method:

1.  Each sample (column) is sorted, and values at each rank are averaged
    across columns

2.  Each sample's values are replaced with the average of their
    respective rank

3.  If `log_trans = TRUE`, applies `log2(QN + 1)` transformation

## Usage

``` r
quantile_normalization(
  x,
  log_trans = FALSE,
  assay_name = NULL,
  new_assay_name = NULL
)
```

## Arguments

- x:

  A numeric `matrix` or `data.frame` of gene counts, or a
  `SummarizedExperiment` containing such counts.

  If a `SummarizedExperiment`,

  :   the function applies normalization to the specified assay (via
      `assay_name`).

  If a `data.frame`/`matrix`,

  :   the normalization is applied directly.

- log_trans:

  Logical. If `TRUE`, apply `log2(... + 1)` transform to the
  quantile-normalized values.

- assay_name:

  If `x` is a SummarizedExperiment, name of the assay to normalize.
  Defaults to the first assay if not specified.

- new_assay_name:

  If `x` is a SummarizedExperiment, name of a new assay in which to
  store the quantile-normalized (or log2-transformed) values. If `NULL`,
  overwrites the original assay.

## Value

A numeric **matrix** of quantile-normalized (or log2-normalized) values
if `x` is a data.frame or matrix. If `x` is a SummarizedExperiment,
returns the modified SummarizedExperiment with the normalized data
placed in the existing or new assay.

## Details

If `x` is a `SummarizedExperiment`, the function will extract the assay
using `assay_name`, apply quantile normalization, and return a new or
updated assay. If `x` is a matrix or data.frame, normalization is
applied directly to the input matrix.

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

## Without log transformation
df_qn <- quantile_normalization(df, log_trans = FALSE)
df_qn[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166      0.000      0.000      0.000        0.0       0.00
#> ENSG00000266931      0.000      0.000      0.000        0.0       0.00
#> ENSG00000104774   2164.875   2319.375   1946.875     1971.5    1610.75
#> ENSG00000267583      0.000      0.750      0.000        0.0       0.00
#> ENSG00000227581      0.625      0.000      0.000        0.0       0.00

## With log transformation
df_qn_log <- quantile_normalization(df, log_trans = TRUE)
df_qn_log[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166  0.0000000  0.0000000    0.00000    0.00000    0.00000
#> ENSG00000266931  0.0000000  0.0000000    0.00000    0.00000    0.00000
#> ENSG00000104774 11.0807343 11.1801423   10.92769   10.94581   10.65441
#> ENSG00000267583  0.0000000  0.8073549    0.00000    0.00000    0.00000
#> ENSG00000227581  0.7004397  0.0000000    0.00000    0.00000    0.00000

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

## Overwrite existing assay
se2 <- quantile_normalization(se)
assay(se2)[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166      0.000      0.000      0.000        0.0       0.00
#> ENSG00000266931      0.000      0.000      0.000        0.0       0.00
#> ENSG00000104774   2164.875   2319.375   1946.875     1971.5    1610.75
#> ENSG00000267583      0.000      0.750      0.000        0.0       0.00
#> ENSG00000227581      0.625      0.000      0.000        0.0       0.00

## Store result in new assay
se3 <- quantile_normalization(se, new_assay_name = "quant_norm")
assay(se3, "quant_norm")[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166      0.000      0.000      0.000        0.0       0.00
#> ENSG00000266931      0.000      0.000      0.000        0.0       0.00
#> ENSG00000104774   2164.875   2319.375   1946.875     1971.5    1610.75
#> ENSG00000267583      0.000      0.750      0.000        0.0       0.00
#> ENSG00000227581      0.625      0.000      0.000        0.0       0.00

## Use specific input assay (simulate new one)
new_matrix <- matrix(
  data = sample(x = seq(1, 100000), size = nrow(se) * ncol(se),
  replace = TRUE),
  nrow = nrow(se),
  ncol = ncol(se)
)
rownames(new_matrix) <- rownames(se)
colnames(new_matrix) <- colnames(se)

## Create a new assay in the SummarizedExperiment
assay(se, "new_counts") <- new_matrix

## Normalize the new assay and store it under a new name
se4 <- quantile_normalization(se, assay_name = "new_counts",
new_assay_name = "quant_new")
assay(se4, "quant_new")[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166   56060.00   94285.88   82647.62   60109.62  32567.250
#> ENSG00000266931   54882.62   76292.25   79960.25   90485.50   3259.625
#> ENSG00000104774     733.25   40996.50    9360.25   58418.00  91851.750
#> ENSG00000267583   19468.50   65204.88   51248.62   78804.12  24111.875
#> ENSG00000227581   61854.62   21166.88   95707.25   40996.50  98769.375
```
