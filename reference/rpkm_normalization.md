# RPKM Normalization

Normalizes read counts by the RPKM (Reads Per Kilobase per Million
mapped reads) method:

1.  Normalize counts by library size (column sums), scaled to millions.

2.  Divide each gene's value by its length in kilobases.

3.  If `log_trans = TRUE`, applies `log2(RPKM + 1)`.

## Usage

``` r
rpkm_normalization(
  x,
  gene_length = NULL,
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

  :   the function retrieves `gene_length` from
      `rowData(x)$gene_length`.

  If a `data.frame`/`matrix`,

  :   the user must provide the `gene_length` argument.

- gene_length:

  A numeric vector of gene lengths (one per row), used only if `x` is a
  data.frame or matrix. Must match the number of rows in `x`. Ignored if
  `x` is a SummarizedExperiment.

- log_trans:

  Logical. If `TRUE`, apply `log2(... + 1)` transform to the
  RPKM-normalized values.

- assay_name:

  If `x` is a SummarizedExperiment, name of the assay to normalize.
  Defaults to the first assay if not specified.

- new_assay_name:

  If `x` is a SummarizedExperiment, name of a new assay in which to
  store the RPKM (or log2-RPKM). If `NULL`, overwrites the assay
  specified in `assay_name`.

## Value

A numeric **matrix** of RPKM or log2(RPKM + 1) values if `x` is a
data.frame or matrix. If `x` is a SummarizedExperiment, returns the
modified SummarizedExperiment with the RPKM data placed in the existing
or new assay.

## Details

If `x` is a `SummarizedExperiment`, the function looks for a numeric
column named `"gene_length"` in `rowData(x)`. That column must have
length equal to the number of rows in the assay being normalized.

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

### Adding a column in rowData regarding the gene_length
rowData(se)$gene_length = rowData(se)$gene_seq_end -
rowData(se)$gene_seq_start

# -------------------------------
# 1) Using a data.frame
# -------------------------------

gene_length = rowData(se)$gene_length
df = assay(se)

## Without log transformation
df = rpkm_normalization(df, gene_length = gene_length)
df[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166   0.000000   0.000000     0.0000     0.0000     0.0000
#> ENSG00000266931   0.000000   0.000000     0.0000     0.0000     0.0000
#> ENSG00000104774 277.787316 342.353267   261.1592   254.7303   226.1147
#> ENSG00000267583   0.000000   0.156504     0.0000     0.0000     0.0000
#> ENSG00000227581   4.667242   0.000000     0.0000     0.0000     0.0000

## With log transformation
df = rpkm_normalization(df, gene_length = gene_length, log_trans = TRUE)
df[1:5, 1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166   0.000000 0.00000000   0.000000   0.000000   0.000000
#> ENSG00000266931   0.000000 0.00000000   0.000000   0.000000   0.000000
#> ENSG00000104774   6.369221 6.63138416   6.230741   6.186639   6.195829
#> ENSG00000267583   0.000000 0.05532084   0.000000   0.000000   0.000000
#> ENSG00000227581   5.514145 0.00000000   0.000000   0.000000   0.000000

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

# If no new_assay_name is provided, then overwrites existing assay
se2 = rpkm_normalization(se, log_trans = FALSE)
head(assay(se2))
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166   0.000000   0.000000     0.0000     0.0000     0.0000
#> ENSG00000266931   0.000000   0.000000     0.0000     0.0000     0.0000
#> ENSG00000104774 277.787316 342.353267   261.1592   254.7303   226.1147
#> ENSG00000267583   0.000000   0.156504     0.0000     0.0000     0.0000
#> ENSG00000227581   4.667242   0.000000     0.0000     0.0000     0.0000
#> ENSG00000227317   0.000000   0.000000     0.0000     0.0000     0.0000
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000260166   0.000000   5.397845     0.0000
#> ENSG00000266931   0.000000   0.000000     0.0000
#> ENSG00000104774 231.402716 263.922352   288.9925
#> ENSG00000267583   0.000000   0.000000     0.0000
#> ENSG00000227581   6.668431   0.000000     0.0000
#> ENSG00000227317   0.000000   0.000000     0.0000

# If new_assay_name is given, normalization stored in a new assay
se2 = rpkm_normalization(se, log_trans = FALSE, new_assay_name =
'rpkm_counts')
head(assay(se2, 'rpkm_counts'))
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166   0.000000   0.000000     0.0000     0.0000     0.0000
#> ENSG00000266931   0.000000   0.000000     0.0000     0.0000     0.0000
#> ENSG00000104774 277.787316 342.353267   261.1592   254.7303   226.1147
#> ENSG00000267583   0.000000   0.156504     0.0000     0.0000     0.0000
#> ENSG00000227581   4.667242   0.000000     0.0000     0.0000     0.0000
#> ENSG00000227317   0.000000   0.000000     0.0000     0.0000     0.0000
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000260166   0.000000   5.397845     0.0000
#> ENSG00000266931   0.000000   0.000000     0.0000
#> ENSG00000104774 231.402716 263.922352   288.9925
#> ENSG00000267583   0.000000   0.000000     0.0000
#> ENSG00000227581   6.668431   0.000000     0.0000
#> ENSG00000227317   0.000000   0.000000     0.0000

# Creating a new assay to test specific input
new_matrix = matrix(data = sample(x = seq(1, 100000),
                                  size = nrow(se) * ncol(se),
                                  replace = TRUE),
                    nrow = nrow(se),
                    ncol = ncol(se))
rownames(new_matrix) = rownames(se)
colnames(new_matrix) = colnames(se)

assay(se, 'new_counts') = new_matrix
se2 = rpkm_normalization(se, new_assay_name = 'rpkm_counts_new',
assay_name = 'new_counts')
head(assay(se2, 'rpkm_counts_new'))
#>                   SRR1039508  SRR1039509  SRR1039512  SRR1039513 SRR1039516
#> ENSG00000260166 2094.2629240  3469.37396  3018.83571  2183.35424 1159.58044
#> ENSG00000266931 7897.3713086 10924.41315 11186.13827 12847.91118  431.24651
#> ENSG00000104774    0.6815632    39.12036    10.59334    56.94763   89.50011
#> ENSG00000267583   15.5723767    55.70475    43.28536    67.78891   20.31845
#> ENSG00000227581 2015.7901382   715.67622  3054.89145  1324.78648 3101.77395
#> ENSG00000227317  342.1502083     1.13400   519.97698   544.10881  226.11313
#>                 SRR1039517  SRR1039520  SRR1039521
#> ENSG00000260166 3563.20745  639.188949  400.674218
#> ENSG00000266931 1277.44493 9424.107099 1061.567274
#> ENSG00000104774   30.94395    3.081551    5.461995
#> ENSG00000267583   57.18676   46.994369    4.034666
#> ENSG00000227581 3095.59814 1252.049146 2844.308539
#> ENSG00000227317  546.65211  430.165035  186.571110
```
