# Counts Per Million normalization

Normalizes a count matrix (or a SummarizedExperiment assay) by the
counts-per-million (CPM) method. Specifically:

1.  If `log_trans = TRUE`, a `log2(x + 1)` transform is applied
    afterward.

## Usage

``` r
cpm_normalization(
  x,
  log_trans = FALSE,
  assay_name = NULL,
  new_assay_name = NULL
)
```

## Arguments

- x:

  A `matrix`, `data.frame`, or a `SummarizedExperiment` object.

- log_trans:

  Logical. If `TRUE`, apply `log2(... + 1)` transform to the
  CPM-normalized values.

- assay_name:

  If `x` is a `SummarizedExperiment`, name of the assay to normalize
  (defaults to the first assay). Ignored otherwise.

- new_assay_name:

  If `x` is a `SummarizedExperiment`, name of a new assay where results
  should be stored (defaults to `NULL`, meaning the existing assay is
  overwritten).

## Value

- If `x` is a `matrix` or `data.frame`, returns a **matrix** of
  CPM-normalized (and optionally `log2`-transformed) counts.

- If `x` is a `SummarizedExperiment`, returns the same
  `SummarizedExperiment` object with the specified assay replaced or a
  new assay created containing the CPM-normalized data.

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

df = assay(se)

## Without log transformation
df1 = cpm_normalization(df, log_trans = FALSE)

df1[1:5,1:5]
#>                  SRR1039508  SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166    0.000000    0.000000      0.000      0.000      0.000
#> ENSG00000266931    0.000000    0.000000      0.000      0.000      0.000
#> ENSG00000104774 5619.915194 6926.148938   5283.511   5153.449   4574.527
#> ENSG00000267583    0.000000    3.633866      0.000      0.000      0.000
#> ENSG00000227581    2.898358    0.000000      0.000      0.000      0.000

## With log transformation
df1 = cpm_normalization(df, log_trans = TRUE)

df1[1:5,1:5]
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166   0.000000   0.000000    0.00000     0.0000    0.00000
#> ENSG00000266931   0.000000   0.000000    0.00000     0.0000    0.00000
#> ENSG00000104774  12.456589  12.758046   12.36755    12.3316   12.15972
#> ENSG00000267583   0.000000   2.212216    0.00000     0.0000    0.00000
#> ENSG00000227581   1.962866   0.000000    0.00000     0.0000    0.00000

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

# If now new_assay_name is provided, then overwrites existing assay
se2 = cpm_normalization(se, log_trans = FALSE)

se2
#> class: RangedSummarizedExperiment 
#> dim: 1000 8 
#> metadata(1): ''
#> assays(1): counts
#> rownames(1000): ENSG00000260166 ENSG00000266931 ... ENSG00000160886
#>   ENSG00000142871
#> rowData names(10): gene_id gene_name ... seq_coord_system symbol
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
head(assay(se2))
#>                  SRR1039508  SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166    0.000000    0.000000      0.000      0.000      0.000
#> ENSG00000266931    0.000000    0.000000      0.000      0.000      0.000
#> ENSG00000104774 5619.915194 6926.148938   5283.511   5153.449   4574.527
#> ENSG00000267583    0.000000    3.633866      0.000      0.000      0.000
#> ENSG00000227581    2.898358    0.000000      0.000      0.000      0.000
#> ENSG00000227317    0.000000    0.000000      0.000      0.000      0.000
#>                  SRR1039517  SRR1039520 SRR1039521
#> ENSG00000260166    0.000000    2.941825      0.000
#> ENSG00000266931    0.000000    0.000000      0.000
#> ENSG00000104774 4681.508353 5339.413106   5846.607
#> ENSG00000267583    0.000000    0.000000      0.000
#> ENSG00000227581    4.141095    0.000000      0.000
#> ENSG00000227317    0.000000    0.000000      0.000

# If new new_assay_name, normalization stored in a new object
se2 = cpm_normalization(se, log_trans = FALSE, new_assay_name = 'cpm_counts')

se2
#> class: RangedSummarizedExperiment 
#> dim: 1000 8 
#> metadata(1): ''
#> assays(2): counts cpm_counts
#> rownames(1000): ENSG00000260166 ENSG00000266931 ... ENSG00000160886
#>   ENSG00000142871
#> rowData names(10): gene_id gene_name ... seq_coord_system symbol
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
head(assay(se2, 'cpm_counts'))
#>                  SRR1039508  SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166    0.000000    0.000000      0.000      0.000      0.000
#> ENSG00000266931    0.000000    0.000000      0.000      0.000      0.000
#> ENSG00000104774 5619.915194 6926.148938   5283.511   5153.449   4574.527
#> ENSG00000267583    0.000000    3.633866      0.000      0.000      0.000
#> ENSG00000227581    2.898358    0.000000      0.000      0.000      0.000
#> ENSG00000227317    0.000000    0.000000      0.000      0.000      0.000
#>                  SRR1039517  SRR1039520 SRR1039521
#> ENSG00000260166    0.000000    2.941825      0.000
#> ENSG00000266931    0.000000    0.000000      0.000
#> ENSG00000104774 4681.508353 5339.413106   5846.607
#> ENSG00000267583    0.000000    0.000000      0.000
#> ENSG00000227581    4.141095    0.000000      0.000
#> ENSG00000227317    0.000000    0.000000      0.000

# A specific assay can also be selected
new_matrix =  matrix(data = sample(x = seq(1, 100000),
                                  size = nrow(se) * ncol(se),
                                  replace = TRUE),
                    nrow = nrow(se),
                    ncol = ncol(se))
rownames(new_matrix) = rownames(se)
colnames(new_matrix) = colnames(se)

## Creating a new assay called new counts
assay(se, 'new_counts') = new_matrix

se2 = cpm_normalization(se, new_assay_name = 'cpm_counts_new', assay_name =
'new_counts')

se2
#> class: RangedSummarizedExperiment 
#> dim: 1000 8 
#> metadata(1): ''
#> assays(3): counts new_counts cpm_counts_new
#> rownames(1000): ENSG00000260166 ENSG00000266931 ... ENSG00000160886
#>   ENSG00000142871
#> rowData names(10): gene_id gene_name ... seq_coord_system symbol
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
head(assay(se2, 'cpm_counts_new'))
#>                 SRR1039508  SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166  1141.3733 1890.808811  1645.2655  1189.9281  631.97134
#> ENSG00000266931  1121.4267 1551.266667  1588.4316  1824.4034   61.23701
#> ENSG00000104774    13.7887  791.444090   214.3138  1152.1075 1810.67663
#> ENSG00000267583   361.5750 1293.408505  1005.0427  1573.9907  471.77409
#> ENSG00000227581  1251.8057  444.434935  1897.0876   822.6924 1926.20162
#> ENSG00000227317  1224.5556    4.058587  1860.9976  1947.3654  809.25890
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000260166  1941.9481  348.35798   218.3674
#> ENSG00000266931   181.3972 1338.22321   150.7426
#> ENSG00000104774   626.0270   62.34286   110.5016
#> ENSG00000267583  1327.8193 1091.16225    93.6809
#> ENSG00000227581  1922.3664  777.52252  1766.3156
#> ENSG00000227317  1956.4679 1539.56066   667.7380
```
