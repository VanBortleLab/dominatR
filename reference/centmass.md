# Compute the "center of mass" for rows of a data frame or SummarizedExperiment

For each row of the numeric data, `centmass()` computes a 2D center of
mass with coordinates (`comx`, `comy`). The `x_coord` and `y_coord`
vectors specify the location for each column's "mass."

The original usage assumes a ternary coordinate system by default, but
this can be generalized to any scenario where columns represent discrete
"masses" at known (x,y) positions.

By default, `x_coord = c(0, 1, 0.5)` and `y_coord = c(0, 0, sqrt(3)/2)`,
which correspond to the corners of an equilateral triangle (often used
in ternary plots).

## Usage

``` r
centmass(
  x,
  x_coord = c(0, 1, 0.5),
  y_coord = c(0, 0, sqrt(3)/2),
  assay_name = NULL
)
```

## Arguments

- x:

  A data.frame (with numeric columns) or a SummarizedExperiment.

- x_coord:

  Numeric vector of length equal to the number of columns in `x`,
  specifying the x-coordinates of each column's mass.

- y_coord:

  Numeric vector of length equal to the number of columns in `x`,
  specifying the y-coordinates of each column's mass.

- assay_name:

  If `x` is a SummarizedExperiment, the name of the assay to use.
  Defaults to the first assay if not specified.

## Value

- If `x` is a data.frame, returns a new `data.frame` with columns `comx`
  and `comy`.

- If `x` is a SummarizedExperiment, returns the same object but with two
  new columns `comx` and `comy` in `rowData(x)`.

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

# Let's subset for the first 3 columns for this example
se = se[,1:3]
# -------------------------------
# 1) Using a data.frame
# -------------------------------


df = assay(se) |> as.data.frame()

df = centmass(df)
head(df)
#>                      comx      comy
#> ENSG00000260166 0.0000000 0.0000000
#> ENSG00000266931 0.0000000 0.0000000
#> ENSG00000104774 0.4973123 0.3236133
#> ENSG00000267583 1.0000000 0.0000000
#> ENSG00000227581 0.0000000 0.0000000
#> ENSG00000227317 0.0000000 0.0000000

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

se2 = centmass(se)

## X and Y coordinates are stored in rowData(se2)
head(rowData(se2))
#> DataFrame with 6 rows and 12 columns
#>                         gene_id      gene_name  entrezid   gene_biotype
#>                     <character>    <character> <integer>    <character>
#> ENSG00000260166 ENSG00000260166  RP11-863P13.6        NA        lincRNA
#> ENSG00000266931 ENSG00000266931 RP11-1252D15.1        NA     pseudogene
#> ENSG00000104774 ENSG00000104774         MAN2B1        NA protein_coding
#> ENSG00000267583 ENSG00000267583  RP11-322E11.5        NA        lincRNA
#> ENSG00000227581 ENSG00000227581   RP13-140E4.1        NA     pseudogene
#> ENSG00000227317 ENSG00000227317          DDAH2        NA protein_coding
#>                 gene_seq_start gene_seq_end       seq_name seq_strand
#>                      <integer>    <integer>    <character>  <integer>
#> ENSG00000260166       88120993     88121538             16          1
#> ENSG00000266931       87282781     87282923              2         -1
#> ENSG00000104774       12757325     12777556             19         -1
#> ENSG00000267583       33023833     33047052             18         -1
#> ENSG00000227581       89294212     89294833              X         -1
#> ENSG00000227317       31676987     31680566 HSCHR6_MHC_DBB         -1
#>                 seq_coord_system         symbol      comx      comy
#>                        <integer>    <character> <numeric> <numeric>
#> ENSG00000260166               NA  RP11-863P13.6  0.000000  0.000000
#> ENSG00000266931               NA RP11-1252D15.1  0.000000  0.000000
#> ENSG00000104774               NA         MAN2B1  0.497312  0.323613
#> ENSG00000267583               NA  RP11-322E11.5  1.000000  0.000000
#> ENSG00000227581               NA   RP13-140E4.1  0.000000  0.000000
#> ENSG00000227317               NA          DDAH2  0.000000  0.000000

```
