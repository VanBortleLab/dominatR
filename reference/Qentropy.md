# Compute Q-Entropy using existing row-normalized data + Entropy

\#' Transform entropy scores into categorical entropy scores \\Q\_{ij} =
\mathrm{Entropy}\_i - \log_2(x\_{ij})\\, or `Inf` if \\x\_{ij} == 0\\.

@details For each row \\i\\ and column \\j\\, \\Q\_{ij}\\ is defined as
\\\mathrm{Entropy}\_i - \log_2\bigl(x\_{ij}\bigr)\\ if \\x\_{ij}\\ is
positive, or `Inf` otherwise.

## Usage

``` r
Qentropy(x, assay_name = "Entropy", new_assay_name = "Qentropy")
```

## Arguments

- x:

  A data.frame (already processed by 'entropy()') or a
  SummarizedExperiment (already processed by 'entropy()').

- assay_name:

  (SummarizedExperiment only) The name of the assay whose row-normalized
  data will be replaced by Q-values. If NULL, uses the first assay.

- new_assay_name:

  If you prefer to store Q-values in a \*new\* assay, provide a name. By
  default 'Qentropy'

## Value

- If `x` is a data.frame: returns the same data.frame with numeric
  columns replaced by \\Q\_{ij}\\ values and `Entropy` column removed.

- If `x` is a SummarizedExperiment: returns the same object with the
  specified assay replaced by \\Q\_{ij}\\ values (or a new assay if
  `new_assay_name` is set) and `rowData(x)$Entropy` removed.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
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

## Entropy needs to be calculated first
df = entropy(df)

## Then you can apply the Qentropy function
df = Qentropy(df)

head(df)
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000260166        Inf        Inf        Inf        Inf        Inf
#> ENSG00000266931        Inf        Inf        Inf        Inf        Inf
#> ENSG00000104774   5.954307   5.979072   5.711755   6.640587   6.005811
#> ENSG00000267583        Inf   0.000000        Inf        Inf        Inf
#> ENSG00000227581   2.503258        Inf        Inf        Inf        Inf
#> ENSG00000227317        Inf        Inf        Inf        Inf        Inf
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000260166        Inf   0.000000        Inf
#> ENSG00000266931        Inf        Inf        Inf
#> ENSG00000104774   5.732659   6.049651   5.941714
#> ENSG00000267583        Inf        Inf        Inf
#> ENSG00000227581   1.503258        Inf        Inf
#> ENSG00000227317        Inf        Inf        Inf

# -------------------------------
# 2) Using a SummarizedExperiment
# -------------------------------

## Calculate Entropy first
se2 = entropy(se, new_assay_name = 'Entropy')

## Transform entropy into Qentropy. new_assay_name specify a new assay
## where data is going to be stored. Assay_name must have Entropy transformed
values
#> new("standardGeneric", .Data = function (x, ...) 
#> standardGeneric("values"), generic = "values", package = "S4Vectors", 
#>     group = list(), valueClass = character(0), signature = "x", 
#>     default = NULL, skeleton = (function (x, ...) 
#>     stop(gettextf("invalid call in method dispatch to '%s' (no default method)", 
#>         "values"), domain = NA))(x, ...))
#> <bytecode: 0x561e8c4c0e78>
#> <environment: 0x561e8c4bc830>
#> attr(,"generic")
#> [1] "values"
#> attr(,"generic")attr(,"package")
#> [1] "S4Vectors"
#> attr(,"package")
#> [1] "S4Vectors"
#> attr(,"group")
#> list()
#> attr(,"valueClass")
#> character(0)
#> attr(,"signature")
#> [1] "x"
#> attr(,"default")
#> `\001NULL\001`
#> attr(,"skeleton")
#> (function (x, ...) 
#> stop(gettextf("invalid call in method dispatch to '%s' (no default method)", 
#>     "values"), domain = NA))(x, ...)
#> attr(,"class")
#> [1] "standardGeneric"
#> attr(,"class")attr(,"package")
#> [1] "methods"
## By default, the function will look for an assay_name 'Entropy' and assign
## a new assay to 'Qentropy'
se2 = Qentropy(se2, new_assay_name = 'Qentropy', assay_name = 'Entropy')

se2
#> class: RangedSummarizedExperiment 
#> dim: 1000 8 
#> metadata(1): ''
#> assays(3): counts Entropy Qentropy
#> rownames(1000): ENSG00000260166 ENSG00000266931 ... ENSG00000160886
#>   ENSG00000142871
#> rowData names(11): gene_id gene_name ... symbol Entropy
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(9): SampleName cell ... Sample BioSample
```
