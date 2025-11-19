# Normalization Functions

## Data

The data used in this article comes from only from the following source:

1.  A `SummarizedExperiment` object from the library `airway`

In this case we are using a single object given that it can be
manipulated as a `SummarizedExperiment` or `data.frame`

``` r
### summarized experiment
data("airway")
se <-airway
```

## Normalization Methods

Normalization is critical for correcting technical biases and enabling
meaningful biological comparisons.

The package contains different normalization methods, some of them
include a `log_transf` attribute that calculates the `log2(x+1)` of the
normalized value if set to `TRUE`:

- `cpm_normalization`
- `minmax_normalization`
- `quantile_normalization`
- `rpkm_normalization`
- `tpm_normalization`

Let’s explore the usage of each normalization method on the `airway`
data.

### Min-Max Normalization

Min-Max normalization is a linear transformation technique that rescales
each gene’s expression values to a specified range (typically \[0, 1\]).
This normalization method is useful when you want to bring the data onto
the same scale.

Function Purpose:

· Rescales each column to fit within a range \[new_min, new_max\].

· Preserves the relative structure of values within each column.

· Useful when different assays or samples have varying scales.

#### Example 1: Normalize a matrix

``` r
# Prepare input matrix
count_mat <- assay(se)

# Apply min-max normalization
se_minmax <- minmax_normalization(count_mat, new_min = 0, new_max = 1)

# Inspect structure
dim(se_minmax)
#> [1] 63677     8
summary(as.vector(se_minmax))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 0.000e+00 0.000e+00 9.679e-04 2.739e-05 1.000e+00
head(se_minmax[, 1:5])
#>                   SRR1039508   SRR1039509   SRR1039512   SRR1039513
#> ENSG00000000003 0.0022792424 0.0017523136 1.699217e-03 0.0014897144
#> ENSG00000000005 0.0000000000 0.0000000000 0.000000e+00 0.0000000000
#> ENSG00000000419 0.0015676086 0.0020143784 1.208721e-03 0.0013327102
#> ENSG00000000457 0.0008727585 0.0008253084 5.119062e-04 0.0005988068
#> ENSG00000000460 0.0002014058 0.0002151278 7.785646e-05 0.0001277941
#> ENSG00000000938 0.0000000000 0.0000000000 3.892823e-06 0.0000000000
#>                   SRR1039516
#> ENSG00000000003 2.860799e-03
#> ENSG00000000005 0.000000e+00
#> ENSG00000000419 1.475649e-03
#> ENSG00000000457 6.159013e-04
#> ENSG00000000460 1.960829e-04
#> ENSG00000000938 2.513883e-06
```

You can set new_min = 10 and new_max = 20 if your downstream application
prefers values in a different scale:

``` r
df_scaled <- minmax_normalization(count_mat, new_min = 10, new_max = 20)
head(df_scaled)  # All columns now range from 10 to 20
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003   10.02279   10.01752   10.01699   10.01490   10.02861
#> ENSG00000000005   10.00000   10.00000   10.00000   10.00000   10.00000
#> ENSG00000000419   10.01568   10.02014   10.01209   10.01333   10.01476
#> ENSG00000000457   10.00873   10.00825   10.00512   10.00599   10.00616
#> ENSG00000000460   10.00201   10.00215   10.00078   10.00128   10.00196
#> ENSG00000000938   10.00000   10.00000   10.00004   10.00000   10.00003
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000000003   10.02607   10.02033   10.01536
#> ENSG00000000005   10.00000   10.00000   10.00000
#> ENSG00000000419   10.01990   10.01101   10.01364
#> ENSG00000000457   10.00824   10.00615   10.00615
#> ENSG00000000460   10.00157   10.00201   10.00161
#> ENSG00000000938   10.00000   10.00000   10.00000
```

#### Example 2: Normalize a SummarizedExperiment

``` r
se <- se

# Option A: Overwrite the default assay
se1 <- minmax_normalization(se)
head(assay(se1))
#>                   SRR1039508   SRR1039509   SRR1039512   SRR1039513
#> ENSG00000000003 0.0022792424 0.0017523136 1.699217e-03 0.0014897144
#> ENSG00000000005 0.0000000000 0.0000000000 0.000000e+00 0.0000000000
#> ENSG00000000419 0.0015676086 0.0020143784 1.208721e-03 0.0013327102
#> ENSG00000000457 0.0008727585 0.0008253084 5.119062e-04 0.0005988068
#> ENSG00000000460 0.0002014058 0.0002151278 7.785646e-05 0.0001277941
#> ENSG00000000938 0.0000000000 0.0000000000 3.892823e-06 0.0000000000
#>                   SRR1039516   SRR1039517   SRR1039520   SRR1039521
#> ENSG00000000003 2.860799e-03 0.0026074678 0.0020325525 0.0015356158
#> ENSG00000000005 0.000000e+00 0.0000000000 0.0000000000 0.0000000000
#> ENSG00000000419 1.475649e-03 0.0019898441 0.0011007460 0.0013637987
#> ENSG00000000457 6.159013e-04 0.0008243284 0.0006150451 0.0006147833
#> ENSG00000000460 1.960829e-04 0.0001568963 0.0002006156 0.0001610786
#> ENSG00000000938 2.513883e-06 0.0000000000 0.0000000000 0.0000000000

# Option B: Write to a new assay slot
se2 <- minmax_normalization(se, new_assay_name = "minmax_counts")
```

By using the option `new_assay_name` it is possible to store the
normalized data in a new assay in the summarizedexperiment object
keeping the count matrix intact. If no name is provided upon
normalization, then the function will overwrite the count matrix

### Quantile Normalization

Quantile normalization makes the distribution of values across all
samples identical. This technique adjusts the data so that the rank
distributions of the data across samples are equal.

#### Example 1: Normalize a matrix

``` r
count_mat <- assay(se)

se_quantile <- quantile_normalization(count_mat)

## Check result
dim((se_quantile))
#> [1] 63677     8
summary(as.vector(se_quantile))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#>      0.00      0.00      0.00    344.36      9.62 361483.12
head(se_quantile[1:5, 1:5])
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003    690.875    504.750    773.875     613.75   1010.000
#> ENSG00000000005      0.000      0.000      0.000       0.00      0.000
#> ENSG00000000419    468.875    582.375    550.625     552.00    516.875
#> ENSG00000000457    257.375    241.375    225.250     254.00    213.125
#> ENSG00000000460     58.000     65.250     31.500      53.75     67.750
```

#### Example 2: Normalize a summarized experiment

``` r
## Apply quantile normalization to a SummarizedExperiment
se_quantile <- quantile_normalization(se)

## Check result
dim(assay(se_quantile))
#> [1] 63677     8
summary(as.vector(assay(se_quantile)))
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#>      0.00      0.00      0.00    344.36      9.62 361483.12
head(assay(se_quantile)[1:5, 1:5])
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003    690.875    504.750    773.875     613.75   1010.000
#> ENSG00000000005      0.000      0.000      0.000       0.00      0.000
#> ENSG00000000419    468.875    582.375    550.625     552.00    516.875
#> ENSG00000000457    257.375    241.375    225.250     254.00    213.125
#> ENSG00000000460     58.000     65.250     31.500      53.75     67.750
```

### CPM Normalization

The cpm_normalization() function rescales raw count data such that each
column sums to one million. This makes count data comparable across
samples of different sequencing depths.

#### Example 1: Normalize a data.frame

``` r
df <- assay(se)
# Normalize without log2-transform
df_cpm <- cpm_normalization(df, log_trans = FALSE)
head(df_cpm[, 1:5])
#>                 SRR1039508 SRR1039509  SRR1039512 SRR1039513  SRR1039516
#> ENSG00000000003  32.900521  23.817776 34.43970525  26.906868 46.54699807
#> ENSG00000000005   0.000000   0.000000  0.00000000   0.000000  0.00000000
#> ENSG00000000419  22.628193  27.379809 24.49834703  24.071095 24.00974329
#> ENSG00000000457  12.598138  11.217747 10.37530639  10.815506 10.02110240
#> ENSG00000000460   2.907263   2.924057  1.57799337   2.308187  3.19039178
#> ENSG00000000938   0.000000   0.000000  0.07889967   0.000000  0.04090246

# Normalize with log2-transform
df_cpm_log <- cpm_normalization(df, log_trans = TRUE)
head(df_cpm_log[, 1:5])
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003   5.083236   4.633302  5.1472947   4.802548 5.57128235
#> ENSG00000000005   0.000000   0.000000  0.0000000   0.000000 0.00000000
#> ENSG00000000419   4.562437   4.826793  4.6723318   4.647953 4.64441834
#> ENSG00000000457   3.765337   3.610906  3.5078335   3.562609 3.46219663
#> ENSG00000000460   1.966158   1.972346  1.3662486   1.726041 2.06708514
#> ENSG00000000938   0.000000   0.000000  0.1095607   0.000000 0.05783488
```

#### Example 2: Normalize a SummarizedExperiment

``` r

# Apply in-place normalization (overwrite assay)
se1 <- cpm_normalization(se, log_trans = FALSE)
head(assay(se1))
#>                 SRR1039508 SRR1039509  SRR1039512 SRR1039513  SRR1039516
#> ENSG00000000003  32.900521  23.817776 34.43970525  26.906868 46.54699807
#> ENSG00000000005   0.000000   0.000000  0.00000000   0.000000  0.00000000
#> ENSG00000000419  22.628193  27.379809 24.49834703  24.071095 24.00974329
#> ENSG00000000457  12.598138  11.217747 10.37530639  10.815506 10.02110240
#> ENSG00000000460   2.907263   2.924057  1.57799337   2.308187  3.19039178
#> ENSG00000000938   0.000000   0.000000  0.07889967   0.000000  0.04090246
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000000003  33.973415  40.259015  27.026857
#> ENSG00000000005   0.000000   0.000000   0.000000
#> ENSG00000000419  25.926226  21.802609  24.002873
#> ENSG00000000457  10.740401  12.182273  10.820193
#> ENSG00000000460   2.044246   3.973617   2.834985
#> ENSG00000000938   0.000000   0.000000   0.000000

# Save to a new assay slot
se2 <- cpm_normalization(se, log_trans = TRUE, new_assay_name = 
                            "cpm_logged")
head(assay(se2, "cpm_logged"))
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003   5.083236   4.633302  5.1472947   4.802548 5.57128235
#> ENSG00000000005   0.000000   0.000000  0.0000000   0.000000 0.00000000
#> ENSG00000000419   4.562437   4.826793  4.6723318   4.647953 4.64441834
#> ENSG00000000457   3.765337   3.610906  3.5078335   3.562609 3.46219663
#> ENSG00000000460   1.966158   1.972346  1.3662486   1.726041 2.06708514
#> ENSG00000000938   0.000000   0.000000  0.1095607   0.000000 0.05783488
#>                 SRR1039517 SRR1039520 SRR1039521
#> ENSG00000000003   5.128187   5.366637   4.808738
#> ENSG00000000005   0.000000   0.000000   0.000000
#> ENSG00000000419   4.750940   4.511127   4.644022
#> ENSG00000000457   3.553410   3.720527   3.563182
#> ENSG00000000460   1.606085   2.314295   1.939221
#> ENSG00000000938   0.000000   0.000000   0.000000
```

### RPKM Normalization

Reads per kilobase per million (RPKM) normalization adjusts for both
gene length and sequencing depth, making it particularly useful for
RNA-Seq data. RPKM helps compare gene expression levels across genes of
different lengths.

#### Example 1: Normalize a data.frame

``` r
df <- assay(se)
length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start
df_rpkm <- rpkm_normalization(df, gene_length = length)

head(df_rpkm[, 1:5])
#>                 SRR1039508 SRR1039509  SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003 2.90614973 2.10385794 3.042108051 2.37672181 4.11156241
#> ENSG00000000005 0.00000000 0.00000000 0.000000000 0.00000000 0.00000000
#> ENSG00000000419 0.95525977 1.15585145 1.034209179 1.01617253 1.01358254
#> ENSG00000000457 0.28224164 0.25131614 0.232442566 0.24230454 0.22450718
#> ENSG00000000460 0.01514389 0.01523137 0.008219743 0.01202331 0.01661870
#> ENSG00000000938 0.00000000 0.00000000 0.003398943 0.00000000 0.00176205
```

#### Example 2: Normalize a SummarizedExperiment

``` r
## Gene length needed
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start

## Apply RPKM normalization
se_rpkm <- rpkm_normalization(se, gene_length, log_trans = TRUE)

## Check the result
dim(assay(se_rpkm)) 
#> [1] 63677     8
head(assay(se_rpkm)[1:5, 1:5])
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003 1.96574725 1.63406253 2.01510789 1.75562333 2.35376434
#> ENSG00000000005 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> ENSG00000000419 0.96736029 1.10825777 1.02446804 1.01161910 1.00976461
#> ENSG00000000457 0.35866816 0.32344632 0.30152042 0.31301889 0.29220123
#> ENSG00000000460 0.02168423 0.02180855 0.01181011 0.01724252 0.02377868
```

### TPM Normalization

Transcripts per million normalization.

#### Example 1: Normalize a data.frame

``` r
df <- assay(se)
length <- sample(c(400:800), nrow(df), replace = TRUE)
df_tpm <- tpm_normalization(df, gene_length = length)

head(df_tpm[, 1:5])
#>                 SRR1039508 SRR1039509  SRR1039512 SRR1039513  SRR1039516
#> ENSG00000000003  33.390121  24.122006 34.96289802  27.202754 47.35784714
#> ENSG00000000005   0.000000   0.000000  0.00000000   0.000000  0.00000000
#> ENSG00000000419  17.352628  20.952835 18.79251602  18.388476 18.45814018
#> ENSG00000000457   9.547481   8.483696  7.86531621   8.165143  7.61347687
#> ENSG00000000460   3.784083   3.798038  2.05453877   2.992829  4.16299091
#> ENSG00000000938   0.000000   0.000000  0.07636036   0.000000  0.03967295
```

#### Example 2: Normalize a SummarizedExperiment

``` r
## Gene length needed
rowData(se)$gene_length <- rowData(se)$gene_seq_end - rowData(se)$gene_seq_start

## Apply RPKM normalization
se_tpm <- tpm_normalization(se, gene_length, log_trans = TRUE)

## Check the result
dim(assay(se_tpm)) 
#> [1] 63677     8
head(assay(se_tpm)[1:5, 1:5])
#>                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
#> ENSG00000000003  4.4294986  3.9577411 4.55434051  4.2249660  4.8647104
#> ENSG00000000005  0.0000000  0.0000000 0.00000000  0.0000000  0.0000000
#> ENSG00000000419  2.9549906  3.1678707 3.11233806  3.0989031  2.9884095
#> ENSG00000000457  1.5828547  1.4524115 1.44301228  1.4877425  1.3427320
#> ENSG00000000460  0.1467549  0.1443756 0.08513068  0.1237197  0.1553897
```

#### Session Info

    #> R version 4.5.2 (2025-10-31)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.3 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    #> [8] base     
    #> 
    #> other attached packages:
    #>  [1] ggplot2_4.0.1               dominatRData_0.99.1        
    #>  [3] airway_1.30.0               SummarizedExperiment_1.40.0
    #>  [5] Biobase_2.70.0              GenomicRanges_1.62.0       
    #>  [7] Seqinfo_1.0.0               IRanges_2.44.0             
    #>  [9] S4Vectors_0.48.0            BiocGenerics_0.56.0        
    #> [11] generics_0.1.4              MatrixGenerics_1.22.0      
    #> [13] matrixStats_1.5.0           dominatR_0.99.5            
    #> [15] knitr_1.50                 
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] sass_0.4.10         SparseArray_1.10.1  lattice_0.22-7     
    #>  [4] digest_0.6.38       magrittr_2.0.4      RColorBrewer_1.1-3 
    #>  [7] evaluate_1.0.5      grid_4.5.2          fastmap_1.2.0      
    #> [10] jsonlite_2.0.0      Matrix_1.7-4        ggnewscale_0.5.2   
    #> [13] scales_1.4.0        tweenr_2.0.3        textshaping_1.0.4  
    #> [16] jquerylib_0.1.4     abind_1.4-8         cli_3.6.5          
    #> [19] rlang_1.1.6         polyclip_1.10-7     XVector_0.50.0     
    #> [22] withr_3.0.2         cachem_1.1.0        DelayedArray_0.36.0
    #> [25] yaml_2.3.10         S4Arrays_1.10.0     tools_4.5.2        
    #> [28] dplyr_1.1.4         vctrs_0.6.5         R6_2.6.1           
    #> [31] lifecycle_1.0.4     fs_1.6.6            MASS_7.3-65        
    #> [34] ragg_1.5.0          pkgconfig_2.0.3     desc_1.4.3         
    #> [37] pkgdown_2.2.0       bslib_0.9.0         pillar_1.11.1      
    #> [40] gtable_0.3.6        glue_1.8.0          ggforce_0.5.0      
    #> [43] systemfonts_1.3.1   xfun_0.54           tibble_3.3.0       
    #> [46] tidyselect_1.2.1    farver_2.1.2        htmltools_0.5.8.1  
    #> [49] rmarkdown_2.30      compiler_4.5.2      S7_0.2.1           
    #> [52] geomtextpath_0.2.0
