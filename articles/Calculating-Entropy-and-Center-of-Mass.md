# Calculating Entropy and Center of Mass

``` r
library(dominatR)
library(airway)
library(dominatRData)
library(ggplot2)
```

## Data

The data used in this article comes from two different sources:

1.  A `SummarizedExperiment` object from the library `airway`
2.  A `data.frame` retrieved from the supplementary library
    `dominatRData`

The purpose of using two different sources is meant to show the capacity
of dominatR to process different objects.

``` r
### summarized experiment
data("airway")
se <-airway

# Only use a random subset of 500 rows
set.seed(123)
idx <- sample(seq_len(nrow(se)), size = min(500, nrow(se)))
se <- se[idx, ]


# dataframe
data("rnapol_score")
df1 = rnapol_score

data("atac_tissue_score")
df2 = atac_tissue_score
```

In the context of information theory, **entropy (Shannon Entropy)** is a
metric used to measure the uncertainty associated with a set of
variables. It’s value ranges from 0 to 𝑙𝑜𝑔(𝑁), where N is the total
number of variables accounted for. Its interpretation is
straightforward:

1.  An observation with high entropy values is not uncertain, that
    observation is expected to occur at a high probability level, across
    all the variables. For example, Gene A if it is expressed in all the
    variables at the same levels will have a high entropy value.

2.  An observation with low entropy values is uncertain, that
    observation is expected to occur at a high probability level, only
    in a set of few variables. For example, Gene B if it is expressed
    only in variable A or variable B and not in the rest of tissues will
    have a low entropy value.

Entropy is calculated by using the following formula:

$H(X) = - \sum p(x)*log2\left( p(x) \right)$

Where $p(x)$ is equal to the relative levels of a gene (g) in a specific
variable (t). From this formula, it is possible to measure the
categorical tissue specificity that is defined as:

$Q_{(g|t)} = H(X) - log2\left( p(x) \right)$

This formula specifies a ‘domination’ value for a specific gene in a
specific tissue. Therefore low values of $Q_{(g|t)}$ occurs when a gene
or observation is mostly present in a small subset of tissues.

In physics, **center of mass** is a way to summarize how a quantity
(typically mass) is distributed across space. It represents the unique
point where the weighted average of all positions, relative to their
associated weights, balances out. For weights $w_{i} \geq 0$ at
coordinates $r_{i}$​, the CoM is

$CoM = \frac{\sum_{i}\hspace{0pt}w_{i}\hspace{0pt}r_{i}}{\sum_{i}\hspace{0pt}w_{i}}\hspace{0pt}\hspace{0pt}\hspace{0pt}$.

For the context of this package, the weights are sample-level
measurements (e.g gene expression) and the positions are user defined
coordinates assigned to samples (For up to three dimensions). CoM
summarizes where the observation is positioned in the coordinate system.
Its interpretation is as follows:

1.  An observation that is dominated by one sample, will be located near
    that sample’s coordinate.
2.  An observation that is uniformly distributed across samples, will be
    located near the geometric center of all sample coordinates.

The package provides formulas that are useful to calculate the Entropy
$H(X)$, ‘Categorical’ Entropy $Q_{(g|t)}$ and Center of Mass $CoM$:

- [`entropy()`](https://vanbortlelab.github.io/dominatR/reference/entropy.md)

- [`Qentropy()`](https://vanbortlelab.github.io/dominatR/reference/Qentropy.md)

- [`centmass()`](https://vanbortlelab.github.io/dominatR/reference/centmass.md)  

### Center of Mass

Center of mass is useful when using no more than 3 dimensions given that
more than 3 dimensions are not representable in a bidimensional plane.
The function requires specification of a set of coordinates, by default,
the function uses `x_coord = c(0, 1, 0.5)` and
`y_coord = c(0, 0, sqrt(3)/2)`, which correspond to an equilateral
triangle. The function returns the pair of coordinates x and y for each
observation representing the spatial location of the CoM as a data frame
or as extra columns in
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
is a SummarizedExperiment is used.

#### Example 1: A data.frame

``` r
df = df1
## calculating center of mass in the numerical functions
df = centmass(df[,5:7])
head(df)
#>        comx       comy
#> 1 0.1022727 0.05248639
#> 2 0.4998737 0.86580669
#> 3 0.1439394 0.00000000
#> 4 0.5002817 0.86510373
#> 5 0.4962568 0.85851978
#> 6 0.4982841 0.85948712
```

#### Example 2: A summarized Experiment

``` r
# Subsetting 3 columns
se2 = se[,1:3]
se2 <- centmass(se2)

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

### Entropy

Entropy offers the benefit of calculating the levels of surprise of an
observation across N-Variables. The function transforms a data frame
into proportions and adds a column with the name Entropy. Nevertheless,
if using a `SummarizedExperiment` object, there is an option for the
user to store the proportions in a new assay and adds a column with the
name Entropy in the `rowdata()`.

#### Example 1: A data.frame

``` r
df = df1
df = entropy(df)
head(df)
#>    Chr     Start      Stop RNA_type         pol1         pol2         pol3
#> 1 chr1 0.4999965 0.5000034     tRNA 7.835873e-08 6.501380e-09 5.474846e-09
#> 2 chr1 0.4999952 0.5000012     tRNA 9.079536e-10 0.000000e+00 3.594286e-06
#> 3 chr1 0.4999970 0.5000030     tRNA 3.417090e-08 5.745550e-09 0.000000e+00
#> 4 chr1 0.4999921 0.4999982     tRNA 2.419044e-09 7.861892e-09 9.649868e-06
#> 5 chr1 0.4999921 0.4999982     tRNA 7.856825e-08 5.741526e-09 9.643649e-06
#> 6 chr1 0.4999959 0.5000019     tRNA 1.208625e-08 4.532346e-09 2.184591e-06
#>    Entropy
#> 1 1.000002
#> 2 1.000067
#> 3 1.000001
#> 4 1.000165
#> 5 1.000167
#> 6 1.000042
```

#### Example 2: A summarized Experiment

``` r
# Subsetting 3 columns
se2 = se
se2 <- entropy(se2, new_assay_name = 'Entropy')

head(rowData(se2))
#> DataFrame with 6 rows and 11 columns
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
#>                 seq_coord_system         symbol   Entropy
#>                        <integer>    <character> <numeric>
#> ENSG00000260166               NA  RP11-863P13.6  0.000000
#> ENSG00000266931               NA RP11-1252D15.1  0.000000
#> ENSG00000104774               NA         MAN2B1  2.979167
#> ENSG00000267583               NA  RP11-322E11.5  0.000000
#> ENSG00000227581               NA   RP13-140E4.1  0.918296
#> ENSG00000227317               NA          DDAH2  0.000000
```

### Categorical Entropy

For calculating Categorical Entropy, it is required to compute first
Entropy in the dataset of interest. If the object is a dataframe, the
Entropy transformed dataframe will be subsequently transformed. If the
object is a SummarizedExperiment, the dataset can be stored in a new
assay.

#### Example 1: A data.frame

``` r
df = df2 
df = entropy(df)
df = Qentropy(df)
head(df[order(df$Heart),])
#>        core_type   Chr     Start       End
#> 9363    Specific  chr1 151022749 151022852
#> 1112    Specific  chr4 155464005 155465338
#> 2850    Specific  chr1 237945359 237946320
#> 8751    Specific chr19  45427468  45427615
#> 4324 25% Tissues  chr9  33657848  33659206
#> 9799 25% Tissues chr21   8253841   8255704
#>                                              Gene  Index     Type    Heart
#> 9363                 U6_snRNA_URS00006B5505&&9098   9098   u6 RNA 2.000000
#> 1112 LSU_rRNA_bacteria_rRNA_URS0000994EE4&&122489 122489 LSU rRNA 3.169925
#> 2850  SSU_rRNA_bacteria_rRNA_URS000096E89F&&14540  14540 SSU rRNA 4.000000
#> 8751     Metazoa_SRP_SRP_RNA_URS00009442C9&&78677  78677      SRP 4.643856
#> 4324 LSU_rRNA_bacteria_rRNA_URS00009A7234&&163885 163885 LSU rRNA 5.169925
#> 9799   SSU_rRNA_eukarya_rRNA_URS0000ABD7D5&&99559  99559 SSU rRNA 5.169925
#>      Spleen    Liver    Colon Adrenal.Gland Lung Pancreas Gallbladder
#> 9363    Inf      Inf      Inf           Inf  Inf      Inf         Inf
#> 1112    Inf 3.169925      Inf      3.169925  Inf      Inf         Inf
#> 2850    Inf      Inf 4.000000           Inf  Inf      Inf         Inf
#> 8751    Inf      Inf      Inf           Inf  Inf      Inf         Inf
#> 4324    Inf      Inf 5.169925      5.169925  Inf      Inf         Inf
#> 9799    Inf      Inf      Inf           Inf  Inf 5.169925         Inf
#>      Urinary.Tract   Breast Fallopian.Tube Psoas.Muscle Vena.Cava
#> 9363           Inf      Inf            Inf     2.000000       Inf
#> 1112           Inf      Inf            Inf          Inf       Inf
#> 2850           Inf 4.000000            Inf          Inf       Inf
#> 8751           Inf      Inf            Inf          Inf  4.643856
#> 4324      5.169925 5.169925            Inf          Inf       Inf
#> 9799           Inf      Inf            Inf     5.169925  5.169925
#>      Gastroesophageal.Sphincter Adipose.Tissue Sciatic.Nerve Ovary...Uterus
#> 9363                        Inf            Inf           Inf            Inf
#> 1112                        Inf            Inf           Inf            Inf
#> 2850                        Inf            Inf           Inf            Inf
#> 8751                        Inf       4.643856      4.643856       4.643856
#> 4324                        Inf       5.169925           Inf            Inf
#> 9799                        Inf            Inf      5.169925            Inf
#>      Stomach    Brain
#> 9363     Inf      Inf
#> 1112     Inf      Inf
#> 2850       4      Inf
#> 8751     Inf      Inf
#> 4324     Inf      Inf
#> 9799     Inf 5.169925
```

From this dataframe, it is possible to see how the observations
`gene_9363`, `gene_1112`, `gene_2850` have low values of Qentropy
emphasizing high score or normalized features in `Heart` when compared
to the other samples.

#### Example 2: A summarized Experiment

``` r
# Subsetting 3 columns
se2 = se
se2 <- entropy(se2, new_assay_name = 'Entropy')
se2 <- Qentropy(se2, assay_name = 'Entropy', new_assay_name = 'Qentropy')

head(assay(se2, 'Qentropy'))
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
    #> [46] tidyselect_1.2.1    knitr_1.50          farver_2.1.2       
    #> [49] htmltools_0.5.8.1   rmarkdown_2.30      compiler_4.5.2     
    #> [52] S7_0.2.1            geomtextpath_0.2.0
