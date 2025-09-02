<img src="man/figures/dominatR_header.png" width="1200px"/>

## Overview

**dominatR** is an R package for quantifying and visualizing feature dominance in datasets. dominatR makes use of Shannon's entropy to identify features that are dominated within a particular category or condition.

### Features

-   **Built-in normalization methods:** If desired, feature counts can be normalized across conditions using built-in functions for quantile normalization, min-max normalization, counts per million (cpm), reads per kb per million (rpkm), transcript per million (tpm), and others.

-   **Categorical entropy calculation:** Quickly identify features that are specific to a particular category or condition, as well as features that are relatively uniform across conditions

-   **Feature dominance plots:** Generate informative and customizable plots that highlight feature dominance, aiding data interpretation and communication.

## Installation

dominatR can be installed from GitHub using the `devtools` package:

``` r
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install dominatR from GitHub
devtools::install_github("VanBortleLab/dominatR")
```

## Usage
