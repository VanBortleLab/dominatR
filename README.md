# <img src="man/figures/dominatR_header.png" width="1200px">

## Overview

**dominatR** is an R package designed for quantifying and visualizing feature dominance in datasets. dominatR leverages Shannon's information entropy to help identify features that are dominated by a particular category or condition. 

### Key Features

- **Built-in normalization methods:** If desired, feature counts can be normalized across conditions using built-in normalization functions (quantile normalization, min-max normalization, counts per million (cpm), reads per kb per million (rpkm), transcript per million (tpm), etc.)

- **Categorical entropy calculation:** Quickly identify features that are specific to a particular category or condition, or features that are relatively uniform across conditions using Shannon's entropy. 

- **Informative plots:** Generate informative and customizable plots that highlight feature dominance, aiding in the interpretation and communication of your data.

## Installation

You can install dominatR from GitHub using the `devtools` package:

```R
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install dominatR from GitHub
devtools::install_github("VanBortleLab/dominatR")
