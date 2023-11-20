# <img src="man/figures/dominatR_header.png" width="1200px">

## Overview

**dominatR** is an R package designed for visualizing feature dominance and quantifying feature uniformity in datasets. dominatR leverages Shannon's information entropy to help identify which features are dominated by a particular category or condition. 

### Key Features

- **Built-in normalization methods:** If desired, feature counts can be normalized across conditions using built-in normalization functions (quantile normalization, min-max normalization, counts per million (cpm), reads per kb per million (rpkm), transcript per million (tpm), etc.)

- **Categorical entropy calculation:** Use Shannon's entropy to measure feature dominance, allowing you to quickly identify features that are specific to a particular category or condition, or features that relatively uniform across conditions

- **Informative plots:** Generate informative plots that highlight feature dominance and uniformity, aiding in the interpretation of your data.
