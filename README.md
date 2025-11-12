README
================

# m6APrediction: m6A Modification Site Prediction in R

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

The `m6APrediction` package provides a straightforward and efficient
tool to predict N6-methyladenosine (m6A) modification sites from various
sequence and genomic features. It uses a pre-trained Random Forest model
to deliver predictions for both single data points and batch datasets,
making it suitable for a wide range of transcriptomic analyses.

The core purpose of this package is to allow researchers to quickly
assess the probabilities of m6A modification on RNA sequences based on a
set of predictive features.

## Model Performance

The underlying Random Forest model demonstrates strong predictive
performance. The Receiver Operating Characteristic (ROC) and
Precision-Recall (PRC) curves below show its ability to distinguish
between true m6A sites and non-m6A sites.

|           ROC Curve and PRC Curve            |
|:--------------------------------------------:|
| ![](man/figures/roc_curve_and_prc_curve.png) |

## Installation

You can install the development version of `m6APrediction` from GitHub
with:

``` r
# install.packages("devtools")
devtools::install_github("Chenyu12161004/m6APrediction")
```

## Usage Example

Here is a minimal example demonstrating how to use the two main
functions of the package: `prediction_single()` and
`prediction_multiple()`.

First, load the package:

``` r
devtools::load_all()
```

    ## ℹ Loading m6APrediction

### Single Prediction

To predict the m6A status for a single site with specific features:

``` r
# Load the pre-trained model included with the package
rf_model <- readRDS(
  system.file("extdata", "rf_fit.rds", package = "m6APrediction")
)

# Predict for a single data point
single_pred <- prediction_single(
  ml_fit = rf_model,
  gc_content = 0.5,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 10,
  distance_to_junction = 8,
  evolutionary_conservation = 0.5,
  DNA_5mer = "GGACA"
)

# Print the result
print(single_pred)
```

    ##   predicted_m6A_prob predicted_m6A_status 
    ##               "0.94"           "Positive"

### Multiple Predictions

To run predictions on a batch of sites from a data frame:

``` r
# Load the example data frame included with the package
example_df <- read.csv(
  system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
)

# Run batch prediction
# Note: We are using the rf_model loaded in the previous step
multi_pred <- prediction_multiple(
  ml_fit = rf_model,
  feature_df = example_df,
  positive_threshold = 0.6
)

# Display the first few rows of the results
head(multi_pred)
```

    ##   gc_content RNA_type RNA_region exon_length distance_to_junction
    ## 1  0.6616915   lncRNA     intron    0.000000            10.668885
    ## 2  0.5223881     mRNA      3'UTR   10.790348             9.330917
    ## 3  0.3980100     mRNA      5'UTR    9.130571             7.592457
    ## 4  0.5870647     mRNA        CDS    7.426265             4.700440
    ## 5  0.6218905     mRNA        CDS    8.445015             5.000000
    ## 6  0.4179104     mRNA      3'UTR   11.477252            11.039605
    ##   evolutionary_conservation DNA_5mer nt_pos1 nt_pos2 nt_pos3 nt_pos4 nt_pos5
    ## 1                0.01641791    GGACC       G       G       A       C       C
    ## 2                0.02537313    TGACC       T       G       A       C       C
    ## 3                0.61592040    GAACA       G       A       A       C       A
    ## 4                0.42238806    GAACC       G       A       A       C       C
    ## 5                0.01094527    GAACT       G       A       A       C       T
    ## 6                0.74875622    AAACA       A       A       A       C       A
    ##   predicted_m6A_prob predicted_m6A_status
    ## 1              0.630             Positive
    ## 2              0.008             Negative
    ## 3              0.146             Negative
    ## 4              0.002             Negative
    ## 5              0.686             Positive
    ## 6              0.010             Negative
