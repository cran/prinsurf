
<!-- README.md is generated from README.Rmd. Please edit that file -->

# prinsurf

<!-- badges: start -->
<!-- badges: end -->

The goal of prinsurf is to construct a principal surface that is
two-dimensional and passes through the middle of a $p$-dimensional
dataset.

## Installation

You can install the development version of prinsurf from
[GitHub](https://github.com/) with:

``` r
library(devtools)
devtools::install_github("RaeesaGaney91/prinsurf")
```

## Example

This is a basic example on a simulated data set:

``` r
library(prinsurf)
surface <- principal.surface(X)
```

<img src="man/figures/README-example-1.png" width="100%" />

<figure>
<img src="man/figures/3d_plot.png" width="800" alt="3D Plot" />
<figcaption aria-hidden="true">3D Plot</figcaption>
</figure>

## Report Bugs and Support

If you encounter any issues or have questions, please open an issue on
the GitHub repository.
