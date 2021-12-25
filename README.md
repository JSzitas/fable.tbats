
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable.tbats

<!-- badges: start -->

[![R-CMD-check](https://github.com/JSzitas/fable.tbats/workflows/R-CMD-check/badge.svg)](https://github.com/JSzitas/fable.tbats/actions)
[![test-coverage](https://github.com/JSzitas/fable.tbats/workflows/test-coverage/badge.svg)](https://github.com/JSzitas/fable.tbats/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fable.tbats)](https://CRAN.R-project.org/package=fable.tbats)
<!-- badges: end -->

fable.tbats is a wrapper around the implementation of **tbats** from the
**forecast** package, to facilitate usage of tbats in the **fable**
framework. Note that this is by no means a reimplementation - and this
is currently not yet tested. Testing is a WIP.

Use at your own risk at this time!

## Installation

Currently only:

``` r
pak::pkg_install("JSzitas/fable.tbats")
```

# TODO:

Test against the forecast package (mainly for multiple seasonality
cases), or against the examples in fasster
