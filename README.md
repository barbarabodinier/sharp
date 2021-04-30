
<!-- README.md is generated from README.Rmd. Please edit that file -->

# focus: Feature selectiOn and Clustering Using Stability

<!-- badges: start -->

<!-- badges: end -->

Implementation of stability-enhanced models for variable selection in
regression, graphical modelling and clustering. These methods are based
on resampling approaches to compute selection proportions. Calibration
of the models is done via maximisation of a stability score measuring
how unlikely it is that the selection procedure is uniform.

## Installation

<!-- You can install the released version of focus from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("focus") -->

<!-- ``` -->

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("barbarabodinier/focus")
```

## Variable selection

This is a basic example which shows you how to solve a common problem:

``` r
library(focus)
## basic example code
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
