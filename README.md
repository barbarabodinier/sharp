
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sharp: Stability-enHanced Approaches using Resampling Procedures <img src="man/figures/logo.png" align="right" width="174" height="200"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sharp)](https://CRAN.R-project.org/package=sharp)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/sharp?color=blue)](https://r-pkg.org/pkg/sharp)
![GitHub last
commit](https://img.shields.io/github/last-commit/barbarabodinier/sharp?logo=GitHub&style=flat-square)
<!-- badges: end -->

## Description

> In stability selection, resampling techniques are used to enhance the
> reliability of the results. In this package, hyper-parameters are
> calibrated by maximising model stability, which is measured by the
> negative log-likelihood under the null hypothesis that all selection
> probabilities are identical. Functions are readily implemented for the
> use of LASSO regression, sparse PCA, sparse (group) PLS or graphical
> LASSO.

## Installation

The released version of the package can be installed from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sharp")
```

The development version can be installed from
[GitHub](https://github.com/):

``` r
remotes::install_github("barbarabodinier/sharp")
```

## Example datasets

To illustrate the use of the main functions implemented in **sharp**,
three artificial datasets are created:

``` r
library(sharp)

# Dataset for regression
set.seed(1)
data_reg <- SimulateRegression(n = 200, pk = 10)
x_reg <- data_reg$xdata
y_reg <- data_reg$ydata

# Dataset for graphical modelling
set.seed(1)
data_ggm <- SimulateGraphical(n = 200, pk = 20)
x_ggm <- data_ggm$data
```

Check out the R package
[**fake**](https://github.com/barbarabodinier/fake) for more details on
these data simulation models.

## Main functions

### Variable selection

In a regression context, stability selection is done using LASSO
regression as implemented in the R package
[**glmnet**](https://CRAN.R-project.org/package=glmnet).

``` r
stab_reg <- VariableSelection(xdata = x_reg, ydata = y_reg)
SelectedVariables(stab_reg)
```

### Graphical modelling

In a graphical modelling context, stability selection is done using the
graphical LASSO as implemented in the R package
[**glassoFast**](https://CRAN.R-project.org/package=glassoFast).

``` r
stab_ggm <- GraphicalModel(xdata = x_ggm)
Adjacency(stab_ggm)
```

## Extraction and visualisation of the results

It is strongly recommended to check the calibration of the
hyper-parameters using the function `CalibrationPlot()` on the output
from any of the main functions listed above. The functions `print()`,
`summary()` and `plot()` can also be used on the outputs from the main
functions.

## References

- Barbara Bodinier, Sarah Filippi, Therese Haugdahl Nost, Julien Chiquet
  and Marc Chadeau-Hyam. Automated calibration for stability selection
  in penalised regression and graphical models: a multi-OMICs network
  application exploring the molecular response to tobacco
  smoking. (2021) arXiv.
  [link](https://doi.org/10.48550/arXiv.2106.02521)

- Nicolai Meinshausen and Peter Bühlmann. Stability selection. (2010)
  Journal of the Royal Statistical Society: Series B (Statistical
  Methodology). [link](https://doi.org/10.1111/j.1467-9868.2010.00740.x)
