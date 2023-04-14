
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sharp: Stability-enHanced Approaches using Resampling Procedures <img src="man/figures/logo.png" align="right" width="174" height="200"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sharp)](https://CRAN.R-project.org/package=sharp)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/sharp?color=blue)](https://r-pkg.org/pkg/sharp)
![GitHub last
commit](https://img.shields.io/github/last-commit/barbarabodinier/sharp?logo=GitHub&style=flat-square)
<!-- badges: end -->

## Description

> In stability selection and consensus clustering, resampling techniques
> are used to enhance the reliability of the results. In this package,
> hyper-parameters are calibrated by maximising model stability, which
> is measured under the null hypothesis that all selection (or
> co-membership) probabilities are identical. Functions are readily
> implemented for the use of LASSO regression, sparse PCA, sparse
> (group) PLS or graphical LASSO in stability selection, and
> hierarchical clustering, partitioning around medoids, K means or
> Gaussian mixture models in consensus clustering.

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

To illustrate the use of the main functions implemented in
[**sharp**](https://CRAN.R-project.org/package=sharp), three artificial
datasets are created:

``` r
library(sharp)

# Dataset for regression
set.seed(1)
data_reg <- SimulateRegression(n = 200, pk = 10)
x_reg <- data_reg$xdata
y_reg <- data_reg$ydata

# Dataset for structural equation modelling
set.seed(1)
data_sem <- SimulateStructural(n = 200, pk = c(5, 2, 3))
x_sem <- data_sem$data

# Dataset for graphical modelling
set.seed(1)
data_ggm <- SimulateGraphical(n = 200, pk = 20)
x_ggm <- data_ggm$data

# Dataset for clustering
set.seed(1)
data_clust <- SimulateClustering(n = c(10, 10, 10))
x_clust <- data_clust$data
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

### Structural equation modelling

In a structural equation modelling context, stability selection is done
using series of LASSO regressions as implemented in the R package
[**glmnet**](https://CRAN.R-project.org/package=glmnet).

``` r
dag <- LayeredDAG(layers = c(5, 2, 3))
stab_sem <- StructuralEquations(xdata = x_sem, adjacency = dag)
LinearSystemMatrix(vect = Stable(stab_sem), adjacency = dag)
```

### Graphical modelling

In a graphical modelling context, stability selection is done using the
graphical LASSO as implemented in the R package
[**glassoFast**](https://CRAN.R-project.org/package=glassoFast).

``` r
stab_ggm <- GraphicalModel(xdata = x_ggm)
Adjacency(stab_ggm)
```

### Clustering

Consensus clustering is done using hierarchical clustering as
implemented in the R package
[**stats**](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html).

``` r
stab_clust <- Clustering(xdata = x_clust)
Clusters(stab_clust)
```

## Extraction and visualisation of the results

It is strongly recommended to check the calibration of the
hyper-parameters using the function `CalibrationPlot()` on the output
from any of the main functions listed above. The functions `print()`,
`summary()` and `plot()` can also be used on the outputs from the main
functions.

## Parametrisation

Stability selection and consensus clustering can theoretically be done
by aggregating the results from any selection (or clustering) algorithm
on subsamples of the data. The choice of the underlying algorithm to use
is specified in argument `implementation` in the main functions.
Consensus clustering using partitioning around medoids, K means or
Gaussian mixture models are also supported in
[**sharp**](https://CRAN.R-project.org/package=sharp):

``` r
stab_clust <- Clustering(xdata = x_clust, implementation = PAMClustering)
stab_clust <- Clustering(xdata = x_clust, implementation = KMeansClustering)
stab_clust <- Clustering(xdata = x_clust, implementation = GMMClustering)
```

Other algorithms can be used by defining a wrapper function to be called
in `implementation`. Check out the documentation of `GraphicalModel()`
for an example using a shrunk estimate of the partial correlation
instead of the graphical LASSO.

## References

- Barbara Bodinier, Sarah Filippi, Therese Haugdahl Nost, Julien Chiquet
  and Marc Chadeau-Hyam. Automated calibration for stability selection
  in penalised regression and graphical models. (2021) arXiv.
  [link](https://doi.org/10.48550/arXiv.2106.02521)

- Nicolai Meinshausen and Peter Bühlmann. Stability selection. (2010)
  Journal of the Royal Statistical Society: Series B (Statistical
  Methodology). [link](https://doi.org/10.1111/j.1467-9868.2010.00740.x)

- Stefano Monti, Pablo Tamayo, Jill Mesirov and Todd Golub. (2003)
  Machine Learning. [link](https://doi.org/10.1023/A:1023949509487)
