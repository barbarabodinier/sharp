
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

> In stability selection (N Meinshausen, P Bühlmann (2010)
> <doi:10.1111/j.1467-9868.2010.00740.x>) and consensus clustering (S
> Monti et al (2003) <doi:10.1023/A:1023949509487>), resampling
> techniques are used to enhance the reliability of the results. In this
> package, hyper-parameters are calibrated by maximising model
> stability, which is measured by the minus log-likelihood under the
> null hypothesis that all selection (or co-membership) probabilities
> are identical (B Bodinier et al (2021) \<arXiv:2106.02521\>).
> Functions are readily implemented for the use of LASSO regression,
> sparse PCA, sparse (group) PLS and graphical LASSO in stability
> selection, or hierarchical clustering, partitioning around medoids, K
> means, gaussian mixture models and DBSCAN in consensus clustering.
> Weighted distances from the COSA algorithm can be used in hierarchical
> clustering and partitioning around medoids.

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

## References

-   Barbara Bodinier, Sarah Filippi, Therese Haugdahl Nost, Julien
    Chiquet and Marc Chadeau-Hyam. Automated calibration for stability
    selection in penalised regression and graphical models: a
    multi-OMICs network application exploring the molecular response to
    tobacco smoking. (2021) arXiv.
    [link](https://doi.org/10.48550/arXiv.2106.02521)

-   Nicolai Meinshausen and Peter Bühlmann. Stability selection. (2010)
    Journal of the Royal Statistical Society: Series B (Statistical
    Methodology).
    [link](https://doi.org/10.1111/j.1467-9868.2010.00740.x)

-   Stefano Monti, Pablo Tamayo, Jill Mesirov and Todd Golub. (2003)
    Machine Learning. [link](https://doi.org/10.1023/A:1023949509487)
