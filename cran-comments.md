## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1111/j.1467-9868.2010.00740.x
    From: man/BiSelection.Rd
          man/GraphicalModel.Rd
          man/PFER.Rd
          man/StabilityMetrics.Rd
          man/VariableSelection.Rd
          man/sharp-package.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1111/j.1467-9868.2011.01034.x
    From: man/BiSelection.Rd
          man/GraphicalModel.Rd
          man/PFER.Rd
          man/StabilityMetrics.Rd
          man/VariableSelection.Rd
          man/sharp-package.Rd
    Status: 503
    Message: Service Unavailable

These URLs can be used.


## Downstream dependencies

There are currently no downstream dependencies for this package.


## Comments

I have split the R package sharp (version 1.1.0) into two packages: sharp (this submission, version 1.2.0) and fake (version 1.0.0).

The R package fake (version 1.0.0) is made of functions from sharp (version 1.1.0). 

The R package fake (version 1.0.0) has been released on CRAN earlier this week (accepted on 9 August 2022).

All functions included in fake (version 1.0.0) have been removed from sharp (version 1.2.0).

The R package sharp (version 1.2.0) depends on fake (version 1.0.0).

I hope this is ok.

Many thanks.
 

