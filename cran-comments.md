## R CMD check results

There were no ERRORs or WARNINGs.

There were 3 NOTEs:

Suggests or Enhances not in mainstream repositories:
  rCOSA
Availability using Additional_repositories specification:
  rCOSA   yes   https://barbarabodinier.github.io/drat

* I have used the R package drat. I also had this note with the previous version of sharp (1.4.1).

Found the following (possibly) invalid URLs:
  URL: http://www.jstor.org/stable/2346178
    From: man/PenalisedRegression.Rd
          man/VariableSelection.Rd
    Status: 403
    Message: Forbidden
  URL: https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2004.02059.x
    From: man/Clustering.Rd
          man/DBSCANClustering.Rd
          man/HierarchicalClustering.Rd
          man/PAMClustering.Rd
    Status: 403
    Message: Forbidden

* These URLs have been checked and are correct. I also had this note with the previous version of sharp (1.4.2). 

Package suggested but not available for checking: 'rCOSA'

* I have used the R package drat. I also had this note with the previous version of sharp (1.4.1).

Package unavailable to check Rd xrefs: 'rCOSA'

* I would like to keep these references to the R package rCOSA which is under "Suggests". I also had this note with the 
previous version of sharp (1.4.3).


## Downstream dependencies

There are currently no downstream dependencies for this package.
