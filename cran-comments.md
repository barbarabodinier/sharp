## R CMD check results

There were no ERRORs or WARNINGs.

There were 3 NOTEs:

Maintainer: 'Barbara Bodinier <b.bodinier@imperial.ac.uk>'

Possibly misspelled words in DESCRIPTION:
  Monti (10:140)
  medoids (10:722)

* These words are not mispelled (surname and statistical term).

Suggests or Enhances not in mainstream repositories:
  rCOSA
Availability using Additional_repositories specification:
  rCOSA   yes   https://barbarabodinier.github.io/drat

* I have used the R package drat and hope this is ok (see Comments below).

Found the following (possibly) invalid URLs:
  URL: http://www.jstor.org/stable/2346178
    From: man/VariableSelection.Rd
    Status: 403
    Message: Forbidden
  URL: https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2004.02059.x
    From: man/Clustering.Rd
          man/DBSCANClustering.Rd
          man/HierarchicalClustering.Rd
          man/PAMClustering.Rd
    Status: 503
    Message: Service Unavailable

* These URLs have been checked and are correct.

Package suggested but not available for checking: 'rCOSA'

* I have used the R package drat and hope this is ok (see Comments below).

Package unavailable to check Rd xrefs: 'rCOSA'

* If that is ok, I would like to keep these references to the R package rCOSA which is under "Suggests" (see Comments below).


## Downstream dependencies

There are currently no downstream dependencies for this package.


## Comments

Some of the new functions may call a function in the R package rCOSA, not available on CRAN. I have used the R package drat to allow for use of rCOSA as stored in my Github Pages. I hope this is ok. 
