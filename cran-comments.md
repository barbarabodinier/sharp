## Resubmission

This is a resubmission. In this version I have: 

* Added references in the Description field of the DESCRIPTION file.

* Added missing \value to GraphComparison.Rd and SelectionPerformanceGraph.Rd.

* Removed examples for unexported functions. Affected functions are: Coefficients() and HugeAdjacency(). 

* Removed or replaced by \donttest all instances of \dontrun.

* Ensured that user's par are reset to original values within functions with on.exit(). Affected functions are: Heatmap() and CalibrationPlot(). 

* Ensured that user's par are reset to original values at the end of examples. Affected functions are: BiSelection(), CalibrationPlot(), GraphicalModel(), Heatmap(), PLS(), SimulateRegression() and VariableSelection(). 


## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

Possibly misspelled words in DESCRIPTION:
  Bodinier (9:446)
  Chadeau (9:489)
  Chiquet (9:478)
  Filippi (9:458)
  Hyam (9:497)
  Meinshausen (9:224)
  Nost (9:470)
  dimensionality (9:117)
  enHanced (3:18)
  hlmann (9:242)

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

Found the following (possibly) invalid DOIs:
  DOI: 10.1111/j.1467-9868.2010.00740.x
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
    
The potentially misspelled words are family names (Bodinier, Chadeau, Chiquet, Filippi, Meinshausen, Buhlmann) or statistical techniques (dimensionality from dimensionality reduction). If acceptable, I would like to keep the word "enHanced" written this way to emphasise on the meaning of the package name: Stability enHanced Approaches using Resampling Procedures (SHARP). 