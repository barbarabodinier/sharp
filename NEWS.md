# sharp version 1.4.7

* Add reference to the publication in the Journal of Statistical Software

# sharp version 1.4.6

* Update the vignette

# sharp version 1.4.5

* Allow for alternative optimisation methods implemented in nloptr
* Update parallelisation, now using the future package
* Fix the formatting of continuous outcome in VariableSelection() 
* Update the vignette

# sharp version 1.4.4

* Update references with published articles

# sharp version 1.4.3

* Add sparse K means from the R package sparcl
* Allow for missing values in proportions for more flexibility

# sharp version 1.4.2

* Remove functions depending on regsem (removed from CRAN)
* Fix the use of packages in Suggests in the examples

# sharp version 1.4.1

* Add package vignette
* Use Ridge regression calibrated by cross validation instead of unpenalised regression in Refit(), ExplanatoryPerformance() and 
Incremental()
* Add new S3 class structural_model
* Fix inclusion of unpenalised predictors in Incremental()
* Fix clustering of rows in Clustering()

# sharp version 1.4.0

* Update the stability score used by default (n_cat=NULL), previous score can be used with n_cat=3
* Add new functions for structural equation modelling including StructuralModel(), PenalisedSEM(), PenalisedOpenMx(), 
PenalisedLinearSystem(), LavaanModel(), LavaanMatrix(), OpenMxModel(), OpenMxMatrix() and LinearSystemMatrix()
* Add new function CART() for classification and regression trees
* Add the option to run randomised or adaptive lasso in PenalisedRegression()
* Fix a bug when running multinomial lasso with predictors with null variance in the subsamples
* Fix a bug where additional parameters in ... were used in glm.control() within Refit()

# sharp version 1.3.0

* Add new functions for consensus clustering including Clustering(), Clusters(), ConsensusMatrix(), ClusteringPerformance() and more
* Add new print(), plot() and summary() functions
* Update plotting functions
* Fix parallelisation using argument n_cores in main functions
* Remove duplicated messages in ExplanatoryPerformance()
* Allow for factor ydata in VariableSelection() and related functions

# sharp version 1.2.1

* Update examples for use with fake 1.3.0
* Fix requirements on input data format in Refitting()
* Add resampling argument in Explanatory()
* Add optional beep at the end of the run in main functions
* Increase igraph vertex size in Graph() and plot()

# sharp version 1.2.0

* Add the functions Ensemble() and EnsemblePredictions() to build and predict from an ensemble model for VariableSelection()
* Add S3 classes including coef() and predict() for VariableSelection()
* Rename Recalibrate() as Refit()
* Fix use of CPSS in GraphicalModel() 
* Fix maximisation of the contrast
* Add simulation functions to the companion R package fake

# sharp version 1.1.0

First release of stability selection methods and simulation models.
