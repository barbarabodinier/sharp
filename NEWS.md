# sharp version 1.4.1

* Added package vignette
* Used Ridge regression calibrated by cross validation instead of unpenalised regression in Refit(), ExplanatoryPerformance() and 
Incremental()
* Added new S3 class structural_model
* Fixed inclusion of unpenalised predictors in Incremental()
* Fixed clustering of rows in Clustering()

# sharp version 1.4.0

* Updated the stability score used by default (n_cat=NULL), previous score can be used with n_cat=3
* Added new functions for structural equation modelling including StructuralModel(), PenalisedSEM(), PenalisedOpenMx(), 
PenalisedLinearSystem(), LavaanModel(), LavaanMatrix(), OpenMxModel(), OpenMxMatrix() and LinearSystemMatrix()
* Added new function CART() for classification and regression trees
* Added the option to run randomised or adaptive lasso in PenalisedRegression()
* Fixed a bug when running multinomial lasso with predictors with null variance in the subsamples
* Fixed a bug where additional parameters in ... were used in glm.control() within Refit()

# sharp version 1.3.0

* Added new functions for consensus clustering including Clustering(), Clusters(), ConsensusMatrix(), ClusteringPerformance() and more
* Added new print(), plot() and summary() functions
* Updated plotting functions
* Fixed parallelisation using argument n_cores in main functions
* Remove duplicated messages in ExplanatoryPerformance()
* Allow for factor ydata in VariableSelection() and related functions

# sharp version 1.2.1

* Updated examples for use with fake 1.3.0
* Fixed requirements on input data format in Refitting()
* Added resampling argument in Explanatory()
* Added optional beep at the end of the run in main functions
* Increased igraph vertex size in Graph() and plot()

# sharp version 1.2.0

* Added the functions Ensemble() and EnsemblePredictions() to build and predict from an ensemble model for VariableSelection()
* Added S3 classes including coef() and predict() for VariableSelection()
* Renamed Recalibrate() as Refit()
* Fixed use of CPSS in GraphicalModel() 
* Fixed maximisation of the contrast
* Simulation functions now added to the companion R package fake

# sharp version 1.1.0

First release of stability selection methods and simulation models.
