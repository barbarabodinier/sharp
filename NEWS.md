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
