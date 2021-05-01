test_that("output of SimulateGraphical() is of correct dimension", {
  n=78
  pk=12
  nlambda=3
  K=5
  pi_list=seq(0.6,0.7,length.out=15)
  simul=SimulateRegression(n=n, pk=pk, family="gaussian")
  stab=VariableSelection(xdata=simul$X, ydata=simul$Y,
                         Lambda_cardinal=nlambda, K=K,
                         pi_list=pi_list,
                         verbose=FALSE)

  ### Checking dimensions of the outputs
  # Group of outputs 1
  expect_equal(dim(stab$S), c(nlambda,1))
  expect_equal(dim(stab$Lambda), c(nlambda,1))
  expect_equal(dim(stab$Q), c(nlambda,1))
  expect_equal(dim(stab$Q_s), c(nlambda,1))
  expect_equal(dim(stab$P), c(nlambda,1))
  expect_equal(dim(stab$PFER), c(nlambda,1))
  expect_equal(dim(stab$FDP), c(nlambda,1))

  # Group of outputs 2
  expect_equal(dim(stab$S_2d), c(nlambda,length(pi_list)))
  expect_equal(dim(stab$PFER_2d), c(nlambda,length(pi_list)))
  expect_equal(dim(stab$FDP_2d), c(nlambda,length(pi_list)))

  # Group of outputs 3
  expect_equal(dim(stab$selprop), c(nlambda,pk))
  expect_equal(dim(stab$Beta), c(nlambda,pk,))

})
