test_that("The optimisation with nloptr is working for VariableSelection().", {
  skip_on_cran()
  set.seed(1)
  simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
  stab_grid <- VariableSelection(
    xdata = simul$xdata,
    ydata = simul$ydata,
    family = "gaussian",
    K = 10,
    optimisation = "grid_search",
    verbose = FALSE
  )
  stab_nloptr <- VariableSelection(
    xdata = simul$xdata,
    ydata = simul$ydata,
    family = "gaussian",
    K = 10,
    optimisation = "nloptr",
    verbose = FALSE
  )
  expect_lt(abs(max(stab_grid$S, na.rm=TRUE) - max(stab_nloptr$S, na.rm=TRUE)), 3)
  expect_lt(length(stab_nloptr$Lambda), 100)
})
