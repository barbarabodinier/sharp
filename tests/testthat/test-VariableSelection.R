test_that("The dimensions of outputs are consistent for VariableSelection().", {
  skip_on_cran()
  set.seed(1)
  simul <- SimulateRegression(n = 100, pk = 50, family = "gaussian")
  stab <- VariableSelection(
    xdata = simul$xdata,
    ydata = simul$ydata,
    family = "gaussian",
    K = 10,
    verbose = FALSE
  )
  expect_equal(class(stab), "variable_selection")
  expect_equal(ncol(simul$xdata), ncol(stab$selprop))
  expect_equal(length(stab$Lambda), nrow(stab$selprop))
  expect_equal(dim(stab$Beta)[1:2], dim(stab$selprop))
  expect_equal(stab$params$K, dim(stab$Beta)[3])
  expect_equal(length(stab$Lambda), nrow(stab$S))
})
