test_that("The dimensions of outputs are consistent for BiSelection().", {
  skip_on_cran()
  set.seed(1)
  simul <- SimulateRegression(n = 100, pk = 15, q = 3, family = "gaussian")
  stab <- BiSelection(
    xdata = simul$xdata,
    ydata = simul$ydata,
    family = "gaussian",
    ncomp = 2,
    LambdaX = seq_len(ncol(simul$xdata) - 1),
    implementation = SparsePLS,
    K = 10,
    verbose = FALSE
  )
  expect_equal(class(stab), "bi_selection")
  expect_equal(ncol(simul$xdata), nrow(stab$selpropX))
  expect_equal(2, ncol(stab$selpropX))
  expect_equal(dim(stab$selpropX), dim(stab$selectedX))
  expect_equal(ncol(simul$ydata), nrow(stab$selpropY))
  expect_equal(2, ncol(stab$selpropY))
  expect_equal(dim(t(stab$coefX[, , 1])), dim(stab$selpropX))
  expect_equal(dim(t(stab$coefY[, , 1])), dim(stab$selpropY))
  expect_equal(stab$params$K, dim(stab$coefX)[3])
  expect_equal(stab$params$K, dim(stab$coefY)[3])
})
