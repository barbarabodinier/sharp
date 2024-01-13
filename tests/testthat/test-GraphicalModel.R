test_that("The dimensions of outputs are consistent for GraphicalModel().", {
  skip_on_cran()
  set.seed(1)
  simul <- SimulateGraphical(n = 100, pk = 20, nu_within = 0.1)
  stab <- GraphicalModel(
    xdata = simul$data,
    K = 10,
    verbose = FALSE
  )
  expect_equal(class(stab), "graphical_model")
  expect_equal(nrow(stab$selprop), ncol(stab$selprop))
  expect_equal(ncol(simul$data), ncol(stab$selprop))
  expect_equal(length(stab$Lambda), dim(stab$selprop)[3])
  expect_equal(length(stab$Lambda), nrow(stab$S))
})
