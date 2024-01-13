test_that("The dimensions of outputs are consistent for Clustering().", {
  skip_on_cran()
  set.seed(1)
  simul <- SimulateClustering(
    n = c(30, 30, 30),
    nu_xc = 1,
    ev_xc = 0.5,
  )
  stab <- Clustering(
    xdata = simul$data,
    K = 100,
    verbose = FALSE
  )
  expect_equal(class(stab), "clustering")
  expect_equal(nrow(stab$coprop), ncol(stab$coprop))
  expect_equal(nrow(simul$data), ncol(stab$coprop))
  expect_equal(nrow(stab$Lambda), nrow(stab$nc))
  expect_equal(nrow(stab$Lambda), dim(stab$coprop)[3])
  expect_equal(nrow(stab$Lambda), nrow(stab$Sc))
})
