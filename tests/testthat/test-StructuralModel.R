test_that("The dimensions of outputs are consistent for StructuralModel().", {
  skip_on_cran()
  set.seed(1)
  pk <- c(3, 2, 3)
  simul <- SimulateStructural(
    n = 50,
    pk = pk,
    nu_between = 0.5,
    v_between = 1,
    v_sign = 1
  )
  dag <- LayeredDAG(layers = pk)
  stab <- StructuralModel(
    xdata = simul$data,
    adjacency = dag,
    K = 10,
    verbose = FALSE
  )
  expect_equal(class(stab), "structural_model")
  expect_equal(sum(dag), ncol(stab$selprop))
  expect_equal(length(stab$Lambda), nrow(stab$selprop))
  expect_equal(dim(stab$Beta)[1:2], dim(stab$selprop))
  expect_equal(stab$params$K, dim(stab$Beta)[3])
  expect_equal(length(stab$Lambda), nrow(stab$S))
})
