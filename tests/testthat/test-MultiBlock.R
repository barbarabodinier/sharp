test_that("multi-block grid is working", {
  # Multi-block
  Lambda <- matrix(c(
    0.8, 0.6, 0.3,
    0.5, 0.4, 0.2,
    0.7, 0.5, 0.1
  ),
  ncol = 3, byrow = TRUE
  )
  mygrid <- BlockLambdaGrid(Lambda, lambda_other_blocks = 0.1)
  expect_equal(dim(mygrid$Lambda), c(3 * nrow(Lambda), 3))

  # Multi-parameter
  Lambda <- matrix(c(
    0.8, 0.6, 0.3,
    0.5, 0.4, 0.2,
    0.7, 0.5, 0.1
  ),
  ncol = 3, byrow = TRUE
  )
  mygrid <- BlockLambdaGrid(Lambda, lambda_other_blocks = NULL)
  expect_equal(Lambda, mygrid$Lambda)
})
