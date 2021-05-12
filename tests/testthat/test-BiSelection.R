test_that("BiSelection() models are working", {
  # Data simulation
  K <- 5
  pk <- 15
  set.seed(1)
  simul <- SimulateRegression(n = 50, pk = pk, family = "gaussian")
  ydata <- cbind(simul$Y, matrix(rnorm(50 * 3), ncol = 3))
  colnames(ydata) <- paste0("outcome", 1:4)
  x <- simul$X
  y <- ydata

  # sPLS: sparsity on both X and Y
  stab <- BiSelection(
    xdata = x, ydata = y,
    family = "gaussian", K = K, ncomp = 2,
    LambdaX = 1:2,
    LambdaY = 1:2,
    implementation = SparsePLS
  )

  # sgPLS: sparsity on both X and Y
  stab <- BiSelection(
    xdata = x, ydata = y,
    group_x = c(10, 5), group_y = c(1, 3),
    family = "gaussian", K = K, ncomp = 2,
    LambdaX = 1:2, AlphaX = c(0.1, 0.3),
    LambdaY = 1:2, AlphaY = c(0.1, 0.3),
    implementation = SparseGroupPLS
  )

  ## Checking with other resampling
  # Data simulation
  simul <- SimulateRegression(family = "binomial")

  # Data simulation for a binary confounder
  conf <- ifelse(runif(n = 100) > 0.5, yes = 1, no = 0)

  # User-defined resampling function
  BalancedResampling <- function(data, tau, Z, ...) {
    s <- NULL
    for (z in unique(Z)) {
      s <- c(s, sample(which((data == "0") & (Z == z)), size = tau * sum((data == "0") & (Z == z))))
      s <- c(s, sample(which((data == "1") & (Z == z)), size = tau * sum((data == "1") & (Z == z))))
    }
    return(s)
  }

  # Resampling keeping proportions by Y and Z
  ids <- Resample(data = simul$Y, family = "binomial", resampling = BalancedResampling, Z = conf)

  # Bi-selection
  stab <- BiSelection(
    xdata = simul$X, ydata = simul$Y, family = "binomial",
    resampling = BalancedResampling, Z = conf
  )
})
