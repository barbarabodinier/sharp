test_that("PLS 'a la carte' is consistent with mixOmics", {
  for (scale in c(TRUE, FALSE)) {
    # Data simulation
    set.seed(1)
    simul <- SimulateRegression(n = 20, pk = c(5, 5, 5), family = "gaussian")
    x <- simul$xdata + 10
    y <- simul$ydata + 10

    # Parameters
    ncomp <- 3

    # Using mixOmics
    mypls1 <- mixOmics::pls(X = x, Y = y, ncomp = ncomp, scale = scale)

    # Using our implementation
    mypls2 <- PLS(xdata = x, ydata = y, ncomp = ncomp, scale = scale)

    # Checking consistency in X-weights (i.e. that Wmat are X-loadings in mixOmics)
    expect_true(min(diag(cor(mypls1$loadings$X, mypls2$Wmat))) > 0.99)

    # Checking consistency in X-loadings (i.e. that Pmat are mat.c in mixOmics)
    expect_true(min(diag(cor(mypls1$mat.c, mypls2$Pmat))) > 0.99)

    # Checking consistency in X-scores (i.e. that Tmat are X-variates in mixOmics)
    expect_true(min(diag(cor(mypls1$variates$X, mypls2$Tmat))) > 0.99)

    # Checking consistency in Y-scores (i.e. that Umat are Y-variates in mixOmics)
    expect_true(min(diag(cor(mypls1$variates$Y, mypls2$Umat))) > 0.99)

    # New data simulation
    set.seed(1)
    simul <- SimulateRegression(n = 20, pk = c(5, 5, 5), family = "gaussian")
    newdata <- simul$xdata + 1

    # Predicting values
    predicted1 <- predict(mypls1, newdata = newdata)$predict
    predicted2 <- PredictPLS(xdata = newdata, pls = mypls2)

    # Checking consistency in predictions
    for (comp in 1:ncomp) {
      expect_true(min(diag(cor(predicted1[, , comp], predicted2[, , comp]))) > 0.99)
    }
    # plot(predicted1[,1,1], y[,1], panel.first=abline(0,1,col="red"))
    # plot(predicted2[,1,1], y[,1], panel.first=abline(0,1,col="red"))
  }
})
