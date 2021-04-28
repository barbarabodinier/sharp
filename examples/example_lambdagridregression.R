# Lambda grid for linear regression
simul=SimulateRegression(n=100, pk=20, family="gaussian") # simulated data
Lambda=LambdaGridRegression(xdata=simul$X, ydata=simul$Y,
                            family="gaussian", Lambda_cardinal=20)

# Grid can be used in VariableSelection()
out=VariableSelection(xdata=simul$X, ydata=simul$Y,
                      family="gaussian", Lambda=Lambda)
print(SelectedVariables(out))

# For use with gglasso (group LASSO)
require(gglasso)
ManualGridGroupLasso=function(x, y, family, ...){
  if (family=="gaussian"){
    return(cv.gglasso(x=x, y=y, pred.loss="L1", ...))
  }
}
Lambda=LambdaGridRegression(xdata=simul$X, ydata=simul$Y,
                            family="gaussian", Lambda_cardinal=20,
                            implementation="ManualGridGroupLasso",
                            group=rep(1:4, each=5))
