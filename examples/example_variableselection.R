# Variable selection in linear regression
simul=SimulateRegression(n=100, pk=20, family="gaussian") # simulated data
out=VariableSelection(xdata=simul$X, ydata=simul$Y, family="gaussian")
