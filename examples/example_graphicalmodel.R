# Data simulation
set.seed(1)
simul=SimulateGraphical(n=100, pk=20, nu=0.1)
out=GraphicalModel(data=simul$data)
