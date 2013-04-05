library(bayesLife)
source('test_functions.R')

test.get.wpp.data()
test.existing.simulation()
test.e0trajectories()
test.get.parameter.traces()


## Time-expensive tests
# test.DLcurve()
# test.plot.density()
# test.plot.map()
# test.estimate.mcmc()
# test.estimate.mcmc(compression='xz')
# test.estimate.mcmc.with.suppl.data()
# test.estimate.mcmc.with.suppl.data(compression='bz')
# test.plot.all()
# test.estimate.mcmc.with.overwrites()
# test.run.mcmc.simulation.auto()
# test.run.mcmc.simulation.auto(compression='gz')
