library(bayesLife)
source('test_functions.R')

cran <- TRUE
wpp <- 2019
test.get.wpp.data(wpp)
test.existing.simulation()
test.e0trajectories()
test.get.parameter.traces()


## Time-expensive tests
if(!cran) {
	test.DLcurve()
	test.plot.density()
	test.plot.map()
	test.estimate.mcmc(compression='xz')
	test.estimate.mcmc.with.suppl.data()
	test.estimate.mcmc.with.suppl.data(compression='bz')
	test.plot.all()
	test.estimate.mcmc.with.overwrites()
	test.run.mcmc.simulation.auto()
	test.run.mcmc.simulation.auto(compression='gz')
	test.run.mcmc.simulation.auto.parallel()
	test.imputation()
	test.subnational.predictions()
	test.my.locations.extra()
	test.reproduce.simulation()
	for (wpp in c(2019, 2024)){
	    test.estimate.mcmc(wpp.year = wpp)
	    test.run.annual.simulation(wpp.year = wpp)
	}
	for(wpp.year in c(2010, 2012, 2015, 2017, 2022, 2024))
	    test.get.wpp.data(wpp.year)
}
