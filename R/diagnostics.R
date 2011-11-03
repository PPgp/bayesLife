e0.raftery.diag <- function(mcmc=NULL, 
							 sim.dir=file.path(getwd(), 'bayesLife.output'),
							 burnin=0, country=NULL,
							 par.names = e0.parameter.names(),
							 par.names.cs = e0.parameter.names.cs(),
							 country.sampling.prop=1,
							 verbose=TRUE, ...) {
return(bayesTFR:::tfr.raftery.diag(mcmc=mcmc, sim.dir=sim.dir, burnin=burnin,
						country=country, par.names=par.names, par.names.cs=par.names.cs,
						country.sampling.prop=country.sampling.prop, verbose=verbose, ...))
}

e0.diagnose <- function(sim.dir, thin=120, burnin=20000, express=FALSE, 
						country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
	invisible(bayesTFR:::.do.diagnose(type='e0', class.name='bayesLife.convergence', 
							sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
							country.sampling.prop=country.sampling.prop, keep.thin.mcmc=keep.thin.mcmc,							verbose=verbose))
	
	
}
