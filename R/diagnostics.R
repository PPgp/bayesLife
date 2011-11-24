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

GoF.dl <- function(sim.dir, pi=c(80,90,95)) {
	pred <- get.e0.prediction(sim.dir)
	meta <- pred$mcmc.set$meta
	nr.countries <- get.nr.countries.est(pred$mcmc.set$meta)
	T.total <- nrow(meta$e0.matrix)
	al.low <- (1 - pi/100)/2
	al.high <- 1 - al.low
	total.GoF <- rep(0, length(pi))
	time.GoF <- matrix(0, nrow=length(pi), ncol=T.total-1)
	country.GoF <- matrix(0, nrow=length(pi), ncol=nr.countries)
	counter.total <- (T.total - 1)*nr.countries
	for(icountry in 1: nr.countries) {
		country.code <- meta$regions$country_code[icountry]
    	observed <- diff(meta$e0.matrix[1:T.total, icountry])
		x <- meta$e0.matrix[1:(T.total - 1), icountry]
		dlc <- .get.dlcurves(x, pred$mcmc.set$mcmc.list, country.code, icountry, burnin=0, 
							nr.curves=2000, predictive.distr=TRUE)
		for (i in 1:length(pi)) {
        	dlpi <- apply(dlc, 2, quantile, c(al.low[i], al.high[i]))
        	country.GoF[i,icountry] <- sum(observed >= dlpi[1,] & observed <= dlpi[2,])
        	total.GoF[i] <- total.GoF[i] + country.GoF[i,icountry]
        	for(time in 1:(T.total-1)) {
        		time.GoF[i,time] <- time.GoF[i,time] + (observed[time] >= dlpi[1,time] & observed[time] <= dlpi[2,time])
        	}
        }
	}
	total.GoF <- total.GoF/counter.total
	pi.names <- paste(pi, '%', sep='')
	names(total.GoF) <- pi.names
	time.GoF <- time.GoF/nr.countries
	dimnames(time.GoF) <- list(pi.names, rownames(meta$e0.matrix)[2:T.total])
	country.GoF <- country.GoF/(T.total-1)
	dimnames(country.GoF) <- list(pi.names, meta$regions$country_code[1:nr.countries])
	return(list(total=total.GoF, time=time.GoF, country=country.GoF))
}