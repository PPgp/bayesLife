get.e0.mcmc <- function(sim.dir=file.path(getwd(), 'bayesLife.output'), chain.ids=NULL, 
						low.memory=TRUE, burnin=0, verbose=FALSE) {
	############
	# Returns an object of class bayesLife.mcmc.set
	############
	mcmc.file.path <- file.path(sim.dir, 'bayesLife.mcmc.meta.rda')
	if(!file.exists(mcmc.file.path)) {
		warning('File ', mcmc.file.path, ' does not exist.')
		return(NULL)
	}
	load(file=mcmc.file.path)
	bayesLife.mcmc.meta$output.dir <- sim.dir
	if (is.null(chain.ids)) {
		mc.dirs.short <- list.files(sim.dir, pattern='^mc[0-9]+', full.names=FALSE)
		chain.ids <- as.integer(substring(mc.dirs.short, 3))
	} else {
		mc.dirs.short <- paste('mc', chain.ids, sep='')
	}
	ord.idx <- order(chain.ids)
	mc.dirs.short <- mc.dirs.short[ord.idx]
	chain.ids <- chain.ids[ord.idx]
	mcmc.chains <- list()
	counter<-1
	for (imc.d in chain.ids) {
		if (verbose)
			cat('Loading chain', imc.d, 'from disk. ')
		load(file=file.path(sim.dir, mc.dirs.short[counter], 'bayesLife.mcmc.rda'))
		bayesLife.mcmc$meta <- bayesLife.mcmc.meta
		if (!low.memory) { # load full mcmc traces
			th.burnin <- bayesTFR:::get.thinned.burnin(bayesLife.mcmc, burnin)
			bayesLife.mcmc$traces <- load.e0.parameter.traces.all(bayesLife.mcmc, burnin=th.burnin)
			bayesLife.mcmc$traces.burnin <- th.burnin
		} else { # traces will be loaded as they are needed
			bayesLife.mcmc$traces <- 0
			bayesLife.mcmc$traces.burnin <- 0
		}
		bayesLife.mcmc$output.dir <- mc.dirs.short[counter]
		if (verbose)
			cat('(mcmc.list[[', counter, ']]).\n')
		mcmc.chains[[counter]] <- bayesLife.mcmc
		counter <- counter+1
	}
	names(mcmc.chains) <- chain.ids
	return(structure(list(meta=bayesLife.mcmc.meta, 
                          mcmc.list=mcmc.chains), class='bayesLife.mcmc.set'))
}

has.e0.mcmc <- function(sim.dir) {
	return(file.exists(file.path(sim.dir, 'bayesLife.mcmc.meta.rda')))
}

e0.mcmc <- function(mcmc.set, chain.id=1) return (mcmc.set$mcmc.list[[chain.id]])

e0.mcmc.list <- function(mcmc.set, chain.ids=NULL) 
	return(bayesTFR:::tfr.mcmc.list(mcmc.set=mcmc.set, chain.ids=chain.ids))

has.e0.prediction <- function(mcmc=NULL, sim.dir=NULL) {
	if (!is.null(mcmc)) sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
	if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
	if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
	return(FALSE)
}

get.e0.prediction <- function(mcmc=NULL, sim.dir=NULL, mcmc.dir=NULL) {
	############
	# Returns an object of class bayesLife.prediction
	# Set mcmc.dir to NA, if the prediction object should not have a pointer 
	# to the corresponding mcmc traces

	############
	if (!is.null(mcmc)) 
		sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
	if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
	output.dir <- file.path(sim.dir, 'predictions')
	pred.file <- file.path(output.dir, 'prediction.rda')
	if(!file.exists(pred.file)) {
		warning('File ', pred.file, ' does not exist.')
		return(NULL)
	}
	load(file=pred.file)
	bayesLife.prediction$output.directory <- output.dir
	
	pred <- bayesLife.prediction
	# re-route mcmcs if necessary
	if(!is.null(mcmc.dir) || !has.e0.mcmc(pred$mcmc.set$meta$output.dir)) {
		if((!is.null(mcmc.dir) && !is.na(mcmc.dir)) || is.null(mcmc.dir)) {
			new.path <- file.path(sim.dir, basename(pred$mcmc.set$meta$output.dir))
			if (has.e0.mcmc(new.path)) pred$mcmc.set <- get.e0.mcmc(new.path)
			else {
				est.dir <- if(is.null(mcmc.dir)) sim.dir else mcmc.dir
				pred$mcmc.set <- get.e0.mcmc(est.dir)
			}
		}
	}
	return(pred)
}

get.e0.convergence.all <- function(sim.dir=file.path(getwd(), 'bayesLife.output')) {
	return(bayesTFR:::.do.get.convergence.all('e0', 'bayesLife', sim.dir=sim.dir))
}

get.e0.convergence <- function(sim.dir=file.path(getwd(), 'bayesLife.output'), 
									thin=120, burnin=20000) {
	file.name <- file.path(sim.dir, 'diagnostics', paste('bayesLife.convergence_', 
							thin, '_', burnin, '.rda', sep=''))
	if(!file.exists(file.name)){
		warning('Convergence diagnostics in ', sim.dir, ' for burnin=', burnin, 
					' and thin=', thin, ' does not exist.')
		return(NULL)
	}
	load(file.name)
	return(bayesLife.convergence)
}

get.e0.parameter.traces <- function(mcmc.list, par.names=e0.parameter.names(), burnin=0,
									thinning.index=NULL, thin=NULL) {
	# get parameter traces either from disk or from memory, if they were already loaded
	mcmc.list <- get.mcmc.list(mcmc.list)
	return(bayesTFR:::do.get.tfr.parameter.traces(is.cs=FALSE, mcmc.list=mcmc.list, par.names=par.names, 
										burnin=burnin, thinning.index=thinning.index, thin=thin))
}

get.e0.parameter.traces.cs <- function(mcmc.list, country.obj, par.names=e0.parameter.names.cs(), 
									   burnin=0, thinning.index=NULL, thin=NULL) {
	# country.obj is result of get.country.object()
	# get traces for country-specific parameters either from disk or from memory, if they were already loaded
	mcmc.list <- get.mcmc.list(mcmc.list)
	return(bayesTFR:::do.get.tfr.parameter.traces(is.cs=TRUE, mcmc.list=mcmc.list, par.names=par.names,
										 country.obj=country.obj,
										burnin=burnin, thinning.index=thinning.index, thin=thin))
}

bdem.parameter.traces.bayesLife.mcmc <- function(mcmc, par.names, ...) {
	# Load traces from the disk
	all.standard.names <- c(e0.parameter.names(), e0.parameter.names.cs())
	return(bayesTFR:::.do.get.traces(mcmc, par.names=par.names, ..., 
							all.standard.names=all.standard.names))
}

load.e0.parameter.traces <- function(mcmc, par.names=e0.parameter.names(), burnin=0, thinning.index=NULL) {
 	return(bdem.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index))
}

load.e0.parameter.traces.cs <- function(mcmc, country, par.names=e0.parameter.names.cs(), burnin=0, 
										thinning.index=NULL) {
 	return(bdem.parameter.traces(mcmc, par.names, paste('_country', country, sep=''),
						par.names.postfix=paste('_c', country, sep=''), burnin=burnin, 
						thinning.index=thinning.index))
}

load.e0.parameter.traces.all <- function(mcmc, par.names=e0.parameter.names(), 
										 par.names.cs=e0.parameter.names.cs(),
										 burnin=0, thinning.index=NULL) {
	result <- load.e0.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index)
	for (country in 1:mcmc$meta$nr.countries) {
		result <- cbind(result, 
						load.e0.parameter.traces.cs(mcmc, 
												    get.country.object(country, 
												         mcmc$meta, index=TRUE)$code, 
												    par.names.cs, burnin=burnin,
												    thinning.index=thinning.index))
	}
	return (result)
}

e0.get.all.parameter.names <- function() {
    # Parameters with the number of its values
	return(list(Triangle=4, k=1, z=1, lambda=4, lambda.k=1, lambda.z=1, omega=1)) 
}				

e0.get.all.parameter.names.cs <- function() {
    # Country-specific parameters with the number of its values
	return(list(Triangle.c=4, k.c=1, z.c=1))
}

e0.get.all.parameter.names.extended <- function(cs=FALSE) {
	pars <- c()
	all.pars <- if(cs) e0.get.all.parameter.names.cs() else e0.get.all.parameter.names()
	for (ipar in 1:length(all.pars)) {
		par <- all.pars[ipar]
		paropt <- unlist(par)
		name <- names(par)
		pars <- c(pars, if (paropt > 1) paste(name, 1:paropt, sep='_') else name)
	}
	return(pars)
}

e0.parameter.names.extended <- function() 
	return(e0.get.all.parameter.names.extended())
	
e0.parameter.names.cs.extended <- function(country.code=NULL) {
	# return a list of cs-parameters extended by i if necessary
	pars <- e0.get.all.parameter.names.extended(cs=TRUE)
	if(!is.null(country.code)) pars <- paste(pars, '_c', country.code, sep='')
	return(pars)
}

e0.parameter.names <- function() return(names(e0.get.all.parameter.names()))
e0.parameter.names.cs <- function() return(names(e0.get.all.parameter.names.cs()))

get.mcmc.list.bayesLife.mcmc.set <- function(mcmc.list, ...) return(mcmc.list$mcmc.list)
get.mcmc.list.bayesLife.mcmc <- function(mcmc.list, ...) return(list(mcmc.list))
get.mcmc.list.bayesLife.prediction <- function(mcmc.list, ...) return(mcmc.list$mcmc.set$mcmc.list)

get.mcmc.meta.bayesLife.mcmc.set <- function(meta, ...) return(meta$meta)
get.mcmc.meta.bayesLife.mcmc.meta <- function(meta, ...) return(meta)
get.mcmc.meta.bayesLife.mcmc <- function(meta, ...) return(meta$meta)
get.mcmc.meta.bayesLife.prediction <- function(meta, ...) return(meta$mcmc.set$meta)


summary.bayesLife.mcmc <- function(object, country=NULL, 
								par.names=e0.parameter.names(), 
								par.names.cs=e0.parameter.names.cs(), ...)
	return (bayesTFR:::summary.bayesTFR.mcmc(object, country=country,
				par.names=par.names, par.names.cs=par.names.cs, ...))

summary.bayesLife.mcmc.set <- function(object, country=NULL, chain.id=NULL, 
								par.names=e0.parameter.names(), 
								par.names.cs=e0.parameter.names.cs(), 
								meta.only=FALSE, thin=1, burnin=0, ...) {
	if(is.null(country) & missing(par.names.cs)) par.names.cs <- NULL
	cat('\nMCMC parameters estimated for', object$meta$nr.countries, 'countries.')
	cat('\nHyperparameters estimated using', object$meta$nr.countries.estimation, 'countries.\n')
	cat('\nWPP:', object$meta$wpp.year)
	cat('\nInput data: e0 for period', object$meta$start.year, '-', object$meta$present.year,'.')
	cat('\n')
	if(meta.only) {
		get.iter <- function(x) x$finished.iter
		res <- list(nr.chains=object$meta$nr.chains, 
					iters=sapply(object$mcmc.list, get.iter),
					thin=object$mcmc.list[[1]]$thin)
		class(res) <- 'summary.bayesLife.mcmc.set.meta'

		return(res)
	} 
	if (!is.null(chain.id))
		return(summary(object$mcmc.list[[chain.id]], country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin, ...))
	if (!is.null(country)) {
		country.obj <- get.country.object(country, object$meta)
		cat('\nCountry:', country.obj$name, '\n')
		country <- country.obj$code
	}
	summary(coda.mcmc.list(object, country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin), ...)
}

print.summary.bayesLife.mcmc.set.meta <- function(x, ...) return(bayesTFR:::print.summary.bayesTFR.mcmc.set.meta(x, ...))

summary.bayesLife.prediction <- function(object, country=NULL, compact=TRUE, ...) {
	res <- bayesTFR:::get.prediction.summary.data(object, 
				unchanged.pars=c('burnin', 'nr.traj'), 
				country=country, compact=compact)
	class(res) <- 'summary.bayesLife.prediction'
	return(bayesTFR:::.update.summary.data.by.shift(res, object, country))
}

print.summary.bayesLife.prediction <- function(x, digits = 3, ...) {
	cat('\nProjections:', length(x$projection.years), '(', x$projection.years[1], '-', 
					x$projection.years[length(x$projection.years)], ')')
	cat('\nTrajectories:', x$nr.traj)
	cat('\nBurnin:', x$burnin)

	if(!is.null(x$country.name)) {
		cat('\nCountry:', x$country.name, '\n')
		cat('\nProjected Life Expectancy:')
		#if(x$manual) cat(' (values have been manually modified)')
		cat('\n')
		print(x$projections, digits=digits, ...)
	}
}

summary.bayesLife.convergence <- function(object, expand=FALSE, ...) {
	return(bayesTFR:::summary.bayesTFR.convergence(object, expand=expand, ...))
}
	
get.thinned.e0.mcmc <- function(mcmc.set, thin=1, burnin=0) {
	dir.name <- file.path(mcmc.set$meta$output.dir, paste('thinned_mcmc', thin, burnin, sep='_'))
	if(file.exists(dir.name)) return(get.e0.mcmc(dir.name))
	return(NULL)
}
	
create.thinned.e0.mcmc <- function(mcmc.set, thin=1, burnin=0, output.dir=NULL, verbose=TRUE) {
	#Return a thinned mcmc.set object with burnin removed and all chanins collapsed into one
	mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
	thin <- max(c(thin, mcthin))
	meta <- mcmc.set$meta
	meta$output.dir <- file.path(
		if(is.null(output.dir)) meta$output.dir else output.dir, 
		paste('thinned_mcmc', thin, burnin, sep='_'))
	if(!file.exists(meta$output.dir)) 
		dir.create(meta$output.dir, recursive=TRUE)
	total.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burnin=burnin)
	meta$is.thinned <- TRUE
	meta$parent.iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
	meta$parent.meta <- mcmc.set$meta
	meta$nr.chains <- 1
	
	if(verbose) cat('\nStoring thinned mcmc:')
	# store the meta object
	store.bayesLife.meta.object(meta, meta$output.dir)

	
	thin.index <- if(thin > mcthin) unique(round(seq(thin, total.iter, by=thin/mcthin))) else 1:total.iter
	nr.points <- length(thin.index)
	
	#create one collapsed mcmc
	thinned.mcmc <- mcmc.set$mcmc.list[[1]]
	thinned.mcmc$meta <- meta
	thinned.mcmc$thin <- 1
	thinned.mcmc$id <- 1
	thinned.mcmc$traces <- 0
	thinned.mcmc$length <- nr.points
	thinned.mcmc$finished.iter <- nr.points
	
	outdir.thin.mcmc <- file.path(meta$output.dir, 'mc1')
	if(!file.exists(outdir.thin.mcmc)) dir.create(outdir.thin.mcmc)
	thinned.mcmc$output.dir <- 'mc1'
	store.bayesLife.object(thinned.mcmc, outdir.thin.mcmc)
	
	if(verbose) cat('\nStoring country-independent parameters ...')
	for (par in e0.parameter.names()) {
		values <- get.e0.parameter.traces(mcmc.set$mcmc.list, par, burnin,
											thinning.index=thin.index)
		bayesTFR:::write.values.into.file.cindep(par, values, outdir.thin.mcmc)
	}
	if(verbose) cat('done.\nStoring country-specific parameters ...')
	par.names.cs <- e0.parameter.names.cs()
	for (country in 1:mcmc.set$meta$nr.countries){
		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
		for (par in par.names.cs) {
			values <- get.e0.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
											burnin=burnin, thinning.index=thin.index)
			bayesTFR:::write.values.into.file.cdep(par, values, outdir.thin.mcmc, country.code=country.obj$code)
		}
	}
	#if (mcmc.set$meta$nr_countries > mcmc.set$meta$nr_countries_estimation) {
	#	.update.thinned.extras(mcmc.set, (mcmc.set$meta$nr_countries_estimation+1):mcmc.set$meta$nr_countries,
	#							burnin=burnin, nr.points=nr.points, dir=outdir.thin.mcmc, verbose=verbose)
	#}
	invisible(structure(list(meta=meta, mcmc.list=list(thinned.mcmc)), class='bayesLife.mcmc.set'))
}

get.e0.trajectories <- function(e0.pred, country) 
	return(bayesTFR:::get.tfr.trajectories(e0.pred, country=country))
	
	
get.nr.countries.bayesLife.mcmc.meta <- function(meta, ...) return (meta$nr.countries)
get.nr.countries.est.bayesLife.mcmc.meta <- function(meta, ...) return (meta$nr.countries.estimation)

get.data.matrix.bayesLife.mcmc.meta <- function(meta, ...) return (meta$e0.matrix)

coda.mcmc.bayesLife.mcmc <- function(mcmc, country=NULL, par.names=e0.parameter.names(), 
						par.names.cs=e0.parameter.names.cs(), ...)
	return(bayesTFR:::coda.mcmc.bayesTFR.mcmc(mcmc, country=country, par.names=par.names, 
												par.names.cs=par.names.cs, ...))
												
e0.coda.mcmc.list <- function(mcmc.list=NULL, country=NULL, chain.ids=NULL,
							sim.dir=file.path(getwd(), 'bayesLife.output'), 
							par.names=e0.parameter.names(), 
							par.names.cs=e0.parameter.names.cs(), 
							low.memory=FALSE, ...) {
	# return a list of mcmc objects that can be analyzed using the coda package
	if (is.null(mcmc.list)) {
		mcmc.list <- get.e0.mcmc(sim.dir, chain.ids=chain.ids, low.memory=low.memory)$mcmc.list
	} else {
		mcmc.list <- get.mcmc.list(mcmc.list)
		if (!is.null(chain.ids)) {
			mcmc.list <- mcmc.list[chain.ids]
		}
	}
	result <- list()
	i <- 1
	for(mcmc in mcmc.list) {
		result[[i]] <- coda.mcmc(mcmc, country=country, par.names=par.names, par.names.cs=par.names.cs, ...)
		i <- i+1
	}
	return(mcmc.list(result))
}

get.countries.index.bayesLife.mcmc.meta  <- function(meta, ...) return (1:meta$nr.countries)

get.countries.table.bayesLife.mcmc.set <- function(object, ...) 
	return(bayesTFR:::get.countries.table.bayesTFR.mcmc(object,...))
get.countries.table.bayesLife.prediction <- function(object, ...) 
	return(bayesTFR:::get.countries.table.bayesTFR.prediction(object,...))