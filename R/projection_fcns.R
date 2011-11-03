

e0.proj.le.SDPropToLoess<-function(x,l.start,kap,n.proj=11, p1=9, p2=9){
  proj<-NULL
  proj[1]<-l.start
  for(a in 2:(n.proj+1)){
   proj[a]<-proj[a-1]+g.dl6(x,proj[a-1], p1=p1, p2=p2)+rnorm(1,mean=0,sd=(kap*loess.lookup(proj[a-1])))
  }
  return(proj)
}

e0.predict <- function(mcmc.set=NULL, end.year=2100, sim.dir=file.path(getwd(), 'bayesLife.output'),
                       replace.output=FALSE, nr.traj = NULL, thin=NULL, burnin=20000, save.as.ascii=1000,
                       output.dir = NULL, low.memory=TRUE, seed=NULL, verbose=TRUE){
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesLife.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesLife.mcmc.set.')
		}
	} else {                
		mcmc.set <- get.e0.mcmc(sim.dir, low.memory=low.memory, verbose=verbose)
	}
	if(!is.null(seed)) set.seed(seed)
	invisible(make.e0.prediction(mcmc.set, end.year=end.year,  
					replace.output=replace.output,  
					nr.traj=nr.traj, thin=thin, burnin=burnin, save.as.ascii=save.as.ascii,
					output.dir=output.dir, verbose=verbose))
}

e0.predict.extra <- function(sim.dir=file.path(getwd(), 'bayesLife.output'), 
					prediction.dir=sim.dir, 
					countries = NULL, save.as.ascii=1000, verbose=TRUE) {
	# Run prediction for given countries/regions (as codes). If they are not given it will be set to countries 
	# for which there are MCMC results but no prediction.
	# It is to be used after running run.e0.mcmc.extra
	
	mcmc.set <- get.e0.mcmc(sim.dir)
	if(is.null(mcmc.set))
		stop('Error in "sim.dir" argument.')
	pred <- get.e0.prediction(sim.dir=prediction.dir)
	if(is.null(pred))
		stop('Error in "prediction.dir" argument.')
	if(length(setdiff(pred$mcmc.set$meta$regions$country_code, mcmc.set$meta$regions$country_code)) > 0)
		stop('Prediction is inconsistent with the mcmc results. Use e0.predict.')
	if(is.null(countries)) {
		countries.idx <- (1:mcmc.set$meta$nr.countries)[!is.element(mcmc.set$meta$regions$country_code, 
												pred$mcmc.set$meta$regions$country_code)]
	} else {
		countries.idx <- (1:mcmc.set$meta$nr.countries)[is.element(mcmc.set$meta$regions$country_code,
												countries)]
	}
	if(length(countries.idx) == 0) {
		cat('\nNothing to be done.\n')
		return(invisible(pred))	
	}
	new.pred <- make.e0.prediction(mcmc.set, end.year=pred$end.year, replace.output=FALSE,
									nr.traj=pred$nr.traj, burnin=pred$burnin,
									countries=countries.idx, save.as.ascii=0, output.dir=prediction.dir,
									write.summary.files=FALSE, verbose=verbose)
									
	# merge the two predictions
	code.other.countries <- setdiff(pred$mcmc.set$meta$regions$country_code, 
									mcmc.set$meta$regions$country_code[countries.idx])
	idx.pred.others <- (1:pred$mcmc.set$meta$nr.countries)[is.element(pred$mcmc.set$meta$regions$country_code, 
												code.other.countries)]
	idx.other.countries <- (1:mcmc.set$meta$nr.countries)[is.element(mcmc.set$meta$regions$country_code,
												code.other.countries)]
												
	prev.pred <- pred
	pred$quantiles <- new.pred$quantiles
	pred$quantiles[idx.other.countries,,] <- prev.pred$quantiles[idx.pred.others,,]
	
	pred$traj.mean.sd <- new.pred$traj.mean.sd
	pred$traj.mean.sd[idx.other.countries,,] <- prev.pred$traj.mean.sd[idx.pred.others,,]
	
	pred$e0.matrix.reconstructed <- new.pred$e0.matrix.reconstructed
	pred$e0.matrix.reconstructed[,idx.other.countries] <- prev.pred$e0.matrix.reconstructed[,idx.pred.others]
	
	pred$mcmc.set <- mcmc.set
	
	# save updated prediction, convert trajectories and create summary files
	bayesLife.prediction <- pred
	prediction.file <- file.path(pred$output.dir, 'prediction.rda')
	save(bayesLife.prediction, file=prediction.file)
	
	bayesTFR:::do.convert.trajectories(pred=bayesLife.prediction, n=save.as.ascii, output.dir=pred$output.dir, 
							verbose=verbose)
	#do.write.projection.summary(pred=bayesLife.prediction, output.dir=pred$output.dir)
	
	cat('\nPrediction stored into', pred$output.dir, '\n')
	invisible(bayesLife.prediction)
}


make.e0.prediction <- function(mcmc.set, end.year=2100, replace.output=FALSE,
								nr.traj = NULL, thin=NULL, burnin=0, countries = NULL,
							    save.as.ascii=1000, output.dir = NULL, write.summary.files=TRUE, 
							    verbose=verbose){
	# if 'countries' is given, it is an index
	nr_project <- ceiling((end.year - mcmc.set$meta$present.year)/5)
	cat('\nPrediction from', mcmc.set$meta$present.year, 
			'(excl.) until', end.year, '(i.e.', nr_project, 'projections)\n\n')
			
	total.iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
	stored.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burnin)
	mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
	if(!is.null(nr.traj) && !is.null(thin)) {
		warning('Both nr.traj and thin are given. Argument thin will be ignored.')
		thin <- NULL
	}
	if(is.null(nr.traj)) nr.traj <- min(stored.iter, 2000)
	if(is.null(thin)) thin <- floor(stored.iter/nr.traj * mcthin)
	if(stored.iter <= 0 || thin == 0)
		stop('The number of simulations is 0. Burnin might be larger than the number of simulated values, or # trajectories is too big.')

	#setup output directory
	if (is.null(output.dir)) output.dir <- mcmc.set$meta$output.dir
	outdir <- file.path(output.dir, 'predictions')

	if(is.null(countries)) {
		if(!replace.output && has.e0.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
				' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
	} 
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	
	thinned.mcmc <- get.thinned.e0.mcmc(mcmc.set, thin=thin, burnin=burnin)
	has.thinned.mcmc <- !is.null(thinned.mcmc) && thinned.mcmc$meta$parent.iter == total.iter
	load.mcmc.set <- if(has.thinned.mcmc) thinned.mcmc
					 else create.thinned.e0.mcmc(mcmc.set, thin=thin, burnin=burnin, 
					 					output.dir=output.dir, verbose=verbose)
	nr_simu <- load.mcmc.set$mcmc.list[[1]]$finished.iter	
	if (verbose) cat('Load world-level parameters.\n')
	var.par.names <- c('omega')
	# load only the first par to check if everything is o.k.
	var.par.values <- get.e0.parameter.traces(load.mcmc.set$mcmc.list, var.par.names, burnin=0)

	prediction.countries <- if(is.null(countries)) 1:mcmc.set$meta$nr.countries else countries
	nr_countries <- mcmc.set$meta$nr.countries
	e0.matrix.reconstructed <- get.e0.reconstructed(mcmc.set$meta$e0.matrix, mcmc.set$meta)
	le0.matrix <- dim(e0.matrix.reconstructed)[1]
	quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- array(NA, c(nr_countries, length(quantiles.to.keep), nr_project+1))
	dimnames(PIs_cqp)[[2]] <- quantiles.to.keep
	proj.middleyears <- seq(max(as.numeric(dimnames(e0.matrix.reconstructed)[[1]])), by=5, length=nr_project+1)
	dimnames(PIs_cqp)[[3]] <- proj.middleyears
	mean_sd <- array(NA, c(nr_countries, 2, nr_project+1))

	var.par.names.cs <- c('Triangle.c', 'k.c', 'z.c')
	
	for (country in prediction.countries){
	#########################################
		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
		if (verbose) {			
 			cat('e0 projection for country', country, country.obj$name, 
 						'(code', country.obj$code, ')\n')
 		}
 		if (is.element(country.obj$code, load.mcmc.set$meta$regions$country_code)) {
			cs.par.values <- get.e0.parameter.traces.cs(load.mcmc.set$mcmc.list, country.obj, 
								var.par.names.cs, burnin=0)
		} else { # there are no thinned traces for this country, use the full traces
			cs.par.values <- get.e0.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, 
								var.par.names.cs, burnin=burnin)
			selected.simu <- bayesTFR:::get.thinning.index(nr_simu, dim(cs.par.values)[1])
			if (length(selected.simu$index) < nr_simu)
				selected.simu$index <- sample(selected.simu$index, nr_simu, replace=TRUE)
			cs.par.values <- cs.par.values[selected.simu$index,]
		}
		missing <- is.na(e0.matrix.reconstructed[,country])
		nmissing <- sum(missing)
		if (verbose && nmissing > 0) 
			cat('\t', nmissing, 'data points reconstructed.\n')

		this.nr_project <- nr_project + nmissing
		this.T_end <- mcmc.set$meta$T.end.c[country]
		trajectories <- matrix(NA, this.nr_project+1, nr_simu)
    	for(j in 1:nr_simu){
           trajectories[,j]<-e0.proj.le.SDPropToLoess(cs.par.values[j,], 
           							mcmc.set$meta$e0.matrix[mcmc.set$meta$T.end.c[country], country], 
           							kap=var.par.values[j,'omega'],n.proj=this.nr_project,
           							p1=mcmc.set$meta$dl.p1, p2=mcmc.set$meta$dl.p2)
    	}
    	if (nmissing > 0) {
    		e0.matrix.reconstructed[(this.T_end+1):le0.matrix,country] <- apply(matrix(trajectories[2:(nmissing+1),],
    											 nrow=nmissing), 1, quantile, 0.5, na.rm = TRUE)
    		trajectories <- trajectories[(nmissing+1):nrow(trajectories),]
    		trajectories[1,] <- quantile(trajectories[1,], 0.5, na.rm = TRUE)
    	}
		save(trajectories, file = file.path(outdir, paste('traj_country', country.obj$code, '.rda', sep='')))
 		PIs_cqp[country,,] = apply(trajectories, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		mean_sd[country,1,] <- apply(trajectories, 1, mean, na.rm = TRUE)
 		mean_sd[country,2,] = apply(trajectories, 1, sd, na.rm = TRUE)
	}
	mcmc.set <- remove.e0.traces(mcmc.set)
	bayesLife.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr_simu,
				e0.matrix.reconstructed = e0.matrix.reconstructed,
				output.directory=outdir,
				mcmc.set=load.mcmc.set,
				nr.projections=nr_project,
				burnin=burnin,
				end.year=end.year),
				class='bayesLife.prediction')
				
	prediction.file <- file.path(outdir, 'prediction.rda')
	save(bayesLife.prediction, file=prediction.file)
	
	bayesTFR:::do.convert.trajectories(pred=bayesLife.prediction, n=save.as.ascii, output.dir=outdir, 
										verbose=verbose)
	#write.summary.files <- FALSE # TODO: remove this after the function is fixed
	if(write.summary.files)
		bayesTFR:::do.write.projection.summary(pred=bayesLife.prediction, output.dir=outdir)
	
	cat('\nPrediction stored into', outdir, '\n')
	invisible(bayesLife.prediction)
}

remove.e0.traces <- function(mcmc.set) {
	for (i in 1:length(mcmc.set$mcmc.list)) 
		mcmc.set$mcmc.list[[i]]$traces <- 0
	invisible(mcmc.set)
}

get.projection.summary.header.bayesLife.prediction <- function(pred, ...) 
		return (list(revision='RevID', variant='VarID', country='LocID', year='TimeID', tfr='e0'))
		
get.UN.variant.names.bayesLife.prediction <- function(pred, ...) 
		return(c('BHM median', 'BHM80 lower',  'BHM80 upper', 'BHM95 lower',  'BHM95 upper', 'Constant mortality'))
	
get.friendly.variant.names.bayesLife.prediction <- function(pred, ...)
	return(c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95', 'constant'))	
convert.e0.trajectories <- function(dir=file.path(getwd(), 'bayesLife.output'), 
								 n=1000, output.dir=NULL, 
								 verbose=FALSE) {
	# Converts all trajectory rda files into UN ascii, selecting n trajectories by equal spacing.
	if(n <= 0) return()
	pred <- get.e0.prediction(sim.dir=dir)
	if (is.null(output.dir)) output.dir <- pred$output.directory
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	cat('Converting trajectories from', dir, '\n')
	bayesTFR:::do.convert.trajectories(pred=pred, n=n, output.dir=output.dir, verbose=verbose)
}

write.e0.projection.summary <- function(dir=file.path(getwd(), 'bayesLife.output'), 
									 output.dir=NULL, revision=14) {
# Writes two prediction summary files, one in a user-friendly format, one in a UN-format.
	pred <- get.e0.prediction(sim.dir=dir)
	if (is.null(output.dir)) output.dir <- pred$output.directory
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	bayesTFR:::do.write.projection.summary(pred, output.dir, revision=revision)
}
		
					
get.traj.ascii.header.bayesLife.mcmc.meta <- function(meta, ...) 
	return (list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='e0'))
	
get.data.imputed.bayesLife.prediction <- function(pred, ...)
	return(get.e0.reconstructed(pred$e0.matrix.reconstructed, pred$mcmc.set$meta))

get.e0.reconstructed <- function(data, meta) {
	return(if(is.null(data)) meta$e0.matrix.all else data)
}

get.e0.shift <- function(country.code, pred) return(bayesTFR:::get.tfr.shift(country.code, pred))

e0.median.shift <- function(sim.dir, country, reset=FALSE, shift=0, from=NULL, to=NULL) {
	invisible(bayesTFR:::.bdem.median.shift(type='e0', sim.dir=sim.dir, country=country, reset=reset, 
				shift=shift, from=from, to=to))
}

e0.median.set <- function(sim.dir, country, values, years=NULL) {
	invisible(bayesTFR:::.bdem.median.set(type='e0', sim.dir=sim.dir, country=country, 
				values=values, years=years))
}