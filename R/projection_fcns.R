

e0.proj.le.SDPropToLoess<-function(x,l.start,kap,n.proj=11, p1=9, p2=9, const.var=FALSE){
  proj<-NULL
  proj[1]<-l.start
  for(a in 2:(n.proj+1)){
  	proj[a]<-proj[a-1]+g.dl6(x,proj[a-1], p1=p1, p2=p2)+rnorm(1,mean=0,
  				sd=(kap*if(const.var) 1 else loess.lookup(proj[a-1])))
  }
  return(proj)
}

e0.predict <- function(mcmc.set=NULL, end.year=2100, sim.dir=file.path(getwd(), 'bayesLife.output'),
                       replace.output=FALSE, nr.traj = NULL, thin=NULL, burnin=20000, 
                       use.diagnostics=FALSE, save.as.ascii=1000,
                       output.dir = NULL, low.memory=TRUE, seed=NULL, verbose=TRUE){
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesLife.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesLife.mcmc.set.')
		}
	} else {                
		mcmc.set <- get.e0.mcmc(sim.dir, low.memory=low.memory, verbose=verbose)
	}
	if(!is.null(seed)) set.seed(seed)
		# Get argument settings from existing convergence diagnostics
	if(use.diagnostics) {
		diag.list <- get.e0.convergence.all(mcmc.set$meta$output.dir)
		ldiag <- length(diag.list)
		if (ldiag == 0) stop('There is no diagnostics available. Use manual settings of "nr.traj" or "thin".')
		use.nr.traj <- use.burnin <- rep(NA, ldiag)
		for(idiag in 1:ldiag) {
			if (bayesTFR:::has.mcmc.converged(diag.list[[idiag]])) {
				use.nr.traj[idiag] <- diag.list[[idiag]]$use.nr.traj
				use.burnin[idiag] <- diag.list[[idiag]]$burnin
			}
		}
		if(all(is.na(use.nr.traj)))
			stop('There is no diagnostics indicating convergence of the MCMCs. Use manual settings of "nr.traj" or "thin".')
		# Try to select those that suggest nr.traj >= 2000 (take the minimum of those)
		traj.is.notna <- !is.na(use.nr.traj)
		larger2T <- traj.is.notna & use.nr.traj>=2000
		nr.traj.idx <- if(sum(larger2T)>0) (1:ldiag)[larger2T][which.min(use.nr.traj[larger2T])] 
						else (1:ldiag)[traj.is.notna][which.max(use.nr.traj[traj.is.notna])]
		nr.traj <- use.nr.traj[nr.traj.idx]
		burnin <- use.burnin[nr.traj.idx]
		if(verbose)
			cat('\nUsing convergence settings: nr.traj=', nr.traj, ', burnin=', burnin, '\n')
	}

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
									force.creating.thinned.mcmc=TRUE,
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
							    force.creating.thinned.mcmc=FALSE,
							    verbose=verbose){
	data(loess_sd)
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
	else {
		if (nr.traj > stored.iter) 
			warning('nr.traj is larger than the available MCMC sample. Only ', stored.iter, ' trajectories will be generated.')
		nr.traj <- min(nr.traj, stored.iter)	
	}
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
		write.to.disk <- TRUE
		if(!file.exists(outdir)) 
			dir.create(outdir, recursive=TRUE)
	} else write.to.disk <- FALSE
	
	thinned.mcmc <- get.thinned.e0.mcmc(mcmc.set, thin=thin, burnin=burnin)
	has.thinned.mcmc <- (!is.null(thinned.mcmc) && thinned.mcmc$meta$parent.iter == total.iter 
							&& mcmc.set$meta$nr.countries == thinned.mcmc$meta$nr.countries)

	load.mcmc.set <- if(has.thinned.mcmc && !force.creating.thinned.mcmc) thinned.mcmc
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
	#for (country in c(23)){
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
		this.T_end <- mcmc.set$meta$Tc.index[[country]][length(mcmc.set$meta$Tc.index[[country]])]
		nmissing <- le0.matrix - this.T_end
		missing <- (this.T_end+1):le0.matrix

		if (verbose && nmissing > 0) 
			cat('\t', nmissing, 'data points reconstructed.\n')

		this.nr_project <- nr_project + nmissing
		#sum.delta <- apply(cs.par.values[,1:4], 1, sum)
		#use.traj <- which(sum.delta <= 110)
		trajectories <- matrix(NA, this.nr_project+1, nr_simu)
		#trajectories <- matrix(NA, this.nr_project+1, length(use.traj))
    	for(j in 1:nr_simu){
    	#for(j in 1:length(use.traj)){
    		#k <- use.traj[j]
           trajectories[,j]<-e0.proj.le.SDPropToLoess(cs.par.values[j,], 
           							mcmc.set$meta$e0.matrix[this.T_end, country], 
           							kap=var.par.values[j,'omega'],n.proj=this.nr_project,
           							p1=mcmc.set$meta$dl.p1, p2=mcmc.set$meta$dl.p2, 
           							const.var=mcmc.set$meta$constant.variance)
    	}
    	if (nmissing > 0) {
    		e0.matrix.reconstructed[(this.T_end+1):le0.matrix,country] <- apply(matrix(trajectories[2:(nmissing+1),],
    											 nrow=nmissing), 1, quantile, 0.5, na.rm = TRUE)
    		trajectories <- trajectories[(nmissing+1):nrow(trajectories),]
    		trajectories[1,] <- quantile(trajectories[1,], 0.5, na.rm = TRUE)
    	}
    	#stop('')
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
		
	if(write.to.disk) {		
		prediction.file <- file.path(outdir, 'prediction.rda')
		save(bayesLife.prediction, file=prediction.file)
	
		bayesTFR:::do.convert.trajectories(pred=bayesLife.prediction, n=save.as.ascii, output.dir=outdir, 
										verbose=verbose)
		if(write.summary.files)
			bayesTFR:::do.write.projection.summary(pred=bayesLife.prediction, output.dir=outdir)
	
		cat('\nPrediction stored into', outdir, '\n')
	}
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

e0.jmale.estimate <- function(mcmc.set, countries.index=NULL, 
								estDof.eq1 = TRUE, start.eq1 = list(dof = 2), 
								min.e0.eq2 = 80, estDof.eq2 = TRUE, start.eq2 = list(dof = 2), 
								my.e0.file=NULL, verbose=FALSE) {
	# Estimate coefficients for joint prediction of female and male e0
	require(hett)
	if (is.null(countries.index)) countries.index <- 1:get.nr.countries.est(mcmc.set$meta)
	e0f.data <- get.data.matrix(mcmc.set$meta)[,countries.index]
	e0m.data <- get.wpp.e0.data.for.countries(mcmc.set$meta, sex='M', my.e0.file=my.e0.file,
					verbose=verbose)$e0.matrix[,countries.index]
	T <- dim(e0f.data)[1] - 1 
	if(verbose) {
		cat('\nEstimating coefficients for joint female and male prediction.')
		cat('\nUsing', length(countries.index), 'countries and', T+1, 'time periods.\n\n')
	}
	G <- e0f.data - e0m.data # observed gap

	e0.1953 <- rep(e0f.data['1953',], each=T)

	data.eq1 <- data.frame(
					G=as.numeric(G[2:nrow(G),]), # dependent variable
					# covariates
					e0.1953=e0.1953, 
					Gprev=as.numeric(G[1:T,]),
					e0=as.numeric(e0f.data[1:T,]),
					e0d75=pmax(0, e0f.data[1:T,]-75) 
				)
	fit.eq1 <- tlm(G~., data=data.eq1, estDof = estDof.eq1, start=start.eq1)
	if(verbose)
		print(summary(fit.eq1))
	errscale.eq1<-as.numeric(exp(fit.eq1$scale.fit$coefficients[1]))
	errsd.eq1<-sqrt(errscale.eq1)

	data.eq2 <- data.eq1[data.eq1$e0 >= min.e0.eq2,]
	if(verbose) 
		cat('\n\nUsing', nrow(data.eq2), 'data points for equation 2.\n\n')
	fit.eq2 <- tlm(G~-1+Gprev, data=data.eq2, start = start.eq2, estDof = estDof.eq2)
	if(verbose)
		print(summary(fit.eq2))
	errscale.eq2<-as.numeric(exp(fit.eq2$scale.fit$coefficients[1]))
	errsd.eq2<-sqrt(errscale.eq2)
	return(list(eq1 = list(coefficients=fit.eq1$loc.fit$coefficients, 
						   sigma=errsd.eq1, dof = fit.eq1$dof, fit=fit.eq1),
				eq2 = list(coefficients=fit.eq2$loc.fit$coefficients,
						   sigma=errsd.eq2, dof = fit.eq2$dof, fit=fit.eq2)
				))
}

e0.jmale.predict <- function(e0.pred, estimates=NULL, gap.lim=c(0,18), my.e0.file=NULL, verbose=TRUE, ...) {
	# Predicting male e0 from female predictions. estimates is the result of 
	# the e0.jmale.estimate function. If it is NULL, the estimation is performed 
	# using the ... arguments
	# If my.e0.file given, it should be a male e0 file. 
	
	meta <- e0.pred$mcmc.set$meta
	if (meta$sex != 'F') stop('The prediction object must be a result of FEMALE projections.')
	if(is.null(estimates)) 
		estimates <- e0.jmale.estimate(e0.pred$mcmc.set, verbose=verbose, ...)
	
	e0f.data <- get.e0.reconstructed(e0.pred$e0.matrix.reconstructed, meta)
	#Tc <- meta$T.end.c

	e0mwpp <- get.wpp.e0.data.for.countries(meta, sex='M', my.e0.file=my.e0.file, verbose=verbose)
	e0m.data <- e0mwpp$e0.matrix
	meta.changes <- list(sex='M', e0.matrix=e0m.data, e0.matrix.all=e0mwpp$e0.matrix.all, suppl.data=e0mwpp$suppl.data)
	maxe0 <- max(e0f.data)

	bayesLife.prediction <- e0.pred

	prediction.file <- file.path(e0.pred$output.directory, 'prediction.rda')
	joint.male <- e0.pred
	joint.male$output.directory <- file.path(e0.pred$output.directory, 'joint_male')
	joint.male$e0.matrix.reconstructed <- e0m.data
	joint.male$fit <- estimates
	joint.male$meta.changes <- meta.changes
	joint.male$mcmc.set <- NULL
	joint.male$joint.male <- NULL
	
	if(file.exists(joint.male$output.directory)) unlink(joint.male$output.directory, recursive=TRUE)
	dir.create(joint.male$output.directory, recursive=TRUE)
	bayesLife.prediction$joint.male <- joint.male

	quantiles <- array(NA, dim(e0.pred$quantiles))
	dimnames(quantiles) <- dimnames(e0.pred$quantiles)
	traj.mean.sd <- array(NA, dim(e0.pred$traj.mean.sd))
	dimnames(traj.mean.sd) <- dimnames(e0.pred$traj.mean.sd)

	quantiles.to.keep <- as.numeric(dimnames(e0.pred$quantiles)[[2]])
	for (icountry in 1:get.nr.countries(meta)) {
		country <- get.country.object(icountry, meta, index=TRUE)
		if(verbose)
			cat('\ne0 male projection for country', icountry, country$name, 
 						'(code', country$code, ')')
		trajectoriesF <- bayesTFR:::get.trajectories(e0.pred, country$code)$trajectories
		Mtraj <- matrix(NA, nrow=nrow(trajectoriesF), ncol=ncol(trajectoriesF))
		#G1 <- e0f.data[Tc[icountry],icountry] - e0m.data[Tc[icountry],icountry]
		Tc <- meta$Tc.index[[icountry]][length(meta$Tc.index[[icountry]])]
		G1 <- e0f.data[Tc,icountry] - e0m.data[Tc,icountry]
		for (itraj in 1:dim(trajectoriesF)[2]) {
			Mtraj[1,itraj] <- e0m.data[Tc,icountry]
			Gprev <- G1
			for(time in 2:dim(trajectoriesF)[1]) {
				if(trajectoriesF[time-1,itraj] <= maxe0) { # 1st part of Equation 3.1
					Gtdeterm <- (estimates$eq1$coefficients[1] + # intercept
					   			 estimates$eq1$coefficients['Gprev']*Gprev +
					   			 estimates$eq1$coefficients['e0.1953']*e0f.data['1953',icountry] +
					   			 estimates$eq1$coefficients['e0']*trajectoriesF[time-1,itraj] +
					   			 estimates$eq1$coefficients['e0d75']*max(0, trajectoriesF[time-1,itraj]-75)					   			)
					Gt <- Gtdeterm + estimates$eq1$sigma*rt(1,estimates$eq1$dof)
					while(Gt < gap.lim[1] || Gt > gap.lim[2]) 
						Gt <- Gtdeterm + estimates$eq1$sigma*rt(1,estimates$eq1$dof)
				} else {  # 2nd part of Equation 3.1
					Gtdeterm <- estimates$eq2$coefficients['Gprev']*Gprev
					Gt <- Gtdeterm + estimates$eq2$sigma*rt(1,estimates$eq2$dof)
					while(Gt < gap.lim[1] || Gt > gap.lim[2]) 
						Gt <- Gtdeterm + estimates$eq2$sigma*rt(1,estimates$eq2$dof)
				}
				Mtraj[time,itraj] <- trajectoriesF[time,itraj] - Gt
				Gprev <- Gt
			}
		}
		quantiles[icountry,,] = apply(Mtraj, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		traj.mean.sd[icountry,1,] <- apply(Mtraj, 1, mean, na.rm = TRUE)
 		traj.mean.sd[icountry,2,] = apply(Mtraj, 1, sd, na.rm = TRUE)
 		trajectories <- Mtraj
		save(trajectories, file = file.path(joint.male$output.directory, 
								paste('traj_country', country$code, '.rda', sep='')))
		bayesLife.prediction$joint.male$quantiles <- quantiles
		bayesLife.prediction$joint.male$traj.mean.sd <- traj.mean.sd
		save(bayesLife.prediction, file=prediction.file)
	}
	cat('\nPrediction stored into', joint.male$output.directory, '\n')
	invisible(bayesLife.prediction)
}

get.e0.jmale.prediction <- function(e0.pred) {
	male.pred <- e0.pred$joint.male
	if(is.null(male.pred)) stop('A joint male prediction does not exist for the given object. Use e0.jmale.predict to simulate male projections from existing female projections.')
	male.pred$mcmc.set <- e0.pred$mcmc.set
	for(item in names(male.pred$meta.changes))
		male.pred$mcmc.set$meta[[item]] <- male.pred$meta.changes[[item]]
	return(male.pred)
}

has.e0.jmale.prediction <- function(e0.pred) return(!is.null(e0.pred$joint.male))
