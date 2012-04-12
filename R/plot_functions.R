e0.trajectories.plot.all <- function(e0.pred, 
									output.dir=file.path(getwd(), 'e0trajectories'),
									output.type="png", verbose=FALSE, ...) {
	# plots e0 trajectories for all countries
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.countries <- country.names(e0.pred$mcmc.set$meta)
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	for (country in all.countries) {
		country.obj <- get.country.object(country, e0.pred$mcmc.set$meta)
		if(verbose)
			cat('Creating e0 graph for', country, '(', country.obj$code, ')\n')

		do.call(output.type, list(file.path(output.dir, 
										paste('e0plot_c', country.obj$code, '.', postfix, sep=''))))
		e0.trajectories.plot(e0.pred, country=country.obj$code, ...)
		dev.off()
	}
	if(verbose)
		cat('\nTrajectory plots stored into', output.dir, '\n')
}

e0.gap.plot <- function(e0.pred, country, e0.pred2=NULL, pi=c(80, 95), nr.traj=0,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='Year', ylab='Gap in life expectancy', main=NULL, ...
								  ) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	country <- get.country.object(country, e0.pred$mcmc.set$meta)
	if(is.null(e0.pred2)) e0.pred2 <- get.e0.jmale.prediction(e0.pred)
	e0.mtx <- e0.pred$e0.matrix.reconstructed
	e0.mtx2 <- e0.pred2$e0.matrix.reconstructed
	T <- nrow(e0.mtx)
	x1 <- as.integer(rownames(e0.mtx))
	x2 <- as.numeric(dimnames(e0.pred$quantiles)[[3]])
	y1 <- e0.mtx[1:T,country$index]	- e0.mtx2[1:T,country$index]
	# get all trajectories for computing the quantiles
	trajobj <- bayesTFR:::get.trajectories(e0.pred, country$code)
	trajobj2 <- bayesTFR:::get.trajectories(e0.pred2, country$code)
	if(is.null(trajobj))
		stop('No trajectories available in the e0.pred object.')
	if(is.null(trajobj2))
		stop('No trajectories available in the e0.pred2 object.')
	if(dim(trajobj$trajectories)[2] != dim(trajobj2$trajectories)[2])
		stop('Trajectories in the two prediction ojects must be of the same shape.')
	thintraj <- bayesTFR:::get.thinning.index(nr.traj, dim(trajobj$trajectories)[2])

	y2all <- trajobj$trajectories - trajobj2$trajectories
	e0.median <- apply(y2all, 1, quantile, 0.5, na.rm = TRUE)
	y2 <- NA
	if (thintraj$nr.points > 0)
		y2 <- y2all[,thintraj$index]
	cqp <- list()
	al <- (1-pi/100)/2
	ylim.loc <- c(min(y2, y1, e0.median, na.rm=TRUE), 
				  max(y2, y1, e0.median, na.rm=TRUE))
	for (i in 1:length(pi)) {
		cqp[[i]] <- apply(y2all, 1, quantile, c(al[i], 1-al[i]), na.rm = TRUE)
		if (is.null(ylim))
			ylim.loc <- c(min(ylim.loc[1], cqp[[i]], na.rm=TRUE), 
						  max(ylim.loc[2], cqp[[i]], na.rm=TRUE))
	}
	if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
	if(is.null(ylim)) ylim <- ylim.loc
	if(is.null(main)) main <- country$name
	# plot historical data: observed
	plot(x1, y1, type=type, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
			panel.first = grid(), ...
					)	
	# plot trajectories
	if(thintraj$nr.points > 0) { 
		for (i in 1:thintraj$nr.points) {
			lines(x2, y2[,i], type='l', col='gray')
		}
	}
	# plot median	
	lines(x2, e0.median, type='l', col='red', lwd=2) 
	# plot given CIs
	lty <- 2:(length(pi)+1)
	for (i in 1:length(pi)) {
		if (!is.null(cqp[[i]])) {
			lines(x2, cqp[[i]][1,], type='l', col='red', lty=lty[i], lwd=2)
			lines(x2, cqp[[i]][2,], type='l', col='red', lty=lty[i], lwd=2)
		}
	}
	legend <- c('median', paste('PI', pi))
	col <- rep('red', length(lty)+1)
	legend <- c(legend, 'observed gap')
	col <- c(col, 'black')
	lty <- c(lty, 1)
	pch <- c(rep(-1, length(legend)-1), 1)
	legend('topleft', legend=legend, lty=c(1,lty), bty='n', col=col, pch=pch)
}

e0.trajectories.plot <- function(e0.pred, country, pi=c(80, 95), both.sexes=FALSE,
								  nr.traj=NULL, typical.trajectory=FALSE,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='Year', ylab='Life expectancy at birth', main=NULL, 
								  lwd=c(2,2,2,2,1), col=c('black', 'green', 'red', 'red', 'gray'),
								  col2=c('gray39', 'greenyellow', 'hotpink', 'hotpink', 'gray'),
								  show.legend=TRUE, add=FALSE, ...
								  ) {
	# lwd/col is a vector of 5 line widths/colors for: 
	#	1. observed data, 2. imputed missing data, 3. median, 4. quantiles, 5. trajectories

	lowerize <- function(str) { # taken from the cwhmisc package
		ff <- function(x) paste(lapply(unlist(strsplit(x, NULL)),lower),collapse="") 
		capply(str,ff)
	}
	capply <- function (str, ff, ...) { # taken from the cwhmisc package
    	x <- strsplit(str, NULL)
    	y <- lapply(x, ff, ...)
    	sapply(y, paste, collapse = "")
	}
	lower <- function (char) { # taken from the cwhmisc package
    	if (any(ind <- LETTERS == char)) letters[ind]
    	else char
	}
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	if(length(lwd) < 5) {
		llwd <- length(lwd)
		lwd <- rep(lwd, 5)
		lwd[(llwd+1):5] <- c(2,2,2,2,1)[(llwd+1):5]
	}
	if(length(col) < 5) {
		lcol <- length(col)
		col <- rep(col, 5)
		col[(lcol+1):5] <- c('black', 'green', 'red', 'red', 'gray')[(lcol+1):5]
	}
	country <- get.country.object(country, e0.pred$mcmc.set$meta)
	pred <- list(e0.pred)
	plotcols <- list(col)
	do.both.sexes <- FALSE
	if(both.sexes == TRUE) {
		do.both.sexes <- TRUE
		if(e0.pred$mcmc.set$meta$sex != 'F') {
			warnings('If both.sexes is TRUE, the given prediction obect must be Female prediction, but is ',  
					get.sex.label(e0.pred$mcmc.set$meta))
			do.both.sexes <- FALSE
		}
		if(!has.e0.jmale.prediction(e0.pred)) {
			warnings('A male prediction does not exist for the given prediction object.')
			do.both.sexes <- FALSE
		}
		if(do.both.sexes) {
			pred <- c(pred, list(get.e0.jmale.prediction(e0.pred)))
			if(length(col2) < 5) {
				lcol <- length(col2)
				col2 <- rep(col2, 5)
				col2[(lcol+1):5] <- c('gray39', 'greenyellow', 'hotpink', 'hotpink', 'gray')[(lcol+1):5]
			}
			plotcols <- c(list(col2), plotcols)
			if(missing(nr.traj)) nr.traj <- 0
			if(missing(pi)) pi <- 95
		}
	}
	lty.all <- cols.all <- legend.all <- pch.all <- lwd.all <- c()
	plot.data <- list()
	ylim.loc <- c(100,0)
	for(ipred in 1:length(pred)) {
		e0pred <- pred[[ipred]]
		meta <- e0pred$mcmc.set$meta
	
		e0.mtx <- meta$e0.matrix
		e0.matrix.reconstructed <- get.e0.reconstructed(e0pred$e0.matrix.reconstructed, meta)
	
		x1 <- as.integer(rownames(e0.matrix.reconstructed))
		x2 <- as.numeric(dimnames(e0pred$quantiles)[[3]])

		if(!is.null(meta$T.end.c)) lpart1 <- meta$T.end.c[country$index]
		else lpart1 <- meta$Tc.index[[country$index]][length(meta$Tc.index[[country$index]])]
	
	
		y1.part1 <- e0.mtx[1:lpart1,country$index]
		y1.part2 <- NULL
		lpart2 <- nrow(e0.mtx) - lpart1
		if (lpart2 > 0) # imputed values
			y1.part2 <- e0.matrix.reconstructed[
				(lpart1+1):nrow(e0.matrix.reconstructed),country$index]
			
		if(!is.null(meta$suppl.data$e0.matrix)) {
    		supp.c.idx <- meta$suppl.data$index.from.all.countries[country$index]
    		if(!is.na(supp.c.idx)) {
    			suppl.data.idx <- which(!is.na(meta$suppl.data$e0.matrix[,supp.c.idx]))
    			lpart1 <- lpart1 + length(suppl.data.idx)
    			y1.part1 <- c(meta$suppl.data$e0.matrix[suppl.data.idx, supp.c.idx], y1.part1)
    			x1 <- c(as.integer(rownames(meta$suppl.data$e0.matrix[suppl.data.idx,,drop=FALSE])), x1)
 			}
    	}
    	plot.data[[ipred]] <- list(obs.x=x1[1:lpart1], obs.y=y1.part1, pred.x=x2)
    	ylim.loc <- c(min(ylim.loc[1], plot.data[[ipred]]$obs.y), max(ylim.loc[2], plot.data[[ipred]]$obs.y))
    	if(lpart2 > 0) {
    		plot.data[[ipred]]$rec.x <- x1[(lpart1+1): length(x1)]
    		plot.data[[ipred]]$rec.y <- y1.part2
    		ylim.loc <- c(min(ylim.loc[1], plot.data[[ipred]]$rec.y), max(ylim.loc[2], plot.data[[ipred]]$rec.y))
    	}
	}
	for(ipred in 1:length(pred)) {
		e0pred <- pred[[ipred]]
		this.col <- plotcols[[ipred]]
		meta <- e0pred$mcmc.set$meta
		trajectories <- bayesTFR:::get.trajectories(e0pred, country$code, nr.traj, typical.trajectory=typical.trajectory)
		e0.median <- bayesTFR:::get.median.from.prediction(e0pred, country$index, country$code)
		cqp <- list()
		if(ipred > 1) add <- TRUE
		if(!add)
			ylim.loc <- c(min(if (!is.null(trajectories$trajectories))
							trajectories$trajectories[,trajectories$index]
					  	else NULL, 
					  	ylim.loc[1], e0.median, na.rm=TRUE), 
				  	max(if (!is.null(trajectories$trajectories))
							trajectories$trajectories[,trajectories$index]
					  else NULL, 
					  ylim.loc[2], e0.median, na.rm=TRUE))
		if(length(pi) > 0) {
			for (i in 1:length(pi)) {
				cqp[[i]] <- bayesTFR:::get.traj.quantiles(e0pred, country$index, 
								country$code, trajectories=trajectories$trajectories, pi=pi[i])
				if (!add && !is.null(cqp[[i]]) && is.null(ylim))
					ylim.loc <- c(min(ylim.loc[1], cqp[[i]], na.rm=TRUE), 
						  		max(ylim.loc[2], cqp[[i]], na.rm=TRUE))
			}
		}
		if (!add) {
			if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
			if(is.null(ylim)) ylim <- ylim.loc
			if(is.null(main)) {
				main <- country$name
				if(!do.both.sexes)
					main <- paste(main, '-', get.sex.label(meta))
			}
		    # plot historical data: observed
			plot(plot.data[[ipred]]$obs.x, plot.data[[ipred]]$obs.y, type=type, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
				panel.first = grid(), lwd=lwd[1], col=this.col[1], ...
					)
		} else # add to an existing plot
			points(plot.data[[ipred]]$obs.x, plot.data[[ipred]]$obs.y, type=type, lwd=lwd[1], col=this.col[1], ...
					)
		if(!is.null(plot.data[[ipred]]$rec.x)) { # plot reconstructed missing data
			lines(plot.data[[ipred]]$rec.x, plot.data[[ipred]]$rec.y, pch=2, type='b', col=this.col[2], lwd=lwd[2])
			lines(c(plot.data[[ipred]]$obs.x[length(plot.data[[ipred]]$obs.x)], plot.data[[ipred]]$rec.x[1]), 
				c(plot.data[[ipred]]$obs.y[length(plot.data[[ipred]]$obs.y)], plot.data[[ipred]]$rec.y[1]), 
				col=this.col[2], lwd=lwd[2]) # connection between the two parts
		}
	
		# plot trajectories
		if(!is.null(trajectories$trajectories)) { 
			for (i in 1:length(trajectories$index)) {
				lines(plot.data[[ipred]]$pred.x, trajectories$trajectories[,trajectories$index[i]], type='l', 
					col=this.col[5], lwd=lwd[5])
			}
		}
		# plot median	
		lines(plot.data[[ipred]]$pred.x, e0.median, type='l', col=this.col[3], lwd=lwd[3])
		legend <- 'median'
		if(do.both.sexes) legend <- paste(lowerize(get.sex.label(meta)), legend)
		lty <- 1
		# plot given CIs
		if(length(pi) > 0) {
			tlty <- 2:(length(pi)+1)
			for (i in 1:length(pi)) {
				if (!is.null(cqp[[i]])) {
					lines(plot.data[[ipred]]$pred.x, cqp[[i]][1,], type='l', col=this.col[4], lty=tlty[i], lwd=lwd[4])
					lines(plot.data[[ipred]]$pred.x, cqp[[i]][2,], type='l', col=this.col[4], lty=tlty[i], lwd=lwd[4])
				}
			}
			legend <- c(legend, paste(pi, '% PI ', if(do.both.sexes) lowerize(get.sex.label(meta)) else '', sep=''))
			lty <- c(lty, tlty)
		}
		legend <- c(legend, paste('observed', if(do.both.sexes) paste(lowerize(get.sex.label(meta)), 'e0') else 'e0'))
		lty <- c(lty, 1)
		pch <- c(rep(-1, length(legend)-1), 1)
		lwds <- c(lwd[3], rep(lwd[4], length(pi)), lwd[1])
		cols <- c(this.col[3], rep(this.col[4], length(pi)), this.col[1])
		if(lpart2 > 0) {
			legend <- c(legend, paste('imputed', if(do.both.sexes) paste(lowerize(get.sex.label(meta)), 'e0') else 'e0'))
			cols <- c(cols, this.col[2])
			lty <- c(lty, 1)
			pch <- c(pch, 2)
			lwds <- c(lwds, lwd[2])
		}
		lty.all <- c(lty.all, lty)
		legend.all <- c(legend.all, legend)
		cols.all <- c(cols.all, cols)
		pch.all <- c(pch.all, pch)
		lwd.all <- c(lwd.all, lwds)
	}
	if(show.legend)
		legend('topleft', legend=legend.all, lty=lty.all, bty='n', col=cols.all, pch=pch.all, lwd=lwd.all)
	#abline(h=1, lty=3)
	#abline(h=1.5, lty=3)
	#abline(h=2.1, lty=3)
}

e0.trajectories.table <- function(e0.pred, ...) {
	return(tfr.trajectories.table(e0.pred, half.child.variant = FALSE, ...))
}

e0.DLcurve.plot.all <- function (mcmc.list = NULL, sim.dir = NULL, 
					output.dir = file.path(getwd(), "DLcurves"), 
					output.type="png",
					burnin = NULL, verbose = FALSE,  ...) {
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	if(is.null(mcmc.list)) mcmc.list <- get.e0.mcmc(sim.dir=sim.dir, verbose=verbose, burnin=burnin)
	mc <- get.mcmc.list(mcmc.list)
	meta <- mc[[1]]$meta
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'

    for (country in 1:meta$nr.countries) {
        country.obj <- get.country.object(country, meta, index=TRUE)
        if (verbose) 
            cat("Creating DL graph for", country.obj$name, '(', country.obj$code, ')\n')
        do.call(output.type, list(file.path(output.dir, 
										paste('DLplot_c', country.obj$code, '.', postfix, sep=''))))
        e0.DLcurve.plot(mcmc.list = mcmc.list, country = country.obj$code, 
            burnin = burnin, ...)
        dev.off()
    }
    if (verbose) 
        cat("\nDL plots stored into", output.dir, "\n")
}

e0.get.dlcurves <- function(x, mcmc.list, country.code, country.index, burnin, nr.curves, predictive.distr=FALSE) {
	dlc <- c()
    nr.curves.from.mc <- if (!is.null(nr.curves)) ceiling(max(nr.curves, 2000)/length(mcmc.list))
    						else NULL
    postfix <- paste('_c', country.code, sep='')
    dl.par.names <- c(paste('Triangle.c_1', postfix,sep=''),
						paste('Triangle.c_2', postfix,sep=''), 
						paste('Triangle.c_3', postfix,sep=''), 
						paste('Triangle.c_4', postfix,sep=''), 
						paste('k.c', postfix,sep=''),
						paste('z.c', postfix,sep=''))
	T <- length(mcmc.list[[1]]$meta$loessSD[,country.index])
	if(predictive.distr) {
		data(loess_sd)
		loessSD <- loess.lookup(x)
		if(!is.null(mcmc.list[[1]]$meta$constant.variance) && mcmc.list[[1]]$meta$constant.variance)
			loessSD[] <- 1
	}
    for (mcmc in mcmc.list) {
    	th.burnin <- bayesTFR:::get.thinned.burnin(mcmc,burnin)
    	thincurves.mc <- bayesTFR:::get.thinning.index(nr.curves.from.mc, 
            all.points=mcmc$length - th.burnin)
        traces <- load.e0.parameter.traces.cs(mcmc, country.code, 
        						burnin=th.burnin, 
								thinning.index=thincurves.mc$index)
		dl.pars <- traces[,dl.par.names, drop=FALSE]
		omegas <- load.e0.parameter.traces(mcmc, par.names='omega', burnin=th.burnin, 
								thinning.index=thincurves.mc$index)
		dl <- t(apply(dl.pars, 1, g.dl6, l=x, p1 = mcmc$meta$dl.p1, p2 = mcmc$meta$dl.p2))
		if(predictive.distr) {
			errors <- matrix(NA, nrow=dim(dl)[1], ncol=dim(dl)[2])
			n <- ncol(errors)
			for(i in 1:nrow(errors))
				errors[i,] <- rnorm(n, mean=0, sd=omegas[i]*loessSD)
        	dlc <- rbind(dlc, dl+errors)
        } else dlc <- rbind(dlc, dl)
    }
	return (dlc)
}

e0.DLcurve.plot <- function (mcmc.list, country, burnin = NULL, pi = 80, e0.lim = NULL, 
    nr.curves = 20, predictive.distr=FALSE, ylim = NULL, xlab = "e(0)", ylab = "5-year gains", 
    main = NULL, ...
    ) 
{	
	if(class(mcmc.list) == 'bayesLife.prediction') {
		if(!is.null(burnin) && burnin != mcmc.list$burnin)
			warning('Prediction was generated with different burnin. Burnin set to ', mcmc.list$burnin)
		burnin <- 0 # because burnin was already cut of the traces
	}
	if(is.null(burnin)) burnin <- 0
    mcmc.list <- get.mcmc.list(mcmc.list)
    meta <- mcmc.list[[1]]$meta
    country <- get.country.object(country, meta)
    data.idx <- which(!is.na(meta$d.ct[,country$index]))
    incr <- meta$d.ct[data.idx, country$index]
    obs.data <- meta$e0.matrix[data.idx, country$index]
    if(!is.null(meta$suppl.data$e0.matrix)) {
    	supp.c.idx <- which(is.element(meta$suppl.data$index.to.all.countries, country$index))
    	if(length(supp.c.idx) > 0) {
    		suppl.data.idx <- which(!is.na(meta$suppl.data$d.ct[,supp.c.idx]))
    		obs.data <- c(meta$suppl.data$e0.matrix[suppl.data.idx, supp.c.idx], obs.data)
    		incr <- c(meta$suppl.data$d.ct[suppl.data.idx, supp.c.idx], incr)    		}
    }
    if (is.null(e0.lim)) e0.lim <- c(min(40, obs.data), max(90, obs.data))
    x <- seq(e0.lim[1], e0.lim[2], length=1000)
    dlc <- e0.get.dlcurves(x, mcmc.list, country$code, country$index, burnin, nr.curves, 
    						predictive.distr=predictive.distr)
    thincurves <- bayesTFR:::get.thinning.index(nr.curves, dim(dlc)[1])
    ltype <- "l"
    if (thincurves$nr.points == 0) {
        ltype <- "n"
        thincurves$index <- 1
    }
    dl50 <- apply(dlc, 2, quantile, 0.5)
    dlpi <- array(0, c(length(pi), 2, ncol(dlc)))
    for (i in 1:length(pi)) {
        al <- (1 - pi[i]/100)/2
        dlpi[i,,] <- apply(dlc, 2, quantile, c(al, 1 - al))
    }
    miny <- min(dlc[thincurves$index, ], dlpi, incr)
    maxy <- max(dlc[thincurves$index, ], dlpi, incr)
    if (is.null(main)) main <- country$name
    if (is.null(ylim)) ylim <- c(miny, maxy)

    plot(dlc[thincurves$index[1], ] ~ x, col = "grey", 
        type = ltype, xlim = c(min(x), max(x)), 
        ylim = ylim, ylab = ylab, xlab = xlab, main = main, ...
        )
    if (thincurves$nr.points > 1) {
        for (i in 2:thincurves$nr.points) {
            lines(dlc[thincurves$index[i], ] ~ x, col = "grey")
        }
    }
    lines(dl50 ~ x, col = "red", lwd = 2)
    lty <- 2:(length(pi) + 1)
    for (i in 1:length(pi)) {
        lines(dlpi[i,1, ] ~ x, col = "red", lty = lty[i], lwd = 2)
        lines(dlpi[i,2, ] ~ x, col = "red", lty = lty[i], lwd = 2)
    }

    points(incr ~ obs.data, pch = 19)
    legend("topright", legend = c("median", paste("PI", pi)), 
        lty = c(1, lty), bty = "n", col = "red")
}

e0.partraces.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesLife.output'),
								chain.ids=NULL, par.names=e0.parameter.names(), 
                                nr.points=NULL, dev.ncol=5, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.e0.mcmc(sim.dir, low.memory=low.memory)
	bayesTFR:::do.plot.tfr.partraces(mcmc.list, load.e0.parameter.traces, chain.ids=chain.ids, 
        						nr.points=nr.points, par.names=par.names, dev.ncol=dev.ncol, ...)
}

e0.partraces.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesLife.output'),
									chain.ids=NULL, par.names=e0.parameter.names.cs(),
                                    nr.points=NULL, dev.ncol=3, low.memory=TRUE, ...) {

	if (is.null(mcmc.list))
		mcmc.list <- get.e0.mcmc(sim.dir, low.memory=low.memory)
	mcmc.list <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.list[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	bayesTFR:::do.plot.tfr.partraces(mcmc.list, load.e0.parameter.traces.cs, 
		main.postfix=paste('(',country.obj$name,')', sep=''), chain.ids=chain.ids, nr.points=nr.points, 
		country=country.obj$code, par.names=par.names, dev.ncol=dev.ncol, ...)
}

e0.pardensity.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesLife.output'), 
									chain.ids=NULL, par.names=e0.parameter.names(), 
									burnin=NULL, dev.ncol=5, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.e0.mcmc(sim.dir, low.memory=low.memory)
	bayesTFR:::do.plot.tfr.pardensity(mcmc.list, get.e0.parameter.traces, chain.ids=chain.ids, par.names=par.names,
			par.names.ext=bayesTFR:::get.full.par.names(par.names, 
						e0.get.all.parameter.names.extended()),
			burnin=burnin, dev.ncol=dev.ncol, ...)
}

e0.pardensity.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesLife.output'), 
									chain.ids=NULL, par.names=e0.parameter.names.cs(), 
									burnin=NULL, dev.ncol=3, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.e0.mcmc(sim.dir, low.memory=low.memory)
	mcmc.l <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.l[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	bayesTFR:::do.plot.tfr.pardensity(mcmc.list, get.e0.parameter.traces.cs, chain.ids=chain.ids, par.names=par.names,
		par.names.ext=bayesTFR:::get.full.par.names.cs(par.names, 
								e0.parameter.names.cs.extended(country.obj$code)),
		main.postfix=paste('(',country.obj$name,')', sep=''),
		func.args=list(country.obj=country.obj),
		burnin=burnin, dev.ncol=dev.ncol, ...)
}

.get.gamma.pars.bayesLife.prediction <- function(pred, ...) {
	# estimated by
	# library(MASS)
	# data <- pred$e0.matrix.reconstructed[12,]
	# gd <- fitdistr(data-min(data)+0.05, densfun='gamma')
	# min(data) is 43.86
	return(list(gamma.pars=list(shape=7.65, rate=0.29), gamma.shift=43.86-0.05, min.value=43.8, 
					max.value=120))
}

get.e0.map.parameters <- function(pred, e0.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, ...) {
	return(bayesTFR:::get.tfr.map.parameters(pred, e0.range, nr.cats=nr.cats, same.scale=same.scale,
							quantile=quantile, ...))
}

.map.main.default.bayesLife.prediction <- function(pred, ...) 
	return(paste(get.sex.label(pred$mcmc.set$meta), 'e0: quantile'))

e0.map <- function(pred, ...) {
	return(bayesTFR:::tfr.map(pred, ...))
}
e0.map.all <- function(pred, output.dir, output.type='png', e0.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, file.prefix='e0wrldmap_', ...) {
	bayesTFR:::bdem.map.all(pred=pred, output.dir=output.dir, type='e0', output.type=output.type, range=e0.range,
						nr.cats=nr.cats, same.scale=same.scale, quantile=quantile, file.prefix=file.prefix, ...)
}

e0.map.gvis <- function(pred, ...)
	bdem.map.gvis(pred, ...)
						
bdem.map.gvis.bayesLife.prediction <- function(pred, ...) {
	bayesTFR:::.do.gvis.bdem.map('e0', paste(get.sex.label(pred$mcmc.set$meta), 'Life Expectancy'), pred, ...)
}

par.names.for.worldmap.bayesLife.prediction <- function(pred, ...) {
	return(e0.parameter.names.cs.extended())
}

get.data.for.worldmap.bayesLife.prediction <- function(pred, ...)
	return(bayesTFR:::get.data.for.worldmap.bayesTFR.prediction(pred, ...))

get.sex.label <- function(meta) return(list(M='Male', F='Female')[[meta$sex]])