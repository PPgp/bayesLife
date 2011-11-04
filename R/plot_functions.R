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


e0.trajectories.plot <- function(e0.pred, country, pi=c(80, 95), 
								  nr.traj=NULL,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='Year', ylab='Life expectancy', main=NULL, ...
								  ) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	country <- get.country.object(country, e0.pred$mcmc.set$meta)
	e0.mtx <- e0.pred$mcmc.set$meta$e0.matrix
	T_end_c <- e0.pred$mcmc.set$meta$T.end.c
	e0.matrix.reconstructed <- get.e0.reconstructed(e0.pred$e0.matrix.reconstructed, 
									e0.pred$mcmc.set$meta)
	x1 <- as.integer(rownames(e0.matrix.reconstructed))
	x2 <- as.numeric(dimnames(e0.pred$quantiles)[[3]])

	lpart1 <- T_end_c[country$index]
	y1.part1 <- e0.mtx[1:T_end_c[country$index],country$index]
	y1.part2 <- NULL
	lpart2 <- nrow(e0.mtx) - T_end_c[country$index]
	if (lpart2 > 0) 
		y1.part2 <- e0.matrix.reconstructed[
			(T_end_c[country$index]+1):nrow(e0.matrix.reconstructed),country$index]
	trajectories <- bayesTFR:::get.trajectories(e0.pred, country$code, nr.traj)
	e0.median <- bayesTFR:::get.median.from.prediction(e0.pred, country$index, country$code)
	cqp <- list()
	ylim.loc <- c(min(trajectories$trajectories, y1.part1, y1.part2, e0.median, na.rm=TRUE), 
				  max(trajectories$trajectories, y1.part1, y1.part2, e0.median, na.rm=TRUE))
	for (i in 1:length(pi)) {
		cqp[[i]] <- bayesTFR:::get.traj.quantiles(e0.pred, country$index, 
					country$code, trajectories=trajectories$trajectories, 
												pi=pi[i])
		if (!is.null(cqp[[i]]) && is.null(ylim))
			ylim.loc <- c(min(ylim.loc[1], cqp[[i]], na.rm=TRUE), 
						  max(ylim.loc[2], cqp[[i]], na.rm=TRUE))
	}
	if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
	if(is.null(ylim)) ylim <- ylim.loc
	if(is.null(main)) main <- country$name
	# plot historical data: observed
	plot(x1[1:lpart1], y1.part1, type=type, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, 
			panel.first = grid(), ...
					)
	if(lpart2 > 0) {
		lines(x1[(lpart1+1): length(x1)], y1.part2, pch=2, type='b', col='green')
		lines(x1[lpart1:(lpart1+1)], c(y1.part1[lpart1], y1.part2[1]), col='green') # connection between the two parts
	}
	
	# plot trajectories
	if(!is.null(trajectories$trajectories)) { 
		for (i in 1:length(trajectories$index)) {
			lines(x2, trajectories$trajectories[,trajectories$index[i]], type='l', col='gray')
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
	legend <- c(legend, 'observed LifeExp')
	col <- c(col, 'black')
	lty <- c(lty, 1)
	pch <- c(rep(-1, length(legend)-1), 1)
	if(lpart2 > 0) {
		legend <- c(legend, 'imputed LifeExp')
		col <- c(col, 'green')
		lty <- c(lty, 1)
		pch <- c(pch, 2)
	}
	legend('topleft', legend=legend, lty=c(1,lty), bty='n', col=col, pch=pch)
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

.get.dlcurves <- function(x, mcmc.list, country.code, country.index, burnin, nr.curves) {
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
		errors <- rnorm(length(mcmc$meta$loessSD[,country.index]), 
						mean=0, sd=omegas*mcmc$meta$loessSD[,country.index])
		dl <- t(apply(dl.pars, 1, g.dl6, l=x, p1 = mcmc$meta$dl.p1, p2 = mcmc$meta$dl.p2))
        dlc <- rbind(dlc, dl+errors)
    }
	return (dlc)
}

e0.DLcurve.plot <- function (mcmc.list, country, burnin = NULL, pi = 80, e0.lim = c(20,90), 
    nr.curves = NULL, ylim = NULL, xlab = "e(0)", ylab = "5-year gains", 
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
    x <- seq(e0.lim[1], e0.lim[2], length=1000)
    dlc <- .get.dlcurves(x, mcmc.list, country$code, country$index, burnin, nr.curves)
    miny <- min(dlc)
    maxy <- max(dlc)
    thincurves <- bayesTFR:::get.thinning.index(nr.curves, dim(dlc)[1])
    ltype <- "l"
    if (thincurves$nr.points == 0) {
        ltype <- "n"
        thincurves$index <- 1
    }
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
    dl50 <- apply(dlc, 2, quantile, 0.5)
    lines(dl50 ~ x, col = "red", lwd = 2)
    lty <- 2:(length(pi) + 1)
    for (i in 1:length(pi)) {
        al <- (1 - pi[i]/100)/2
        dlpi <- apply(dlc, 2, quantile, c(al, 1 - al))
        lines(dlpi[1, ] ~ x, col = "red", lty = lty[i], 
            lwd = 2)
        lines(dlpi[2, ] ~ x, col = "red", lty = lty[i], 
            lwd = 2)
    }
    T.total <- nrow(meta$e0.matrix)
    incr <- diff(meta$e0.matrix[1:T.total, country$index])
    points(incr ~ meta$e0.matrix[1:(T.total - 
        1), country$index], pch = 19)
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

.map.main.default.bayesLife.prediction <- function(pred, ...) return('e0: quantile')

e0.map <- function(pred, ...) return(bayesTFR:::tfr.map(pred, ...))

e0.map.all <- function(pred, output.dir, output.type='png', e0.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, file.prefix='e0wrldmap_', ...) {
	bayesTFR:::bdem.map.all(pred=pred, output.dir=output.dir, type='e0', output.type=output.type, range=e0.range,
						nr.cats=nr.cats, same.scale=same.scale, quantile=quantile, file.prefix=file.prefix, ...)
}

e0.map.gvis <- function(pred, ...)
	bdem.map.gvis(pred, ...)
						
bdem.map.gvis.bayesLife.prediction <- function(pred, ...) {
	sex.label <- list(M='Male', F='Female')
	bayesTFR:::.do.gvis.bdem.map('e0', paste(sex.label[[pred$mcmc.set$meta$sex]], 'Life Expectancy'), pred, ...)
}

par.names.for.worldmap.bayesLife.prediction <- function(pred, ...) {
	return(e0.parameter.names.cs.extended())
}

get.data.for.worldmap.bayesLife.prediction <- function(pred, ...)
	return(bayesTFR:::get.data.for.worldmap.bayesTFR.prediction(pred, ...))
