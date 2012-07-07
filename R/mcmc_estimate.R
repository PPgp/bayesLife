e0.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {
	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	niter <- mcmc$iter
	meta <- mcmc$meta
	#lex.mtx <- meta$e0.matrix
	T <- nrow(meta$e0.matrix) - 1
	C <- meta$nr.countries
	mcmc$thin <- thin
	
	delta.sq <- mcmc$meta$delta^2
	#dlf <- matrix(NA, nrow=T, ncol=C)
	#psi.shape <- (sum(meta$T.end.c-1)-1)/2.
	
	if (start.iter > niter) return(mcmc)
	
	mcenv <- new.env() # Create an environment for the mcmc stuff in order to avoid 
						# copying of the mcmc list 
    for (item in names(mcmc)) mcenv[[item]] <- mcmc[[item]]
    Triangle.prop <- rep(0,4)
    
    dlf <- list()
    DLdata <- get.DLdata.for.estimation(meta, 1:C)

	for(iter in start.iter:niter) {
		if(verbose.iter > 0 && (iter %% verbose.iter == 0))
			cat('\nIteration:', iter, '--', date())
		unblock.gtk('bDem.e0mcmc')

		# Update Triangle, k, z using Gibbs sampler
		###########################################
		sum.Trkz.c <- rowSums(mcenv$Triangle.c)
		sum.Trkz.c <- c(sum.Trkz.c, sum(mcenv$k.c), sum(mcenv$z.c))
		lambdas <- c(mcenv$lambda, mcenv$lambda.k, mcenv$lambda.z)
		Tr.var <- 1./(1./delta.sq + C*lambdas)
		Tr.sd <- sqrt(Tr.var)
		Tr.mean <- (meta$a/delta.sq + sum.Trkz.c*lambdas)*Tr.var
		ntries <- 1
		while(ntries <= 50) {
			for (i in 1:4) {
				# Triangle is truncated normal in [0,100]
				Triangle.prop[i] <- rnorm.trunc(mean=Tr.mean[i], sd=Tr.sd[i], 
									low=meta$Triangle.prior.low[i], high=meta$Triangle.prior.up[i])
			}
			sT <- sum(Triangle.prop)
			# discard samples for which the sum(Triangle) is outside of given interval.
			if(sT <= meta$sumTriangle.lim[2] && sT >= meta$sumTriangle.lim[1]) break
			ntries <- ntries + 1
		}
		mcenv$Triangle <- Triangle.prop
		# k is truncated normal in [0,10]
		mcenv$k <- rnorm.trunc(mean=Tr.mean[5], sd=Tr.sd[5], low=meta$k.prior.low, high=meta$k.prior.up)
		# z is truncated normal in [0,1.15]
		mcenv$z <- rnorm.trunc(mean=Tr.mean[6], sd=Tr.sd[6], low=meta$z.prior.low, high=meta$z.prior.up)
		
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in 1:C) {
			Triangle.k.z.c.update(mcenv, country, DLdata=DLdata)
			dlf[[country]] <- g.dl6(c(mcenv$Triangle.c[,country], 
								mcenv$k.c[country], mcenv$z.c[country]), 
								DLdata[[country]]['e0',], meta$dl.p1, meta$dl.p2)
		}
		# Update omega 
		###########################################
		omega.update(mcenv, dlf, DLdata=DLdata)
		
		# Update lambdas using MH-algorithm
		###########################################
		lambdas.update(mcenv)

		################################################################### 
        # write samples simu/thin to disk
        ##################################################################
		mcenv$finished.iter <- mcenv$finished.iter+1
        mcenv$rng.state <- .Random.seed         
        if (iter %% thin == 0){
        	mcenv$length <- mcenv$length + 1
        	flush.buffer <- FALSE
            if (iter + 1 > niter) flush.buffer <- TRUE                
            store.e0.mcmc(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
         }
	}
	resmc <- as.list(mcenv)
	class(resmc) <- class(mcmc)
	return(resmc)
}

unblock.gtk <- function(...) bayesTFR:::unblock.gtk(...)

e0.mcmc.sampling.extra <- function(mcmc, mcmc.list, countries, posterior.sample,
											 iter=NULL, burnin=2000, 
											 verbose=FALSE, verbose.iter=100) {
	#run mcmc sampling for countries given by the index 'countries'
	niter <- mcmc$iter
	if (is.null(iter))
    	niter <- mcmc$length	
	# get values of the hyperparameters (sample from the posterior)
    hyperparameter.names <- e0.parameter.names()
    hyperparameters <- list()
    sampled.index <- sample(posterior.sample, niter, replace=TRUE)
    th.burnin <- bayesTFR:::get.thinned.burnin(mcmc, burnin)
    for (par in hyperparameter.names) {
    	hyperparameters[[par]] <- c()
    	for(mc in mcmc.list) {
    		if (bayesTFR:::no.traces.loaded(mc)  || th.burnin < mc$traces.burnin) {
    			traces <- bdem.parameter.traces(mc, par, burnin=th.burnin)
        	} else {
          		traces <- bayesTFR:::get.burned.tfr.traces(mc, par, th.burnin)
       		}
       		hyperparameters[[par]] <- rbind(hyperparameters[[par]], traces)
       	}
    	hyperparameters[[par]] <- hyperparameters[[par]][sampled.index,]
    }
	mcmc.orig <- mcmc
	mcenv <- new.env() # Create an environment for the mcmc stuff in order to avoid 
						# copying of the mcmc list 
	for (item in names(mcmc)) mcenv[[item]] <- mcmc[[item]]
	updated.var.names <- c('Triangle.c', 'k.c', 'z.c')
	DLdata <- get.DLdata.for.estimation(mcmc$meta, countries)

	for(iter in 1:niter) {
		if(verbose.iter > 0 && (iter %% verbose.iter == 0))
			cat('\nIteration:', iter, '--', date())
		unblock.gtk('bDem.e0mcmcExtra')
		# set hyperparameters for this iteration
        for (par in hyperparameter.names) {
        	if(is.null(dim(hyperparameters[[par]]))) {
        		mcenv[[par]] <- hyperparameters[[par]][iter]
        	} else {
        		mcenv[[par]] <- hyperparameters[[par]][iter,]
        	}
        }
						
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in countries) {
			Triangle.k.z.c.update(mcenv, country, DLdata=DLdata)
		}

		################################################################### 
        # write samples simu/thin to disk
        ##################################################################
         #update the original mcmc with the new values
         for(var in updated.var.names) {
         	mcmc.orig[[var]] <- mcenv[[var]]
         }
		 flush.buffer <- FALSE
         append <- TRUE
		 if (iter == 1) {
			append <- FALSE
			flush.buffer <- TRUE
		 } else {
			if (iter + 1 > niter) flush.buffer <- TRUE
		 }           
         store.e0.mcmc(mcmc.orig, append=append, flush.buffer=flush.buffer, countries=countries, 
         				verbose=verbose)
	}
	return(mcmc.orig)
}


get.DLdata.for.estimation <- function(meta, countries) {
	DLdata <- list()
    T.suppl.end <- if(!is.null(meta$suppl.data$e0.matrix)) nrow(meta$suppl.data$e0.matrix) else 0
    for(country in countries) {
    	idx <- which(!is.na(meta$d.ct[, country]))
    	DLdata[[country]] <- matrix(NA, nrow=3, ncol=length(idx), 
    							dimnames=list(c('e0', 'dct', 'loess'), 
    										rownames(meta$d.ct[idx,country])))
    	DLdata[[country]][1,] <- meta$e0.matrix[idx,country]
    	DLdata[[country]][2,] <- meta$d.ct[idx,country]
    	DLdata[[country]][3,] <- meta$loessSD[idx,country]
    }
    if(T.suppl.end > 0) {
    	for(country in 1:ncol(meta$suppl.data$e0.matrix)) {
    		cidx <- meta$suppl.data$index.to.all.countries[country]
    		if (!is.element(cidx, countries)) next
    		idx <- which(!is.na(meta$suppl.data$d.ct[, country]))
    		if(length(idx) <= 0) next
    		start.col <- ncol(DLdata[[cidx]]) + 1
    		DLdata[[cidx]] <- cbind(DLdata[[cidx]], 
    								matrix(NA, nrow=3, ncol=length(idx),
    										dimnames=list(rownames(DLdata[[cidx]]), 
    											rownames(meta$suppl.data$d.ct[idx,country]))))
    		end.col <- ncol(DLdata[[cidx]])
    		DLdata[[cidx]][1,start.col:end.col] <- meta$suppl.data$e0.matrix[idx,country]
       		DLdata[[cidx]][2,start.col:end.col] <- meta$suppl.data$d.ct[idx,country]
          	DLdata[[cidx]][3,start.col:end.col] <- meta$suppl.data$loessSD[idx,country]  
  		}
    }
	return(DLdata=DLdata)	
}

slice.sampling <- function(x0, fun, width,  ..., low, up) {
	gx0 <- fun(x0, ..., low=low, up=up)
	z <- gx0 - rexp(1) # slice S={x: z < gx0}
	u <- runif(2)
	L <- max(x0 - width*u[1], low)
	R <- min(L + width, up)
	i<-1
	maxit <- 50
	while (TRUE) {
		u <- runif(1)
		x1 <- L + u*(R-L)
		if(z < fun(x1,  ..., low=low, up=up)) return(x1)
		if (x1 < x0) L <- x1
		else R <- x1
		i <- i+1
		if(i>maxit) stop('Problem in slice sampling')
	}
}

Triangle.k.z.c.update <- function(mcmc, country, DLdata) {
	# Update Triangle pars using slice sampling
	sigmas <- 1/sqrt(mcmc$lambda)
	Triangle.c.low <- c(mcmc$meta$country.bounds$Triangle_1.c.prior.low[country], 
					  mcmc$meta$country.bounds$Triangle_2.c.prior.low[country],
					  mcmc$meta$country.bounds$Triangle_3.c.prior.low[country],
					  mcmc$meta$country.bounds$Triangle_4.c.prior.low[country])
	Triangle.c.up <- c(mcmc$meta$country.bounds$Triangle_1.c.prior.up[country], 
					  mcmc$meta$country.bounds$Triangle_2.c.prior.up[country],
					  mcmc$meta$country.bounds$Triangle_3.c.prior.up[country],
					  mcmc$meta$country.bounds$Triangle_4.c.prior.up[country])
	Triangle.prop <- rep(0,4)
	dlx <- c(mcmc$Triangle.c[,country], mcmc$k.c[country], mcmc$z.c[country])
	ntries <- 1
	while(ntries <= 50) {
		for (i in 1:4) {
			Triangle.prop[i] <- slice.sampling(mcmc$Triangle.c[i, country],
										logdensity.Triangle.k.z.c, mcmc$meta$Triangle.c.width[i], 
										mean=mcmc$Triangle[i], 
										sd=sigmas[i], dlx=dlx,
										low=Triangle.c.low[i], 
										up=Triangle.c.up[i], par.idx=i, 
										p1=mcmc$meta$dl.p1, p2=mcmc$meta$dl.p2, omega=mcmc$omega,
										DLdata=DLdata[[country]])
			dlx[i] <- Triangle.prop[i]
		}
		sT <- sum(Triangle.prop)
		if(sT <= mcmc$meta$sumTriangle.lim[2] && sT >= mcmc$meta$sumTriangle.lim[1]) break
		dlx <- c(mcmc$Triangle.c[,country], mcmc$k.c[country], mcmc$z.c[country])
		ntries <- ntries + 1
	}
	mcmc$Triangle.c[, country] <- Triangle.prop
	mcmc$k.c[country] <- slice.sampling(mcmc$k.c[country],
										logdensity.Triangle.k.z.c, mcmc$meta$k.c.width, mean=mcmc$k, 
										sd=1/sqrt(mcmc$lambda.k), dlx=dlx,
										low=mcmc$meta$country.bounds$k.c.prior.low[country], 
										up=mcmc$meta$country.bounds$k.c.prior.up[country], 
										par.idx=5, p1=mcmc$meta$dl.p1, p2=mcmc$meta$dl.p2, 
										omega=mcmc$omega, DLdata=DLdata[[country]])
	dlx[5] <- mcmc$k.c[country]

	mcmc$z.c[country] <- slice.sampling(mcmc$z.c[country],
										logdensity.Triangle.k.z.c, mcmc$meta$z.c.width, mean=mcmc$z, 
										sd=1/sqrt(mcmc$lambda.z), dlx=dlx,
										low=mcmc$meta$country.bounds$z.c.prior.low[country], 
										up=mcmc$meta$country.bounds$z.c.prior.up[country], 
										par.idx=6, p1=mcmc$meta$dl.p1, p2=mcmc$meta$dl.p2, 
										omega=mcmc$omega, DLdata=DLdata[[country]])
	return()
}

logdensity.Triangle.k.z.c <- function(x, mean, sd, dlx, low, up, par.idx, p1, p2, omega, DLdata) {
	logdens <- 0.0
	res <- .C("dologdensityTrianglekz", as.double(x), as.double(mean), as.double(sd), as.double(low), 
				as.double(up), as.integer(par.idx), as.double(dlx), 
				as.double(p1), as.double(p2), as.double(DLdata['e0',]), as.integer(ncol(DLdata)), 
				as.double(DLdata['dct',]), as.double(omega*DLdata['loess',]), logdens=logdens)
	return(res$logdens)	
}

lambdas.update <- function(mcmc) {
	# Update lambdas using MH-algorithm
	tau.sq <- mcmc$meta$tau^2
	for(i in 1:4) {
		prop <- proposal.lambda(mcmc$meta$nu, tau.sq[i], mcmc$Triangle[i], mcmc$Triangle.c[i,], 
								mcmc$meta$nr.countries)
		lpx0 <- logdensity.lambda(mcmc$lambda[i], mcmc, mcmc$meta$Triangle.prior.low[i], 
						mcmc$meta$Triangle.prior.up[i], mcmc$Triangle[i], mcmc$Triangle.c[i,], tau.sq[i])
		lpx1 <- logdensity.lambda(prop, mcmc, mcmc$meta$Triangle.prior.low[i], 
						mcmc$meta$Triangle.prior.up[i], mcmc$Triangle[i], mcmc$Triangle.c[i,], tau.sq[i])
		prob_accept <- min(exp(lpx1 - lpx0), 1)
		if (runif(1) < prob_accept){
			mcmc$lambda[i] <- prop	
		}
	}
	# update lambda.k 
	prop <- proposal.lambda(mcmc$meta$nu, tau.sq[5], mcmc$k, mcmc$k.c, mcmc$meta$nr.countries)
	lpx0 <- logdensity.lambda(mcmc$lambda.k, mcmc, mcmc$meta$k.prior.low, mcmc$meta$k.prior.up, 
					mcmc$k, mcmc$k.c, tau.sq[5])
	lpx1 <- logdensity.lambda(prop, mcmc, mcmc$meta$k.prior.low, mcmc$meta$k.prior.up, mcmc$k, mcmc$k.c, tau.sq[5])
	prob_accept <- min(exp(lpx1 - lpx0), 1)
	if (runif(1) < prob_accept) mcmc$lambda.k <- prop	
	# update lambda.z
	prop <- proposal.lambda(mcmc$meta$nu, tau.sq[6], mcmc$z, mcmc$z.c, mcmc$meta$nr.countries)
	lpx0 <- logdensity.lambda(mcmc$lambda.z, mcmc, mcmc$meta$z.prior.low, mcmc$meta$z.prior.up, mcmc$z, mcmc$z.c, tau.sq[6])
	lpx1 <- logdensity.lambda(prop, mcmc, mcmc$meta$z.prior.low, mcmc$meta$z.prior.up, mcmc$z, mcmc$z.c, tau.sq[6])
	prob_accept <- min(exp(lpx1 - lpx0), 1)
	if (runif(1) < prob_accept)	mcmc$lambda.z <- prop

	return()
}

logdensity.lambda <- function(lambda, mcmc, low, high, Triangle, Triangle.c, tau.sq) {
	nu <- mcmc$meta$nu
	return(log(dgamma(lambda, nu/2, rate=tau.sq)) + 
				sum(log(pmax(dnorm.trunc(Triangle.c, mean=Triangle, sd=1/sqrt(lambda), low, high), 
						.Machine$double.xmin))))
}

proposal.lambda <- function(nu, tau.sq, Triangle, Triangle.c, C) {
	return(rgamma(1, (nu+C)/2, rate=tau.sq+0.5*sum((Triangle.c-Triangle)^2)))
}

omega.update <- function(mcmc, dlf, DLdata) {
	omega.width <- 1		
	#print('Slice sampling for omega')
	mcmc$omega <- slice.sampling(mcmc$omega, logdensity.omega, omega.width, 
								mcmc=mcmc, dlf=dlf, DLdata=DLdata, low=0, up=10)
	return()
}

logdensity.omega <- function(x, mcmc, dlf, DLdata, low, up) {
	log.d.cDens <- 0
	for(country in 1:mcmc$meta$nr.countries) 
		log.d.cDens <- log.d.cDens + sum(log(pmax(dnorm(DLdata[[country]]['dct',], dlf[[country]], 
										sd=x*DLdata[[country]]['loess',]),.Machine$double.xmin)))
	return(log(dunif(x, low, up)) + log.d.cDens)
}

