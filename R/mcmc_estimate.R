e0.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE) {
	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	niter <- mcmc$iter
	meta <- mcmc$meta
	#lex.mtx <- meta$e0.matrix
	T <- nrow(meta$e0.matrix) - 1
	C <- meta$nr.countries
	mcmc$thin <- thin
	
	delta.sq <- mcmc$meta$delta^2
	#dlf <- matrix(NA, nrow=T, ncol=C)
	psi.shape <- (sum(meta$T.end.c-1)-1)/2.
	
	if (start.iter > niter) return(mcmc)
	
	mcenv <- new.env() # Create an environment for the mcmc stuff in order to avoid 
						# copying of the mcmc list 
    for (item in names(mcmc)) mcenv[[item]] <- mcmc[[item]]
    Triangle.prop <- rep(0,4)
    
    dlf <- lex.mtx <- dct <- loess.sd <- list()
    T.suppl.end <- if(!is.null(meta$suppl.data$e0.matrix)) nrow(meta$suppl.data$e0.matrix) else 0
    for(country in 1:C) {
    	lex.mtx[[country]] <- meta$e0.matrix[1:meta$T.end.c[country],country]
    	dct[[country]] <- meta$d.ct[1:(meta$T.end.c[country]-1),country]
    	loess.sd[[country]] <- meta$loessSD[1:(meta$T.end.c[country]-1),country]
    }
    if(T.suppl.end > 0) {
    	for(country in 1:ncol(meta$suppl.data$e0.matrix)) {
    		idx <- meta$suppl.data$index.to.all.countries[country]
    		lex.mtx[[idx]] <- c(meta$suppl.data$e0.matrix[meta$suppl.data$T.start.c[country]:T.suppl.end,country],
    								lex.mtx[[idx]])
    		dct[[idx]] <- c(meta$suppl.data$d.ct[meta$suppl.data$T.start.c[country]:T.suppl.end,country],
    							dct[[idx]])
    		loess.sd[[idx]] <- c(meta$suppl.data$loessSD[meta$suppl.data$T.start.c[country]:T.suppl.end,country],
    								loess.sd[[idx]])
    	}
    }
	for(iter in start.iter:niter) {
		if(verbose || (iter %% 10 == 0))
			cat('\nIteration:', iter, '--', date())
				
		# Update Triangle, k, z using Gibbs sampler
		###########################################
		sum.Trkz.c <- rowSums(mcenv$Triangle.c)
		sum.Trkz.c <- c(sum.Trkz.c, sum(mcenv$k.c), sum(mcenv$z.c))
		lambdas <- c(mcenv$lambda, mcenv$lambda.k, mcenv$lambda.z)
		Tr.var <- 1./(1./delta.sq + C*lambdas)
		Tr.sd <- sqrt(Tr.var)
		Tr.mean <- (meta$a/delta.sq + sum.Trkz.c*lambdas)*Tr.var
		while(TRUE) {
			for (i in 1:4) {
				# Triangle is truncated normal in [0,100]
				Triangle.prop[i] <- rnorm.trunc(mean=Tr.mean[i], sd=Tr.sd[i], 
									low=meta$Triangle.prior.low[i], high=meta$Triangle.prior.up[i])
			}
			if(sum(Triangle.prop) < 110) break # discard samples for which the sum(Triangle) > 110.
		}
		mcenv$Triangle <- Triangle.prop
		# k is truncated normal in [0,10]
		mcenv$k <- rnorm.trunc(mean=Tr.mean[5], sd=Tr.sd[5], low=meta$k.prior.low, high=meta$k.prior.up)
		# z is truncated normal in [0,1.15]
		if(meta$vary.z.over.countries)
			mcenv$z <- rnorm.trunc(mean=Tr.mean[6], sd=Tr.sd[6], low=meta$z.prior.low, high=meta$z.prior.up)
		else {
			mcenv$z <- z.gibbs(mcenv, iter)
		}
		
		#ratesum <- 0		
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in 1:C) {
			Triangle.k.z.c.update(mcenv, country, lex.mtx=lex.mtx[[country]], dct=dct, loess.sd=loess.sd)
			#idx <- 1:(meta$T.end.c[country]-1)
			dlf[[country]] <- g.dl6(c(mcenv$Triangle.c[,country], 
								mcenv$k.c[country], mcenv$z.c[country]), 
								lex.mtx[[country]][2:length(lex.mtx[[country]])], meta$dl.p1, meta$dl.p2)
			#ratesum <- ratesum + ((meta$d.ct[idx,country]-dlf[idx,country])^2)/(meta$loessSD[idx,country]^2)
		}
		# Update omega using Gibbs sampler
		###########################################
		#psi <- rgamma.ltrunc(shape=psi.shape, rate=0.5*ratesum, low=0.01)
		#mcmc$omega <- 1/sqrt(psi)
		omega.update(mcenv, dlf, dct=dct, loess.sd=loess.sd)
		
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

e0.mcmc.sampling.extra <- function(mcmc, mcmc.list, countries, posterior.sample,
											 iter=NULL, burnin=2000, verbose=FALSE) {
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

	for(iter in 1:niter) {
		if(verbose || (iter %% 10 == 0))
			cat('\nIteration:', iter, '--', date())
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
			Triangle.k.z.c.update(mcenv, country)
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

Triangle.k.z.c.update <- function(mcmc, country, lex.mtx, dct, loess.sd) {
	# Update Triangle pars using slice sampling
	sigmas <- 1/sqrt(mcmc$lambda)
	Triangle.c.width <- mcmc$meta$Triangle.c.width
	k.c.width <- mcmc$meta$k.c.width
	z.c.width <- mcmc$meta$z.c.width
	idx <- 1:(mcmc$meta$T.end.c[country]-1)
	lex <- lex.mtx[2:length(lex.mtx)]
	Triangle.prop <- rep(0,4)
	while(TRUE) {
		for (i in 1:4) {
			Triangle.prop[i] <- slice.sampling(mcmc$Triangle.c[i, country],
										logdensity.Triangle.k.z.c, Triangle.c.width[i], 
										mean=mcmc$Triangle[i], 
										sd=sigmas[i], low=mcmc$meta$Triangle.c.prior.low[i], 
										up=mcmc$meta$Triangle.c.prior.up[i], par.idx=i, 
										mcmc=mcmc, 
										country=country, lex=lex, idx=idx,
										dct=dct, loess.sd=loess.sd)
		}
		if(sum(Triangle.prop) <= 110) break
	}
	mcmc$Triangle.c[, country] <- Triangle.prop
	mcmc$k.c[country] <- slice.sampling(mcmc$k.c[country],
										logdensity.Triangle.k.z.c, k.c.width, mean=mcmc$k, 
										sd=1/sqrt(mcmc$lambda.k),
										low=mcmc$meta$k.c.prior.low, up=mcmc$meta$k.c.prior.up, 
										par.idx=5, mcmc=mcmc, 
										country=country, lex=lex, idx=idx,
										dct=dct, loess.sd=loess.sd)
	if(mcmc$meta$vary.z.over.countries)
		mcmc$z.c[country] <- slice.sampling(mcmc$z.c[country],
										logdensity.Triangle.k.z.c, z.c.width, mean=mcmc$z, 
										sd=1/sqrt(mcmc$lambda.z),
										low=mcmc$meta$z.c.prior.low, up=mcmc$meta$z.c.prior.up, 
										par.idx=6, mcmc=mcmc, 
										country=country, lex=lex, idx=idx,
										dct=dct, loess.sd=loess.sd)
	else mcmc$z.c[country] <- mcmc$z
	return()
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
	if (runif(1) < prob_accept){
		mcmc$lambda.k <- prop	
	}
	# update lambda.z
	if(mcmc$meta$vary.z.over.countries) {
		prop <- proposal.lambda(mcmc$meta$nu, tau.sq[6], mcmc$z, mcmc$z.c, mcmc$meta$nr.countries)
		lpx0 <- logdensity.lambda(mcmc$lambda.z, mcmc, mcmc$meta$z.prior.low, mcmc$meta$z.prior.up, mcmc$z, mcmc$z.c, tau.sq[6])
		lpx1 <- logdensity.lambda(prop, mcmc, mcmc$meta$z.prior.low, mcmc$meta$z.prior.up, mcmc$z, mcmc$z.c, tau.sq[6])
		prob_accept <- min(exp(lpx1 - lpx0), 1)
		if (runif(1) < prob_accept){
			mcmc$lambda.z <- prop	
		}
	} else mcmc$lambda.z <- 1
	return()
}

logdensity.lambda <- function(lambda, mcmc, low, high, Triangle, Triangle.c, tau.sq, ...) {
	nu <- mcmc$meta$nu
	return(log(dgamma(lambda, nu/2, rate=tau.sq)) + 
				sum(log(dnorm.trunc(Triangle.c, mean=Triangle, sd=1/sqrt(lambda), low, high))))
}

proposal.lambda <- function(nu, tau.sq, Triangle, Triangle.c, C) {
	return(rgamma(1, (nu+C)/2, rate=tau.sq+0.5*sum((Triangle.c-Triangle)^2)))
}

omega.update <- function(mcmc, dlf, dct, loess.sd) {
	omega.width <- 1		
	#print('Slice sampling for omega')
	mcmc$omega <- slice.sampling(mcmc$omega, logdensity.omega, omega.width, 
								mcmc=mcmc, dlf=dlf, dct=dct, loess.sd=loess.sd, low=0, up=10)
	return()
}

logdensity.omega <- function(x, mcmc, dlf, dct, loess.sd, low, up) {
	log.d.cDens <- 0
	for(country in 1:mcmc$meta$nr.countries) 
		log.d.cDens <- log.d.cDens + sum(log(pmax(dnorm(dct[[country]], dlf[[country]], 
										sd=x*loess.sd[[country]]),1e-100)))
	return(log(dunif(x, low, up)) + log.d.cDens)
}


logdensity.Triangle.k.z.c <- function(x, mean, sd, low, up, par.idx, mcmc, country, lex, idx, dct, loess.sd) {
	dlx <- c(Triangle_1=mcmc$Triangle.c[1,country], Triangle_2=mcmc$Triangle.c[2,country],
			Triangle_3=mcmc$Triangle.c[3,country], Triangle_4=mcmc$Triangle.c[4,country],
			k=mcmc$k.c[country], z=mcmc$z.c[country])
	logdens <- 0.0
	res <- .C("dologdensityTrianglekz", x, mean, sd, low, up, par.idx, dlx, 
				mcmc$meta$dl.p1, mcmc$meta$dl.p2,
				lex, idx, length(dct[[country]]), dct[[country]], mcmc$omega*loess.sd[[country]], 
				logdens=logdens)

	return(res$logdens)	
}

z.gibbs.use <- function(mcmc, it) {
	library(msm)
	z <- 0
	sd.d <- d <- 0
	res <- .C("do_z_gibbs", it, mcmc$meta$nr.countries, mcmc$meta$T.end.c, 
					as.numeric(as.matrix(mcmc$meta$e0.matrix)), 
					as.numeric(as.matrix(mcmc$meta$loessSD)), 
					as.numeric(as.matrix(mcmc$Triangle.c)), mcmc$k.c, mcmc$meta$dl.p1, 
					mcmc$meta$dl.p2, mcmc$omega, mcmc$meta$alpha[6], mcmc$meta$delta[6],
						0, 1.15, zmean=d, sd_d=sd.d, z=z)
	#print(res$z)
	z <- rtnorm(1, res$zmean, res$sd_d, 0, 1.15)
	#print(c(res$zmean, res$sd_d, z))
	#stop('')
	#return(res$z)
	return(z)
}

z.gibbs <- function(mcmc, it) {
	library(msm)
	z <- 0
	sd.d <- d <- 0
	test.mc <- get.e0.mcmc('/Users/hana/bayespop/R/LE/wpp2010Fauto')
	mcmc <- test.mc$mcmc.list[[2]]
	#ncountries <- 1
	ncountries <- test.mc$meta$nr.countries
	zmeans <- zs <- zcs <- c()
	s <- summary(test.mc, par.names='omega', par.names.cs=NULL, burnin=20000)
	delta <- k <- c()
	
	for (country in 1:test.mc$meta$nr.countries) {
		country.obj <- get.country.object(country, test.mc$meta, index=TRUE)
		sc <- summary(test.mc, country=country.obj$code, par.names=NULL, 
						par.names.cs=c('Triangle.c', 'k.c', 'z.c'), burnin=20000)
		delta <- cbind(delta, sc$quantiles[1:4, '50%'])
		k <- c(k, sc$quantiles[5, '50%'])
	
#		res <- .C("do_z_gibbs", as.integer(country.obj$code), as.integer(ncountries), 
#					as.integer(mcmc$meta$T.end.c[country]), 
#					as.numeric(as.matrix(mcmc$meta$e0.matrix[,country])), 
#						as.numeric(as.matrix(mcmc$meta$loessSD[,country])), 
#						as.numeric(as.matrix(sc$quantiles[1:4, '50%'])), sc$quantiles[5, '50%'],
#						mcmc$meta$dl.p1, mcmc$meta$dl.p2, s$quantiles['50%'], mcmc$meta$alpha[6],
#						mcmc$meta$delta[6],
#						0, 1.15, zmean=d, sd_d=sd.d)
#		zmeans <- c(zmeans, res$zmean)
#		zs <- c(zs, z)
#		zcs <- c(zcs, sc$quantiles[6, '50%'])
#		cat("\n")
#		print(c(country.obj$name, res$zmean, res$sd_d, z, mcmc$z.c[country]))
	}
	res <- .C("do_z_gibbs", #as.integer(country.obj$code),
					it, 
					as.integer(ncountries), 
					as.integer(mcmc$meta$T.end.c), 
					as.numeric(as.matrix(mcmc$meta$e0.matrix)), 
						as.numeric(as.matrix(mcmc$meta$loessSD)), 
						as.numeric(as.matrix(delta)), k,
						mcmc$meta$dl.p1, mcmc$meta$dl.p2, s$quantiles['50%'], mcmc$meta$alpha[6],
						mcmc$meta$delta[6],
						0, 1.15, zmean=d, sd_d=sd.d, z=z)
		z <- rtnorm(1, res$zmean, res$sd_d, 0, 1.15)

		cat("\n")
		print(c(res$zmean, res$sd_d, z))
	
	#data <- cbind(test.mc$meta$regions$country_code, zmeans, zs, zcs)
#	for (country in 1:test.mc$meta$nr.countries) {
#		country.obj <- get.country.object(country, test.mc$meta, index=TRUE)
#		res <- .C("do_z_gibbs", as.integer(country.obj$code), as.integer(ncountries), as.integer(mcmc$meta$T.end.c[country]), 
#					as.numeric(as.matrix(mcmc$meta$e0.matrix[,country])), 
#						as.numeric(as.matrix(mcmc$meta$loessSD[,country])), 
#						as.numeric(as.matrix(mcmc$Triangle.c[,country])), mcmc$k.c[country],
#						mcmc$meta$dl.p1, mcmc$meta$dl.p2, mcmc$omega, 
#						0, 1.15, zmean=d, sd_d=sd.d)
#		z <- rtnorm(1, res$zmean, res$sd_d, 0, 1.15)
#		zmeans <- c(zmeans, res$zmean)
#		zs <- c(zs, z)
#		cat("\n")
#		print(c(country.obj$name, res$zmean, res$sd_d, z, mcmc$z.c[country]))
#	}
#	data <- cbind(test.mc$meta$regions$country_code, zmeans, zs, mcmc$z.c)
	#colnames(data) <- c('country', 'd', 'z', 'z_c')
	#write.table(data, 'z_c_comparison_median.txt', row.names=FALSE, col.names=TRUE)	
#	res <- .C("do_z_gibbs", it, mcmc$meta$nr.countries, mcmc$meta$T.end.c, as.numeric(as.matrix(mcmc$meta$e0.matrix)), 
#						as.numeric(as.matrix(mcmc$meta$loessSD)), 
#						as.numeric(as.matrix(mcmc$Triangle.c)), mcmc$k.c, mcmc$meta$dl.p1, mcmc$meta$dl.p2, mcmc$omega, 
#						0, 1.15, zmean=d, sd_d=sd.d)
	#print(res$z)
	#z <- rtnorm(1, res$zmean, res$sd_d, 0, 1.15)
	#print(c(res$zmean, res$sd_d, z))
	stop('')
	return(z)
}
