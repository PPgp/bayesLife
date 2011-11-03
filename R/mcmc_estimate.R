e0.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE) {
	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	niter <- mcmc$iter
	meta <- mcmc$meta
	lex.mtx <- meta$e0.matrix
	T <- nrow(lex.mtx) - 1
	C <- meta$nr.countries
	mcmc$thin <- thin
	
	delta.sq <- mcmc$meta$delta^2
	dlf <- matrix(NA, nrow=T, ncol=C)
	psi.shape <- (sum(meta$T.end.c-1)-1)/2.
	
	if (start.iter > niter) return(mcmc)
	
	mcenv <- new.env() # Create an environment for the mcmc stuff in order to avoid 
						# copying of the mcmc list 
    for (item in names(mcmc)) mcenv[[item]] <- mcmc[[item]]
    
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
		for (i in 1:4) {
			# Triangle is truncated normal in [0,100]
			mcenv$Triangle[i] <- rnorm.trunc(mean=Tr.mean[i], sd=Tr.sd[i], low=0, high=100)
		}
		# k is truncated normal in [0,10]
		mcenv$k <- rnorm.trunc(mean=Tr.mean[5], sd=Tr.sd[5], low=0, high=10)
		# z is truncated normal in [0,1.15]
		mcenv$z <- rnorm.trunc(mean=Tr.mean[6], sd=Tr.sd[6], low=0, high=1.15)
		
		#ratesum <- 0		
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in 1:C) {
			Triangle.k.z.c.update(mcenv, country)
			idx <- 1:(meta$T.end.c[country]-1)
			dlf[idx,country] <- g.dl6(c(mcenv$Triangle.c[,country], 
								mcenv$k.c[country], mcenv$z.c[country]), 
								lex.mtx[idx+1, country], meta$dl.p1, meta$dl.p2)
			#ratesum <- ratesum + ((meta$d.ct[idx,country]-dlf[idx,country])^2)/(meta$loessSD[idx,country]^2)
		}
		# Update omega using Gibbs sampler
		###########################################
		#psi <- rgamma.ltrunc(shape=psi.shape, rate=0.5*ratesum, low=0.01)
		#mcmc$omega <- 1/sqrt(psi)
		update.omega(mcenv, dlf)
		
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
	return(as.list(mcenv))
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
	updated.var.names <- c('Triangle.c', 'k.c', 'z.c')

	for(iter in 1:niter) {
		if(verbose || (iter %% 10 == 0))
			cat('\nIteration:', iter, '--', date())
		# set hyperparameters for this iteration
        for (par in hyperparameter.names) {
        	if(is.null(dim(hyperparameters[[par]]))) {
        		mcmc[[par]] <- hyperparameters[[par]][iter]
        	} else {
        		mcmc[[par]] <- hyperparameters[[par]][iter,]
        	}
        }
						
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in countries) {
			mcmc <- Triangle.k.z.c.update(mcmc, country)
		}

		################################################################### 
        # write samples simu/thin to disk
        ##################################################################
         #update the original mcmc with the new values
         for(var in updated.var.names) {
         	mcmc.orig[[var]] <- mcmc[[var]]
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

Triangle.k.z.c.update <- function(mcmc, country) {
	# Update Triangle pars using slice sampling
	sigmas <- 1/sqrt(mcmc$lambda)
	Triangle.c.width <- mcmc$meta$Triangle.c.width
	k.c.width <- mcmc$meta$k.c.width
	z.c.width <- mcmc$meta$z.c.width
	idx <- 1:(mcmc$meta$T.end.c[country]-1)
	lex <- mcmc$meta$e0.matrix[2:mcmc$meta$T.end.c[country], country]
	for (i in 1:4) {
		mcmc$Triangle.c[i, country] <- slice.sampling(mcmc$Triangle.c[i, country],
										logdensity.Triangle.k.z.c, Triangle.c.width[i], 
										mean=mcmc$Triangle[i], 
										sd=sigmas[i], low=0, up=100, par.idx=i, 
										mcmc=mcmc, 
										country=country, lex=lex, idx=idx)
	}
	mcmc$k.c[country] <- slice.sampling(mcmc$k.c[country],
										logdensity.Triangle.k.z.c, k.c.width, mean=mcmc$k, 
										sd=1/sqrt(mcmc$lambda.k),
										low=0, up=10, par.idx=5, mcmc=mcmc, 
										country=country, lex=lex, idx=idx)
	mcmc$z.c[country] <- slice.sampling(mcmc$z.c[country],
										logdensity.Triangle.k.z.c, z.c.width, mean=mcmc$z, 
										sd=1/sqrt(mcmc$lambda.z),
										low=0, up=1.15, par.idx=6, mcmc=mcmc, 
										country=country, lex=lex, idx=idx)
	return()
}

lambdas.update <- function(mcmc) {
	# Update lambdas using MH-algorithm
	tau.sq <- mcmc$meta$tau^2
	for(i in 1:4) {
		prop <- proposal.lambda(mcmc$meta$nu, tau.sq[i], mcmc$Triangle[i], mcmc$Triangle.c[i,], 
								mcmc$meta$nr.countries)
		lpx0 <- logdensity.lambda(mcmc$lambda[i], mcmc, 100, mcmc$Triangle[i], mcmc$Triangle.c[i,], tau.sq[i])
		lpx1 <- logdensity.lambda(prop, mcmc, 100, mcmc$Triangle[i], mcmc$Triangle.c[i,], tau.sq[i])
		prob_accept <- min(exp(lpx1 - lpx0), 1)
		if (runif(1) < prob_accept){
			mcmc$lambda[i] <- prop	
		}
	}
	# update lambda.k 
	prop <- proposal.lambda(mcmc$meta$nu, tau.sq[5], mcmc$k, mcmc$k.c, mcmc$meta$nr.countries)
	lpx0 <- logdensity.lambda(mcmc$lambda.k, mcmc, 10, mcmc$k, mcmc$k.c, tau.sq[5])
	lpx1 <- logdensity.lambda(prop, mcmc, 10, mcmc$k, mcmc$k.c, tau.sq[5])
	prob_accept <- min(exp(lpx1 - lpx0), 1)
	if (runif(1) < prob_accept){
		mcmc$lambda.k <- prop	
	}
	# update lambda.z
	prop <- proposal.lambda(mcmc$meta$nu, tau.sq[6], mcmc$z, mcmc$z.c, mcmc$meta$nr.countries)
	lpx0 <- logdensity.lambda(mcmc$lambda.z, mcmc, 1.15, mcmc$z, mcmc$z.c, tau.sq[6])
	lpx1 <- logdensity.lambda(prop, mcmc, 1.15, mcmc$z, mcmc$z.c, tau.sq[6])
	prob_accept <- min(exp(lpx1 - lpx0), 1)
	if (runif(1) < prob_accept){
		mcmc$lambda.z <- prop	
	}
	return()
}

logdensity.lambda <- function(lambda, mcmc, high, Triangle, Triangle.c, tau.sq, ...) {
	nu <- mcmc$meta$nu
	return(log(dgamma(lambda, nu/2, rate=tau.sq)) + 
				sum(log(dnorm.trunc(Triangle.c, mean=Triangle, sd=1/sqrt(lambda), 0, high))))
}

proposal.lambda <- function(nu, tau.sq, Triangle, Triangle.c, C) {
	return(rgamma(1, (nu+C)/2, rate=tau.sq+0.5*sum((Triangle.c-Triangle)^2)))
}

update.omega <- function(mcmc, dlf) {
	omega.width <- 1		
	#print('Slice sampling for omega')
	mcmc$omega <- slice.sampling(mcmc$omega, log.density.omega, omega.width, 
								mcmc=mcmc, dlf=dlf, low=0, up=10)
	return()
}

log.density.omega <- function(x, mcmc, dlf, low, up) {
	log.d.cDens <- 0
	for(country in 1:mcmc$meta$nr.countries) 
		log.d.cDens <- log.d.cDens + sum(log(pmax(dnorm(mcmc$meta$d.ct[,country], dlf[,country], 
										sd=x*mcmc$meta$loessSD[,country]),1e-100)))
	return(log(dunif(x, low, up)) + log.d.cDens)
}


logdensity.Triangle.k.z.c <- function(x, mean, sd, low, up, par.idx, mcmc, country, lex, idx) {
	dlx <- c(Triangle_1=mcmc$Triangle.c[1,country], Triangle_2=mcmc$Triangle.c[2,country],
			Triangle_3=mcmc$Triangle.c[3,country], Triangle_4=mcmc$Triangle.c[4,country],
			k=mcmc$k.c[country], z=mcmc$z.c[country])
	logdens <- 0.0
	res <- .C("dologdensityTrianglekz", x, mean, sd, low, up, par.idx, dlx, 
				mcmc$meta$dl.p1, mcmc$meta$dl.p2,
				lex, idx, length(idx), mcmc$meta$d.ct[idx,country], mcmc$omega*mcmc$meta$loessSD[idx,country], 
				logdens=logdens)

	return(res$logdens)	
}

