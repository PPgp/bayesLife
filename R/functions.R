if(getRversion() >= "2.15.1") utils::globalVariables("loess.sd")


g.dl6<-function(x,l, p1, p2){
	dlvalue <- rep(0.0, length(l))
	res <- .C("doDL", x, l, p1, p2, length(l), dl_values=dlvalue)
	return(res$dl_values)
}

loess.lookup<-function(look){
   # call data(loess_sd) before using this function
   idx <- cut(look, loess.sd$x, labels=FALSE, include.lowest = TRUE)
   return(loess.sd$y[idx])
}


dnorm.trunc<-function(x,mean,sd,low,high){
  out<-dnorm(x,mean=mean,sd=sd)/(pnorm(high,mean=mean,sd=sd)-pnorm(low,mean=mean,sd=sd))
  out[x<low]<-0
  out[x>high]<-0
  return(out)
}

rnorm.trunc<-function(mean,sd,low,high){
  temp<--999
  maxit <- 10
  i <- 1
  while((temp<low || temp>high) && i <= maxit) {
     temp<-rnorm(1,mean=mean,sd=sd)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- if(temp<low) low else high
  	warning(paste('Maximum iterations reached in rnorm.trunc(', 
  				mean, ',', sd, '). Value truncated to ', temp, '.', sep=''), immediate.=TRUE)
  }
  return(temp)
}

rgamma.ltrunc<-function(shape,rate,low){
  temp<--999
  maxit <- 10
  i <- 1
  while(temp<low && i <= maxit) {
     temp<-rgamma(1,shape=shape,rate=rate)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- low
  	warning(paste('Maximum iterations reached in rgamma.ltrunc(', shape, ',', rate, ').', sep=''), immediate.=TRUE)
  }
  return(temp)
}

compute.residuals <- function(sim.dir, burnin = 1000) {
    # Compute residuals from an existing simulation.
    # The simulation should have run with constant.variance = TRUE.
    mc <- get.e0.mcmc(sim.dir)
    meta <- mc$meta
    cs.par.names <- c('Triangle.c', 'k.c', 'z.c')
    nT <- dim(meta$e0.matrix)[1]-1
    nC <- get.nrest.countries(meta)
    residuals <- obsdata <- matrix(NA, nrow = nC, ncol = nT)
    dlfunc <- function(t) {
            mean(apply(cs.pars, 1, 
                       function(pars) g.dl6(pars[1:6], l = e0.obs[t], 
                                            p1 = meta$mcmc.options$dl.p1, p2 = meta$mcmc.options$dl.p2)))
    }
    # Compute residuals
    for(cntry in 1:nC) {
        country.obj <- get.country.object(cntry, meta, index = TRUE)
        cs.pars <- get.e0.parameter.traces.cs(mc$mcmc.list, country.obj, 
                                              cs.par.names, burnin = burnin)
        e0.obs <- meta$e0.matrix[, cntry]
        dl <- sapply(1:nT, dlfunc)
        residuals[cntry, ] <- abs(e0.obs[-length(e0.obs)] + dl - e0.obs[-1])
        obsdata[cntry, ] <- e0.obs[-1]
    }
    return(na.omit(data.frame(x = as.vector(obsdata), y = as.vector(residuals))))
}

do.compute.loess <- function(df) {
    lfit <- lowess(df)
    lfun <- approxfun(lfit)
    x <- sort(unique(df$x))
    lws <- list(x=x, y=lfun(x))
    # add additional point at the end to be able to bin it
    lws$x <- c(lws$x, 999)
    lws$y <- c(lws$y, lws$y[length(lws$y)])
    if(x[1] > 15) { # set the first point to 15 so that extreme cases also work (Cambodia has 18.12)
        lws$x <- c(15, lws$x)
        lws$y <- c(lws$y[1], lws$y)
    }
    lws
}

compute.loess <- function(sim.dir = NULL, burnin = 1000, residuals = NULL) {
    # Compute the loess curves from either 
    # an existing simulation (sim.dir) or from given residuals.
    # If sim.dir is used, the simulation should have run with constant.variance=TRUE.
    # If residuals are given, it should have columns x and y.
    
    if(is.null(residuals))
        residuals <- compute.residuals(sim.dir, burnin = burnin)
    
    return(do.compute.loess(residuals))
}

