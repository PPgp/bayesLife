if(getRversion() >= "2.15.1") utils::globalVariables("loess.sd")


g.dl6<-function(x,l, p1, p2){
	dlvalue <- rep(0.0, length(l))
	res <- .C("doDL", x, l, p1, p2, length(l), dl_values=dlvalue)
	return(res$dl_values)
}

loess.lookup <- function(look) {
   # call data(loess_sd) before using this function
   idx <- cut(look, loess.sd$x, labels=FALSE, include.lowest = TRUE)
   loess.sd$y[idx]
}

loess.lookup.hiv.model <- function(look, is.hiv) {
    find.look <- function(value, hiv) {
        look.in <- if(hiv) loess.sd$hiv else loess.sd
        idx <- cut(value, look.in$x, labels=FALSE, include.lowest = TRUE)
        look.in$y[idx]
    }
    mapply(find.look, look, is.hiv)
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


