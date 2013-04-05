\name{get.e0.parameter.traces}
\alias{get.e0.parameter.traces}
\alias{get.e0.parameter.traces.cs}

\title{
Accessing MCMC Parameter Traces
}
\description{
Functions for accessing traces of the MCMC parameters, either country-independent or country-specific.
}
\usage{
get.e0.parameter.traces(mcmc.list, par.names = e0.parameter.names(), 
    burnin = 0, thinning.index = NULL, thin = NULL)
    
get.e0.parameter.traces.cs(mcmc.list, country.obj, 
    par.names = e0.parameter.names.cs(), 
    burnin = 0, thinning.index = NULL, thin = NULL)
}
\arguments{
  \item{mcmc.list}{List of \code{\link{bayesLife.mcmc}} objects.}
  \item{country.obj}{Country object list (see \code{\link[bayesTFR]{get.country.object}}).}
  \item{par.names}{Names of country-independent parameters (in case of \code{get.e0.parameter.traces}) or country-specific parameters (in case of \code{get.e0.parameter.traces.cs}) to be included.}
  \item{burnin}{Burnin indicating how many iterations should be removed from the beginning of each chain.}
  \item{thinning.index}{Index of the traces for thinning. If it is \code{NULL}, \code{thin} is used. \code{thinning.index} does not include \code{burnin} and should be flattened over all chains. For example, if there are two MCMC chains of length 1000, \code{burnin=200} and we want a sample of length 400, then the value should be \code{thinning.index=seq(1,1600, length=400)}.}
  \item{thin}{Alternative to \code{thinning.index}. The above example is equivalent to \code{thin=4}.}
}

\value{
Both functions return a matrix with columns being the parameters and rows being the MCMC values, attached to one another in case of multiple chains. \code{get.e0.parameter.traces} returns country-independent parameters, \code{get.e0.parameter.traces.cs} returns country-specific parameters.
}


\author{
Hana Sevcikova
}


\seealso{
\code{\link{e0.coda.list.mcmc}} for another way of retrieving parameter traces. 
}

\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
m <- get.e0.mcmc(sim.dir)
e0.values <- get.e0.parameter.traces(m$mcmc.list, burnin=10, par.names="z")
hist(e0.values, main=colnames(e0.values))

e0.values.cs <- get.e0.parameter.traces.cs(m$mcmc.list, 
                    get.country.object("Canada", meta=m$meta),
                    burnin=10, par.names="z.c")
hist(e0.values.cs, main=colnames(e0.values.cs))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

