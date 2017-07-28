\name{summary.bayesLife.mcmc.set}
\Rdversion{1.1}
\alias{summary.bayesLife.mcmc.set}
\alias{print.summary.bayesLife.mcmc.set}
\alias{print.summary.bayesLife.mcmc.set.meta}
\alias{summary.bayesLife.mcmc}
\title{
Summary Statistics for Life Expectancy MCMCs
}
\description{
Summary of an object \code{\link{bayesLife.mcmc.set}} or \code{\link{bayesLife.mcmc}}, computed via \code{\link{run.e0.mcmc}}.  It can be obtained either for all countries or for a specific country, and either for all parameters or for specific parameters.  The function uses the \code{\link[coda]{summary.mcmc}} function of the \pkg{coda} package.
}
\usage{
\method{summary}{bayesLife.mcmc.set}(object, country = NULL, chain.id = NULL, 
    par.names = e0.parameter.names(), 
    par.names.cs = e0.parameter.names.cs(), 
    meta.only = FALSE, thin = 1, burnin = 0, \dots)
	
\method{summary}{bayesLife.mcmc}(object, country = NULL,
    par.names = e0.parameter.names(), 
    par.names.cs = e0.parameter.names.cs(), \dots)
}

\arguments{
  \item{object}{Object of class \code{\link{bayesLife.mcmc.set}} or \code{\link{bayesLife.mcmc}}.}
  \item{country}{Country name or code if a country-specific summary is desired.}
  \item{chain.id}{Identifiers of MCMC chains. By default, all chains are considered.}
  \item{par.names}{Country independent parameters to be included in the summary.}
  \item{par.names.cs}{Country-specific parameters to be included in the summary.}
  \item{meta.only}{Logical. If it is \code{TRUE}, only meta information of the simulation is included.}
  \item{thin}{Thinning interval. Only used if larger than the \code{thin} argument used in \code{\link{run.e0.mcmc}}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of each chain before computing the summary.}
  \item{\dots}{For \code{summary.bayesLife.mcmc}, arguments \code{thin} and \code{burnin} can be passed here. In addition, both functions accept arguments passed to the \code{\link[coda]{summary.mcmc}} function of the \pkg{coda} package.}
}

\author{
Hana Sevcikova
}


\seealso{
\code{\link{bayesLife.mcmc.set}}, \code{\link[coda]{summary.mcmc}}
}
\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
m <- get.e0.mcmc(sim.dir)
summary(m, country="Czechia", burnin=20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }

