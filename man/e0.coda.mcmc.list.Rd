\name{e0.coda.list.mcmc}
\alias{e0.coda.list.mcmc}
\alias{coda.mcmc.bayesLife.mcmc}

\title{
Convertion to coda's Objects
}
\description{
The functions convert MCMC traces (simulated using \code{\link{run.e0.mcmc}}) into objects that can be used with the \pkg{coda} package.
}
\usage{
e0.coda.list.mcmc(mcmc.list = NULL, country = NULL, chain.ids = NULL, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    par.names = e0.parameter.names(), 
    par.names.cs = e0.parameter.names.cs(), low.memory = FALSE, \dots)
    
\method{coda.mcmc}{bayesLife.mcmc}(mcmc, country = NULL, par.names = e0.parameter.names(), 
    par.names.cs = e0.parameter.names.cs(), \dots)
}

\arguments{
  \item{mcmc.list}{List of \code{\link{bayesLife.mcmc}} objects, or an object of class \code{\link{bayesLife.mcmc.set}} or \code{\link{bayesLife.prediction}}. If it is \code{NULL}, the MCMCs are loaded from \code{sim.dir}. Either \code{mcmc.list} or \code{sim.dir} must be given.}
  \item{mcmc}{Object of class \code{\link{bayesLife.mcmc}}.}
  \item{country}{Country name or code. It is used in connection with the \code{par.names.cs} argument (see below).}
  \item{chain.ids}{Vector of chain identifiers. By default, all chains available in the \code{mcmc.list} object are included.
}
  \item{sim.dir}{Directory with the MCMC simulation results. Only used if \code{mcmc.list} is \code{NULL}.}
  \item{par.names}{Names of country-independent parameters to be included.}
  \item{par.names.cs}{Names of country-specific parameters to be included. The argument \code{country} is used to filter out traces that correspond to a specific country. If \code{country} is not given, for each parameter, traces for all countries are included.}
  \item{low.memory}{Logical indicating if the function should run in a memory-efficient mode.}
  \item{\dots}{Additional arguments passed to the \pkg{coda}'s \code{\link[coda]{mcmc}} function, such as \code{burnin} and \code{thin}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The function \code{e0.coda.list.mcmc} returns an object of class \dQuote{mcmc.list}. The function \code{coda.mcmc} returns an object of class \dQuote{mcmc}, both defined in the \pkg{coda} package.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{e0.partraces.plot}} for plotting the MCMC traces and \code{\link{summary.bayesLife.mcmc.set}}.
}
\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
coda.list <- e0.coda.list.mcmc(sim.dir=sim.dir, country="France", burnin=30)
summary(coda.list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
