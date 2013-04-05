\name{e0.partraces.plot}
\alias{e0.partraces.plot}
\alias{e0.partraces.cs.plot}

\title{
Plotting MCMC Parameter Traces
}
\description{
Functions for plotting the MCMC parameter traces.
}
\usage{
e0.partraces.plot(mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    chain.ids = NULL, par.names = e0.parameter.names(), 
    nr.points = NULL, dev.ncol = 5, low.memory = TRUE, \dots)
    
e0.partraces.cs.plot(country, mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    chain.ids = NULL, par.names = e0.parameter.names.cs(), 
    nr.points = NULL, dev.ncol = 3, low.memory = TRUE, \dots)
}
\arguments{
  \item{country}{Name or numerical code of a country.}
  \item{mcmc.list}{List of \code{bayesLife.mcmc} objects, or an object of class \code{bayesLife.mcmc.set} or of class \code{bayesLife.prediction}. If it is \code{NULL}, the traces are loaded from \code{sim.dir}.}
  \item{sim.dir}{Directory with the MCMC simulation results. It is only used if \code{mcmc.list} is \code{NULL}.}
  \item{chain.ids}{List of MCMC identifiers to be plotted. If it is \code{NULL}, all chains found in \code{mcmc.list} or \code{sim.dir} are plotted.}
  \item{par.names}{Names of parameters for which traces should be plotted. By default all country-independent parameters are plotted if used within \code{e0.partraces.plot}, or country-specific parameters are plotted if used within \code{e0.partraces.cs.plot}.}
  \item{nr.points}{Number of points to be plotted. If \code{NULL}, all (stored) points are plotted, otherwise the traces are thinned evenly.}
  \item{dev.ncol}{Number of column for the graphics device. If the number of parameters is smaller than \code{dev.ncol}, the number of columns is automatically decreased.}
  \item{low.memory}{Logical indicating if the processing should run in a low-memory mode. If it is \code{FALSE}, traces of all available parameters are loaded into memory. Otherwise, parameters are loaded as they are needed and are not kept in the memory.}
  \item{\dots}{Additional graphical arguments. It can also contain the arguments \code{burnin} or \code{thin}.}
}
\details{
The functions plot MCMC traces either for country-independent parameters (\code{e0.partraces.plot}) or for country-specific parameters (\code{e0.partraces.cs.plot}), one graph per parameter.  One can restrict it to specific chains by setting the \code{chain.ids} argument, and to specific parameters by setting the \code{par.names} argument.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{e0.coda.list.mcmc}} and \code{\link{get.e0.parameter.traces}} for retrieving the raw values of traces.
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
e0.partraces.plot(sim.dir=sim.dir)
e0.partraces.cs.plot(country="Ireland", sim.dir=sim.dir)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
