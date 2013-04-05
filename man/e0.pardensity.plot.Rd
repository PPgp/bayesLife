\name{e0.pardensity.plot}
\alias{e0.pardensity.plot}
\alias{e0.pardensity.cs.plot}

\title{
Plotting MCMC Parameter Density
}
\description{
Functions for plotting density of the posterior distribution of the MCMC parameters.
}
\usage{
e0.pardensity.plot(mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    chain.ids = NULL, par.names = e0.parameter.names(), 
    burnin = NULL, dev.ncol = 5, low.memory = TRUE, ...)
    
e0.pardensity.cs.plot(country, mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    chain.ids = NULL, par.names = e0.parameter.names.cs(), 
    burnin = NULL, dev.ncol = 3, low.memory = TRUE, ...)

}
\arguments{
  \item{country}{Name or numerical code of a country.}
  \item{mcmc.list}{List of \code{\link{bayesLife.mcmc}} objects, or an object of class \code{\link{bayesLife.mcmc.set}} or of class \code{\link{bayesLife.prediction}}. If it is \code{NULL}, the parameter values are loaded from \code{sim.dir}.}
  \item{sim.dir}{Directory with the MCMC simulation results. It is only used if \code{mcmc.list} is \code{NULL}.}
  \item{chain.ids}{List of MCMC identifiers to be plotted. If it is \code{NULL}, all chains found in \code{mcmc.list} or \code{sim.dir} are plotted.}
  \item{par.names}{Names of parameters for which density should be plotted. By default all country-independent parameters are plotted if used within \code{e0.pardensity.plot}, or all country-specific parameters are plotted if used within \code{e0.pardensity.cs.plot}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of each chain.}
  \item{dev.ncol}{Number of columns for the graphics device. If the number of parameters is smaller than \code{dev.ncol}, the number of columns is automatically decreased.}
  \item{low.memory}{Logical indicating if the processing should run in a memory-efficient mode.}
  \item{\dots}{Further arguments passed to the \code{density} function.}
}
\details{
The functions plot the density of the posterior distribution either for country-independent parameters (\code{e0.pardensity.plot}) or for country-specific parameters (\code{e0.pardensity.cs.plot}), one graph per parameter.  One can restrict it to specific chains by setting the \code{chain.ids} argument and to specific parameters by setting the \code{par.names} argument. 

If \code{mcmc.list} is an object of class \code{\link{bayesLife.prediction}}, thinned traces are used instead of the full chains. In such a case, \code{burnin} and \code{chain.ids} cannot be modified - their value is set to the one used when the thinned traces were created, namely when running \code{\link{e0.predict}}.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{e0.partraces.plot}}
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
e0.pardensity.plot(sim.dir=sim.dir, burnin=10)
e0.pardensity.cs.plot(country="Ireland", sim.dir=sim.dir, burnin=10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

