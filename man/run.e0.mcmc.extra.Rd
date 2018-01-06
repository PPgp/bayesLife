\name{run.e0.mcmc.extra}
\alias{run.e0.mcmc.extra}

\title{
Run MCMC for Extra Countries, Areas or Regions
}
\description{
Run MCMC for extra countries, areas or regions. It uses the posterior distribution of model hyperparameters from an existing simulation to generate country-specific parameters. 
}
\usage{
run.e0.mcmc.extra(sim.dir = file.path(getwd(), "bayesLife.output"), 
    countries = NULL, my.e0.file = NULL, 
    iter = NULL, thin = 1, burnin = 2000, country.overwrites = NULL, 
    parallel = FALSE, nr.nodes = NULL,  my.locations.file = NULL,
    verbose = FALSE, verbose.iter = 100, ...)
}
\arguments{
  \item{sim.dir}{
Directory with an existing simulation.
}
  \item{countries}{Vector of country codes. These include codes of areas and regions (see column \code{country_code} in \code{\link[wpp2017]{UNlocations}}).
}
  \item{my.e0.file}{
File name containing user-specified time series of life expectancy for countries for which the simulation should run (see Details below).
}
  \item{iter}{
Number of iterations to be used for sampling from the posterior distribution of the hyperparameters. By default, the number of (possibly thinned) iterations used in the existing simulation is taken.
}
  \item{thin}{
Thinning interval for sampling from the posterior distribution of the hyperparameters.
}
  \item{burnin}{
Number of iterations discarded before sampling from the posterior distribution of the hyperparameters.
}
\item{country.overwrites}{This argument allows to overwrite some of the prior parameters for specific countries. If it is not \code{NULL} it should be a data frame with an obligatory column \sQuote{country_code}. Each row then corresponds to one country. Other columns can be \sQuote{k.c.prior.low}, \sQuote{k.c.prior.up}, \sQuote{z.c.prior.low}, \sQuote{z.c.prior.up}, \sQuote{Triangle_\eqn{x}.c.prior.low} and \sQuote{Triangle_\eqn{x}.c.prior.up} where \eqn{x} can be an integer from 1 to 4. Rows corresponding to countries that are not processed in this function are ignored.}
  \item{parallel}{
Logical determining if the simulation should run multiple chains in parallel.
}
  \item{nr.nodes}{
Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. By default it equals to the number of chains contained in the existing simulation.
}
\item{my.locations.file}{File name containing user-specified locations. See Details below.}
  \item{verbose}{
Logical switching log messages on and off.
}
\item{verbose.iter}{Integer determining how often (in number of iterations) log messages are outputted during the estimation.}
\item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
}
\details{
The function can be used to make predictions for countries, areas or regions (further denoted as \sQuote{countries}) that were not included in the MCMC estimation (invoked by \code{\link{run.e0.mcmc}}). It creates MCMC traces for country-specific parameters. The purpose of this function is to have country-specific parameters available in order to be able to generate projections for additional countries or their aggregations, without having to re-run the often time-expensive MCMC simulation.

The set of countries to be considered by this function can be given either by their codes, using the argument \code{countries}, in which case the countries must be included in the UN WPP \code{\link[wpp2017]{e0}} dataset. Or, it can be given by a user-specific file, using the argument \code{my.e0.file}. The function considers a union of both arguments. The function will ignore all countries that were used in the existing MCMC simulation for estimating the hyperparameters. Countries that already own country-specific parameters (e.g. because they were included in \code{my.e0.file} passed to \code{\link{run.e0.mcmc}}) get their parameters recomputed. Note that all countries must be included in the \code{\link[wpp2017]{UNlocations}} dataset, but unlike in \code{\link{run.e0.mcmc}}, their  \code{include_code} is ignored. As in the case of \code{\link{run.e0.mcmc}}, the default dataset of locations \code{\link[wpp2017]{UNlocations}} can be overwritten using a file of the same structure as \code{\link[wpp2017]{UNlocations}} passed via the \code{my.locations.file} argument. This file should be especially used, if e0 is simulated for new locations that are not included in \code{\link[wpp2017]{UNlocations}}.

}
\note{If there is an existing projection for the directory \code{sim.dir}, use \code{\link{e0.predict.extra}} to obtain projections for the extra countries used in this function.}
\value{
An object of class \code{\link{bayesLife.mcmc.set}}.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{run.e0.mcmc}}, \code{\link{e0.predict.extra}}
}
\examples{
\dontrun{
m <- run.e0.mcmc(nr.chains=1, iter=20, thin=1, verbose=TRUE)
m <- run.e0.mcmc.extra(countries=c(908,924), burnin=10, verbose=TRUE)
summary(m, country=924)
pred <- e0.predict(burnin=10, verbose=TRUE)
summary(pred, country=908)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}

