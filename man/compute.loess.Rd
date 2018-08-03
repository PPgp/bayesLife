\name{compute.loess}
\alias{compute.loess}
\alias{compute.residuals}

\title{
Computing Residuals and Fitting Local Polynomial Regression
}
\description{
Functions for computing residuals from the observed life expectancy and MCMC estimation, and fitting a local polynomial regression. 
}
\usage{
compute.residuals(sim.dir, burnin = 1000)

compute.loess(sim.dir = NULL, burnin = 1000, residuals = NULL)
}

\arguments{
  \item{sim.dir}{Directory with the MCMC estimation. In \code{compute.loess}, it is only used if \code{residuals} is not given.}
  \item{burnin}{Number of (unthinned) iterations to be discarded. In \code{compute.loess}, it is only used if \code{residuals} is not given.}
  \item{residuals}{Residuals can be computed outside of the \code{compute.loess} function and passed as the \code{residuals} argument. If not given, the \code{compute.loess} function calls \code{compute.residuals} internally, in which case \code{sim.dir} must be given.}
}

\details{
The Bayesian hierarchical model for life expectancy uses a lowess curve as a multiplier of the variance. The dataset is stored in the package as the \code{loess_sd} dataset. These functions can be used to re-compute this \code{loess_sd} dataset. In such a case, the simulation should be run with the argument \code{constant.variance = TRUE} (in \code{\link{run.e0.mcmc}}).

The residuals are computed for each country as the absolute differences between the observed life expectancy increases and the mean of the estimated double logistic function at the corresponding life expectancy level.
}
\value{
\code{compute.residuals} returns a data frame with columns \sQuote{x} (life expectancy levels) and \sQuote{y} (absolute residuals). 

\code{compute.loess} also returns a data frame with columns \sQuote{x} and \sQuote{y}, where \sQuote{x} is the same as before (with added a minimum and maximum) and \sQuote{y} is the local polynomial fit with constant tails.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{run.e0.mcmc}}
}

\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
resid <- compute.residuals(sim.dir, burnin = 30)
lws <- compute.loess(residuals = resid)

# plot residuals and loess
plot(resid$x, resid$y, ylim = c(0, 4))
lines(lws$x, lws$y, col = "red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ loess }% use one of  RShowDoc("KEYWORDS")

