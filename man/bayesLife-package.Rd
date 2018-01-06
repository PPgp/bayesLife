\name{bayesLife-package}
\alias{bayesLife-package}
\alias{bayesLife}
\docType{package}
\title{
Bayesian Projection of the Life Expectancy}
\description{
Collection of functions for making probabilistic projections of the life expectancy for all countries of the world, 
using a Bayesian hierarchical model and the United Nations demographic time series.

}
\details{
\tabular{ll}{
Package: \tab bayesLife\cr
Type: \tab Package\cr
Version: \tab 3.2-0\cr
Date: \tab 2017-10-13\cr
License: \tab GPL (>= 2)\cr
URL: \tab \url{https://bayespop.csss.washington.edu}\cr
}
The projection follows a method developed by Chunn et al (2010).  It uses historical data provided by the United Nations to simulate a posterior distribution of the life expectancy for all countries in the world simultaneously. 

The package is implemented in a similar way as the \pkg{bayesTFR} package and thus, many functions have their equivalents in \pkg{bayesTFR}. The main functions of the \pkg{bayesLife} package are:
\itemize{
\item \link{run.e0.mcmc}: Runs a Markov Chain Monte Carlo (MCMC) simulation for one or more chains, possibly in parallel.  It results in a posterior sample of the MCMC parameters.  Existing simulation runs can be resumed using \link{continue.e0.mcmc}.
\item \link{e0.predict}: Using the posterior parameter samples it derives posterior trajectories of the life expectancy for all countries.
\item \link{e0.jmale.predict}: Given existing outputs of \link{e0.predict} for female life expectancy, this function estimates and predicts a joint male life expectancy as described in Raftery et al (2014).
}
A number of functions analyzing results are included in the package:
\itemize{
\item \link{e0.trajectories.plot}: Shows the posterior trajectories for a given country, including their median and given probability intervals.
\item \link{e0.trajectories.table}: Shows the posterior trajectories for a given country in a tabular form.
\item \link{e0.map}: Shows a world map of life expectancy for a given projection period.
\item \link{e0.DLcurve.plot}: Shows the posterior curves of the double logistic function used in the simulation, including their median and given probability intervals.
\item \link{e0.partraces.plot} and \link{e0.partraces.cs.plot}: Plot the MCMC traces of country-independent parameters and country-specific parameters, respectively.
\item \link{e0.pardensity.plot} and \link{e0.pardensity.cs.plot}: Plot the posterior density of the MCMCs for country-independent parameters and country-specific parameters, respectively.
\item \link{summary.bayesLife.mcmc.set}: Summary function for the MCMC results.
\item \link{summary.bayesLife.prediction}: Summary function for the prediction results.
}
For MCMC diagnostics, function \link{e0.coda.list.mcmc} creates an object of type \dQuote{mcmc.list} that can be used with the \pkg{coda} package.  Furthermore, function \link{e0.diagnose} analyzes the MCMCs using the Raftery diagnostics implemented in the \pkg{coda} package and gives information about parameters that did not converge. Function \link{e0.dl.coverage} computes a goodness of fit of the double logistic function.

Existing simulation results can be accessed using the \link{get.e0.mcmc} function.  An existing prediction can be accessed via \link{get.e0.prediction}.

Historical data are taken from one of the packages \pkg{wpp2017} (default), \pkg{wpp2015}, \pkg{wpp2012} or \pkg{wpp2010}, depending on users settings.

}

\note{
There is a directory \code{ex-data} shipped with the package which contains results from an example simulation, containing one chain with 60 iterations.  The Example section below shows how these results were created.  These data are used in Example sections throughout the manual.  The user can either reproduce the data in her/his local directory, or use the ones from the package.

%A \sQuote{real} (converged) simulation including predictions can be downloaded from \url{http://bayespop.csss.washington.edu}.   
}

\author{
Hana Sevcikova, Adrian Raftery, Jennifer Chunn

Maintainer: Hana Sevcikova <hanas@uw.edu>
}
\references{
J. L. Chunn, A. E. Raftery, P. Gerland (2010): Bayesian Probabilistic Projections of Life Expectancy for All Countries.  Working Paper nr. 105, Center for Statistics and the Social Sciences, University of Washington. \url{http://www.csss.washington.edu/Papers}

A. E. Raftery,  N. Li, H. Sevcikova , P. Gerland, G. K. Heilig (2012). Bayesian probabilistic population projections for all countries. Proceedings of the National Academy of Sciences 109:13915-13921.

A. E. Raftery,  N. Lalic, P. Gerland (2014). Joint Probabilistic Projection of Female and Male Life Expectancy. Demographic Research, 30:795-822.

%Lalic, N. (2011). Master's thesis at the Department of Statistics, University of Washington.
}

\keyword{ package }
\seealso{
  \code{\link[bayesTFR:bayesTFR-package]{bayesTFR}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.e0.mcmc(sex='F', nr.chains=1, iter=60, seed=1, thin=1, 
     output.dir=sim.dir, verbose=TRUE)
pred <- e0.predict(m, burnin=30, verbose=TRUE)
summary(pred, country="Canada")
unlink(sim.dir, recursive=TRUE)}
}
%# exact commands to create data in ex-data directory
% m <- run.e0.mcmc(sex='F', nr.chains=1, iter=60, seed=1, thin=1, verbose=TRUE, compression.type='bz', start.year=1950)
% pred <- e0.predict(m, burnin=30, verbose=TRUE, predict.jmale=FALSE, save.as.ascii=0) 
%# to save disk space, remove thinned_* directory and csv files from predictions dir