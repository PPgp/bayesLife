\name{e0.raftery.diag}
\alias{e0.raftery.diag}

\title{
Raftery Diagnostics for Parameters of the Life Expectancy
}
\description{
The function computes the Raftery diagnostics for each parameter in the same way as \code{\link[bayesTFR]{tfr.raftery.diag}} of the \pkg{bayesTFR} package.}
\usage{
e0.raftery.diag(mcmc = NULL, sim.dir = file.path(getwd(), "bayesLife.output"), 
    burnin = 0, country = NULL, par.names = e0.parameter.names(), 
    par.names.cs = e0.parameter.names.cs(),
    country.sampling.prop = 1, verbose = TRUE, \dots)
}

\arguments{
  \item{mcmc}{
A \code{bayesLife.mcmc} or \code{bayesLife.mcmc.set} object.
}
  \item{sim.dir}{
Directory with the MCMC simulation results. Only used if \code{mcmc} is \code{NULL}.
}
  \item{burnin}{
Burnin.
}
  \item{country}{Name or code of a country. If it is given, only country-specific parameters parameters of that country are considered.}
  \item{par.names}{
Names of country-independent parameters for which the Raftery diagnostics should be computed.
}
  \item{par.names.cs}{
Names of country-specific parameters for which the Raftery diagnostics should be computed.
}
\item{country.sampling.prop}{Proportion of countries that are included in the diagnostics. It should be between 0 and 1. If it is smaller than 1, the countries are randomly sampled. It is only relevant if \code{par.names.cs} is not \code{NULL}. }
\item{verbose}{Logical switching log messages on and off.}
\item{\dots}{Additional arguments passed to the \code{\link{e0.coda.list.mcmc}} function.}
}
\details{
See \code{\link[bayesTFR]{tfr.raftery.diag}} for details. This function is called from \code{\link{e0.diagnose}}. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Adrian Raftery
}

\seealso{
\code{\link[bayesTFR]{tfr.raftery.diag}}, \code{\link[coda]{raftery.diag}}, \code{\link{e0.diagnose}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
