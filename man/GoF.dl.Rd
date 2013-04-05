\name{e0.dl.coverage}
\alias{e0.dl.coverage}
\title{
Goodness of Fit of the Double Logistic Function
}
\description{
The function computes coverage, i.e. the ratio of observed data fitted within the given probability intervals of the predictive posterior distribution of the double logistic function, as well as the root mean square error of the simulation.
}
\usage{
e0.dl.coverage(sim.dir, pi = c(80, 90, 95), burnin = 10000, verbose = TRUE)
}

\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results. If a prediction and its corresponding thinned mcmcs are available in the simulation directory, those are taken for assessing the goodness of fit.}
  \item{pi}{Probability interval. It can be a single number or an array.}
    \item{burnin}{Burnin. Only relevant if \code{sim.dir} does not contain thinned chains.}
  \item{verbose}{Logical switching log messages on and off.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{List with the same components as \code{\link[bayesTFR]{tfr.dl.coverage}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
\note{
To see the fit visually per country, use \code{\link{e0.DLcurve.plot}(\dots, predictive.distr=TRUE,\dots)}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e0.DLcurve.plot}}
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
e0 <- get.e0.mcmc(sim.dir)
# Note that this simulation is a toy example and thus has not converged.
gof <- e0.dl.coverage(sim.dir)
gof$country.coverage
e0.DLcurve.plot(e0, country=608, predictive.distr=TRUE, pi=c(80, 90, 95))
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

