\name{e0.GoF.dl}
\alias{e0.GoF.dl}
\title{
Goodness of Fit of the Double Logistic Function
}
\description{
The function computes the ratio of observed data fitted within the given probability intervals of the predictive posterior distribution of the double logistic function. 
}
\usage{
e0.GoF.dl(sim.dir, pi = c(80, 90, 95), burnin = 20000, verbose = TRUE)
}

\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results. If a prediction and its corresponding thinned mcmcs are available in the simulation directory, those are taken for assessing the goodness of fit.}
  \item{pi}{Probability interval. It can be a single number or an array.}
    \item{burnin}{Burnin. Only relevant if \code{sim.dir} does not contained thinned chains.}
  \item{verbose}{Logical switching log messages on and off.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{List with the following components:
  \item{total}{Vector of the total goodness of fit ratios, one element per probability interval.}
  \item{time}{Matrix corresponding to the goodness of fit computed per time period. (Rows correspond to probability intervals, columns correspond to time.)}
  \item{country}{Matrix corresponding to the goodness of fit computed per country. (Rows correspond to probability intervals, columns correspond to countries.)}
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
sim.dir <- file.path(.find.package("bayesLife"), "ex-data", "bayesLife.output")
e0 <- get.e0.mcmc(sim.dir)
# Note that this simulation is a toy example and thus has not converged.
e0.GoF.dl(sim.dir)
e0.DLcurve.plot(e0, country=608, predictive.distr=TRUE, pi=c(80, 90, 95))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

