\name{e0.DLcurve.plot}
\alias{e0.DLcurve.plot}
\alias{e0.DLcurve.plot.all}

\title{
Plotting  Posterior Distribution of the Double Logistic Function of Life Expectancy
}
\description{
The functions plot the posterior distribution of the double logistic function used in the simulation, including their median and given probability intervals.
}
\usage{
e0.DLcurve.plot(mcmc.list, country, burnin = NULL, pi = 80, 
    e0.lim = c(20, 90), nr.curves = NULL, ylim = NULL, 
    xlab = "e(0)", ylab = "5-year gains", main = NULL, \dots)
    
e0.DLcurve.plot.all(mcmc.list=NULL, sim.dir=NULL, 
    output.dir=file.path(getwd(), 'DLcurves'),
    output.type="png", burnin=NULL, verbose=FALSE, \dots)
}
\arguments{
  \item{mcmc.list}{List of \code{\link{bayesLife.mcmc}} objects, an object of class \code{\link{bayesLife.mcmc.set}} or of class \code{\link{bayesLife.prediction}}. In case of \code{e0.DLcurve.plot.all} if it si \code{NULL}, it is loaded from \code{sim.dir}.}
  \item{country}{Name or numerical code of a country.}
  \item{burnin}{Number of iterations to be discarded from the beginning of parameter traces.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{e0.lim}{Tuple of the minimum and maximum life expectancy to be shown in the plot.}
  \item{nr.curves}{Number of curves to be plotted. If \code{NULL}, all curves are plotted.}
  \item{ylim, xlab, ylab, main}{Graphical parameters passed to the \code{plot} function.}
  \item{\dots}{Additional graphical parameters. In addition, for \code{e0.DLcurve.plot.all}, 
  	\code{\dots} contains arguments \code{pi}, \code{e0.lim} and \code{nr.curves}.}
  \item{sim.dir}{Directory with the simulation results. Only relevant, if \code{mcmc.list} is \code{NULL}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}

}
\details{\code{e0.DLcurve.plot} plots double logistic curves for the given country.  \code{e0.DLcurve.plot.all} creates such plots for all countries and stores them in \code{output.dir}. Parameters passed to the double logistic function are either thinned traces created by the \code{\link{e0.predict}} function (if \code{mcmc.list} is an object of class \code{\link{bayesLife.prediction}}), or they are selected by equal spacing from the MCMC traces. In the former case, \code{burnin} is set automatically; in the latter case, \code{burnin} defaults to 0. If \code{nr.curves} is smaller than 2000, the median and probability intervals are computed on a sample of 2000 equally spaced data points, otherwise on all plotted curves.}

\author{
Hana Sevcikova
}

\examples{
\dontrun{
sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
mcmc.set <- get.e0.mcmc(sim.dir=sim.dir)
e0.DLcurve.plot(country='Japan', mcmc.set, burnin=40)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

