\name{e0.DLcurve.plot}
\alias{e0.DLcurve.plot}
\alias{e0.DLcurve.plot.all}
\alias{e0.world.dlcurves}
\alias{e0.country.dlcurves}
\alias{e0.parDL.plot}

\title{
Plotting  Posterior Distribution of the Double Logistic Function of Life Expectancy
}
\description{
The functions plot the posterior distribution of the double logistic function used in the simulation, including their median and given probability intervals.
}
\usage{
e0.DLcurve.plot(mcmc.list, country, burnin = NULL, pi = 80, 
    e0.lim = NULL, nr.curves = 20, predictive.distr = FALSE, ylim = NULL, 
    xlab = "e(0)", ylab = "5-year gains", main = NULL, show.legend=TRUE, 
    col=c('black', 'red', "#00000020"), \dots)
    
e0.DLcurve.plot.all(mcmc.list = NULL, sim.dir = NULL, 
    output.dir = file.path(getwd(), "DLcurves"),
    output.type = "png", burnin = NULL, verbose = FALSE, \dots)
    
e0.parDL.plot(mcmc.set, country = NULL, burnin = NULL, lty = 2, 
	ann = TRUE, \dots)
	
e0.world.dlcurves(x, mcmc.list, burnin=NULL, \dots)

e0.country.dlcurves(x, mcmc.list, country, burnin=NULL, \dots)
}
\arguments{
  \item{mcmc.list}{List of \code{\link{bayesLife.mcmc}} objects, an object of class \code{\link{bayesLife.mcmc.set}} or of class \code{\link{bayesLife.prediction}}. In case of \code{e0.DLcurve.plot.all} if it si \code{NULL}, it is loaded from \code{sim.dir}.}
  \item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}} or \code{\link{bayesLife.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{burnin}{Number of iterations to be discarded from the beginning of parameter traces.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{e0.lim}{It can be a tuple of the minimum and maximum life expectancy to be shown in the plot. If \code{NULL}, it takes the minimum of observed data and 40, and the maximum of observed data and 90.}
  \item{nr.curves}{Number of curves to be plotted. If \code{NULL}, all curves are plotted.}
  \item{predictive.distr}{Logical. If \code{TRUE}, an error term is added to each trajectory.}
  \item{ylim, xlab, ylab, main, lty}{Graphical parameters passed to the \code{plot} function.}
  \item{show.legend}{Logical determining if the legend should be shown.}
  \item{col}{Vector of colors in this order: 1. observed data points, 2. quantiles, 3. trajectories}
  \item{\dots}{Additional graphical parameters. In addition, any arguments from \code{e0.DLcurve.plot} except \code{country} can be passed to \code{e0.DLcurve.plot.all}.}
  \item{sim.dir}{Directory with the simulation results. Only relevant, if \code{mcmc.list} is \code{NULL}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{x}{e0 values for which the double logistic should be computed.}
  \item{ann}{Logical if parameters should be annotated.}

}
\details{\code{e0.DLcurve.plot} plots double logistic curves for the given country.  \code{e0.DLcurve.plot.all} creates such plots for all countries and stores them in \code{output.dir}. Parameters passed to the double logistic function are either thinned traces created by the \code{\link{e0.predict}} function (if \code{mcmc.list} is an object of class \code{\link{bayesLife.prediction}}), or they are selected by equal spacing from the MCMC traces. In the former case, \code{burnin} is set automatically; in the latter case, \code{burnin} defaults to 0 since such object has already been \dQuote{burned}. If \code{nr.curves} is smaller than 2000, the median and probability intervals are computed on a sample of 2000 equally spaced data points, otherwise on all plotted curves.

Function \code{e0.parDL.plot} draws the means of the DL parameters as vertical and horizontal lines. The lines are added to the current graphical device and annotated if \code{ann} is \code{TRUE}. If country is \code{NULL}, the mean of world parameters are drawn.
 
Function \code{e0.world.dlcurves} returns the DL curves of the hierarchical distribution. Function \code{e0.country.dlcurves} returns DL curves for a given country. If \code{mcmc.list} is a prediction object, \code{burnin} should not be given, as such object has already been \dQuote{burned}.  
}
\value{
	\code{e0.world.dlcurves} and \code{e0.country.dlcurves} return a matrix of size \eqn{N \times M} where \eqn{N} is the number of trajectories and \eqn{M} is the number of values of \eqn{x}.
}

\author{
Hana Sevcikova
}

\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
mcmc.set <- get.e0.mcmc(sim.dir=sim.dir)
e0.DLcurve.plot(mcmc.set, country="Japan", burnin=40)
e0.parDL.plot(mcmc.set, "Japan")

# add the median of the hierarchical DL curves
x <- seq(40, 90, length=100)
world <- e0.world.dlcurves(x, mcmc.set, burnin=40)
qw <- apply(world, 2, median) 
lines(x, qw, col='blue')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

