\name{e0.trajectories.plot}
\alias{e0.trajectories.plot}
\alias{e0.trajectories.plot.all}
\alias{e0.trajectories.table}

\title{
Posterior Distribution of Trajectories of Life Expectancy
}
\description{
The functions plot/tabulate the posterior distribution of trajectories of the life expectancy for a given country, or for all countries, including their median and given probability intervals. 
}
\usage{
e0.trajectories.plot(e0.pred, country, pi = c(80, 95), 
    nr.traj = NULL, xlim = NULL, ylim = NULL, 
    type = "b", xlab = "Year", ylab = "Life expectancy", 
    main = NULL, \dots)
    
e0.trajectories.plot.all(e0.pred, 
    output.dir = file.path(getwd(), 'e0trajectories'),
    output.type = "png", verbose = FALSE, \dots)
	
e0.trajectories.table(e0.pred, \dots)
}

\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{xlim, ylim, type, xlab, ylab, main}{Graphical parameters passed to the \code{plot} function.}
  \item{\dots}{Additional graphical parameters. In addition, for \code{e0.trajectories.plot.all}, \code{\dots} 
  	contains the arguments \code{pi} and \code{nr.traj}, and for \code{e0.trajectories.table}, \code{\dots} contains 
  	the \code{pi} and \code{country} arguments.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
\code{e0.trajectories.plot} plots posterior distribution of trajectories of life expectancy for a given country. \code{e0.trajectories.table} gives the same output in a tabular format. \cr
\code{e0.trajectories.plot.all} creates a set of such graphs (one per country) that are stored in \code{output.dir}. 

The median and given probability intervals are computed using all available trajectories. Thus, \code{nr.traj} does not influence those values - it is used only to control the number of trajectories plotted.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesLife.prediction}}
}

\examples{

sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
pred <- get.e0.prediction(sim.dir)
e0.trajectories.table(pred, country='Japan', pi=c(80, 95))
\dontrun{
e0.trajectories.plot(pred, country='Japan', pi=c(80, 95))
}
}

\keyword{ hplot }

