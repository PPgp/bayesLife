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
e0.trajectories.plot(e0.pred, country, pi = c(80, 95), both.sexes = FALSE,
    nr.traj = NULL, adjusted.only = TRUE, typical.trajectory = FALSE, 
    xlim = NULL, ylim = NULL, type = "b", 
    xlab = "Year", ylab = "Life expectancy at birth", main = NULL, 
    lwd = c(2, 2, 2, 2, 1), col = c('black', 'green', 'red', 'red', '#00000020'),
    col2 = c('gray39', 'greenyellow', 'hotpink', 'hotpink', '#00000020'),
    show.legend = TRUE, add = FALSE, \dots)
    
e0.trajectories.plot.all(e0.pred, 
    output.dir = file.path(getwd(), 'e0trajectories'),
    output.type = "png", verbose = FALSE, \dots)
	
e0.trajectories.table(e0.pred, country, pi = c(80, 95), 
    both.sexes = FALSE, \dots)
}

\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{pi}{Probability interval. It can be a single number or an array. If \code{both.sexes} is \code{TRUE} the default is 95.}
  \item{both.sexes}{Logical or the character \dQuote{A}. If \code{TRUE} the distribution of both sexes is plotted into one graphics (or tabulated), provided the \code{e0.pred} is a female prediction and contains a joint male prediction as a result of running the function \code{\link{e0.jmale.predict}}. For \dQuote{A} it plots/tabulates the distribution of the average life expectancy over both sexes.}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly. If \code{both.sexes} is \code{TRUE} the default is zero.}
  \item{adjusted.only}{Logical. By default, if the projection median is adjusted using e.g. \code{\link{e0.median.set}}, the function plots the adjusted median. If \code{adjusted.only=FALSE} the original (non-adjusted) median is plotted as well.}
  \item{typical.trajectory}{Logical. If \code{TRUE} one trajectory is shown  for which the median absolute deviation from the median e0 projection is the median among all the trajectories.}
  \item{xlim, ylim, type, xlab, ylab, main}{Graphical parameters passed to the \code{plot} function.}
  \item{lwd, col, col2}{Vector of five elements giving the line width and color for: 1. observed data, 2. imputed missing data, 3. median, 4. quantiles, 5. trajectories. \code{col2} is only used if \code{both.sexes} is \code{TRUE}. In such a case, \code{col2} is used for female lines and \code{col} is used for male lines, which in this case defaults to 
  	\code{c('black', 'green', 'darkgreen', 'darkgreen', 'gray')}.}
  \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{add}{Logical controlling whether the trajectories should be plotted into a new graphic device (\code{FALSE}) or into an existing device (\code{TRUE}). One can use this argument to plot trajectories from multiple countries into one graphics.}
  \item{\dots}{Additional graphical parameters. In addition, for \code{e0.trajectories.plot.all}, \code{\dots} 
  	contains any of the arguments of  \code{e0.trajectories.plot}, and for \code{e0.trajectories.table}, \code{\dots} contains 
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
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
pred <- get.e0.prediction(sim.dir)
e0.trajectories.table(pred, country="Japan", pi=c(80, 95))
e0.trajectories.plot(pred, country="Japan", pi=c(80, 95))

# plot multiple countries into one plot
e0.trajectories.plot(pred, "Japan", col=rep("green", 5), nr.traj=0, pi=c(80), 
        show.legend=FALSE, main="")
e0.trajectories.plot(pred, "United States of America", col=rep("blue", 5),  
        add=TRUE, nr.traj=0, pi=c(80), show.legend=FALSE)
legend("topleft", legend=c("Japan", "USA"), col=c("green", "blue"), lty=1, bty="n")
}}

\keyword{ hplot }

