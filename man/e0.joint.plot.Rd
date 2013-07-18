\name{e0.joint.plot}
\alias{e0.joint.plot}
\alias{e0.joint.plot.all}

\title{
Joint Posterior Distribution of Female and Male Life Expectancy
}
\description{
The functions plot the joint posterior distribution of female and male life expectancy, modeled and predicted using the joint model described in Lalic (2011). This can be done for one country (\code{e0.joint.plot}) or for all countries (\code{e0.joint.plot.all}).}

\usage{
e0.joint.plot(e0.pred, country, pi = 95, years, nr.points = 500, 
    obs.pch = 17, obs.cex=1, xlim = NULL, 
    ylim = NULL, xlab = "Female life expectancy", ylab = "Male life expectancy", 
    main = NULL, col = NULL, show.legend = TRUE, add = FALSE, \dots)
    
e0.joint.plot.all(e0.pred, output.dir = file.path(getwd(), "e0joint"),
    output.type = "png", verbose = FALSE, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}} containing female projections.}
  \item{country}{Name or numerical code of a country.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{years}{Array of future years for which to plot the distribution.}
  \item{nr.points}{Number of points shown in the plot for each year.}
  \item{obs.pch, obs.cex}{Graphical parameters used for displaying observed data or data without variation.}
  \item{xlim, ylim, xlab, ylab, main}{Graphical parameters passed to the \code{\link{plot}} function.}
  \item{col}{Array of colors, one for each year.}
  \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{add}{Logical controlling whether the distribution should be plotted into a new graphic device (\code{FALSE}) or into an existing device (\code{TRUE}).}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Additional graphical parameters passed to the \code{\link{lines}} function drawing ellipses. In addition, for \code{e0.joint.plot.all}, \code{\dots} contains any of the arguments of  \code{e0.joint.plot}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\references{
Lalic, N. (2011). Master's thesis at the Department of Statistics, University of Washington.
}
\author{
Hana Sevcikova, Adrian Raftery
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e0.gap.plot}}, \code{\link{e0.trajectories.plot}}, \code{\link{e0.jmale.predict}}
}
\examples{
# See example for e0.jmale.predict
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

