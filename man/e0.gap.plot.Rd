\name{e0.gap.plot}
\alias{e0.gap.plot}
\alias{e0.gap.plot.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Distribution of Gaps in Female and Male Life Expectancy
}
\description{
The functions plot the posterior distribution of the gaps between female and male life expectancy, modeled and predicted using a model described in Lalic (2011). This can be done for one country (\code{e0.gap.plot}) or for all countries (\code{e0.gap.plot.all}).
}
\usage{
e0.gap.plot(e0.pred, country, e0.pred2 = NULL, pi = c(80, 95), 
    nr.traj = 0, xlim = NULL, ylim = NULL, type = "b", 
    xlab = "Year", ylab = "Gap in life expectancy", main = NULL, 
    show.legend = TRUE, \dots)
    
e0.gap.plot.all(e0.pred, output.dir = file.path(getwd(), "e0gaps"),
    output.type = "png", verbose = FALSE, \dots)
}

\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}} containing female projections. If \code{e0.pred2} is not given, then this object must contain the \sQuote{joint.male} component, generated using the \code{\link{e0.jmale.predict}} function.}
  \item{country}{Name or numerical code of a country.}
  \item{e0.pred2}{Object of class \code{\link{bayesLife.prediction}} containing male projections. By default the male projections is taken from the joint female-male projections in \code{e0.pred}, see above.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{nr.traj}{Number of trajectories to be plotted.}
  \item{xlim, ylim, type, xlab, ylab, main}{Graphical parameters passed to the \code{plot} function.}
  \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Additional graphical parameters. In addition, for \code{e0.gap.plot.all}, \code{\dots} 
  	contains any of the arguments of  \code{e0.gap.plot}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\references{
Lalic, N. (2011). Master's thesis at the Department of Statistics, University of Washington.
}
\author{
Hana Sevcikova
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e0.joint.plot}}, \code{\link{e0.jmale.estimate}}, \code{\link{e0.jmale.predict}}, \code{\link{get.e0.jmale.prediction}}
}
\examples{
# See example for e0.jmale.predict
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

