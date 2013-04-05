\name{get.e0.trajectories}
\alias{get.e0.trajectories}

\title{
Accessing Trajectories of Life Expectancy
}
\description{
Function for accessing trajectories of the life expectancy.
}
\usage{
get.e0.trajectories(e0.pred, country)
}

\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}}.}
  \item{country}{Name or numerical code of a country.}
}
\details{
The function loads trajectories of life expectancy for the given country from disk
%, offsets it if needed (see \code{\link{e0.median.shift}}) 
and returns it as a matrix.
}
\value{
Array of size the number of projection periods (including the present year) times the number of trajectories. The row names correspond to the mid-years of the prediction periods.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesLife.prediction}}, \code{\link{get.e0.prediction}}, \code{\link{e0.trajectories.table}}
%, \code{\link{e0.median.shift}}
}
\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output") 
pred <- get.e0.prediction(sim.dir=sim.dir)
get.e0.trajectories(pred, "Germany")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
