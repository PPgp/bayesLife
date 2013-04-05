\name{summary.bayesLife.prediction}
\Rdversion{1.1}
\alias{summary.bayesLife.prediction}
\alias{print.summary.bayesLife.prediction}

\title{
Summary of a Prediction of the Life Expectancy
}
\description{
Country-specific summary of an object of class \code{\link{bayesLife.prediction}}, created using the function \code{\link{e0.predict}}. The summary contains the mean, standard deviation and several commonly used quantiles of the simulated trajectories.
}
\usage{
\method{summary}{bayesLife.prediction}(object, country = NULL, compact = TRUE, \dots)
}
\arguments{
  \item{object}{Object of class \code{\link{bayesLife.prediction}}.}
  \item{country}{Country name or code.}
  \item{compact}{Logical switching between a smaller and larger number of displayed quantiles.}
  \item{\dots}{Not used.}
}

\author{
Hana Sevcikova
}


\seealso{
\code{\link{bayesLife.prediction}}
}
\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
pred <- get.e0.prediction(sim.dir=sim.dir)                                                    
summary(pred, country="Iceland")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ print }
