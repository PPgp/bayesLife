\name{get.e0.prediction}
\alias{get.e0.prediction}
\alias{has.e0.prediction}

\title{
Accessing a Prediction Object
}
\description{
Function \code{get.e0.prediction} retrieves results of a prediction and creates an object of class \code{\link{bayesLife.prediction}}. Function \code{has.e0.prediction} checks an existence of such results.
}
\usage{
get.e0.prediction(mcmc = NULL, sim.dir = NULL, mcmc.dir = NULL)

has.e0.prediction(mcmc = NULL, sim.dir = NULL)
}

\arguments{
  \item{mcmc}{Object of class \code{\link{bayesLife.mcmc.set}} used to make the prediction.  If it is \code{NULL}, the prediction is loaded from directory given by \code{sim.dir}.}
  \item{sim.dir}{Directory where the prediction is stored.  It should correspond to the value of the \code{output.dir} argument used in the \code{\link{e0.predict}} function.  Only relevant if \code{mcmc} is \code{NULL}.}
  \item{mcmc.dir}{Optional argument to be used only in a special case when the mcmc object contained in the prediction object was estimated in different directory than in the one to which it points to (for example due to moving or renaming the original directory). The argument causes that the mcmc is redirected to the given directory.}
}
\details{If \code{mcmc} is not \code{NULL}, the search directory is set to \code{mcmc$meta$output.dir}. This approach assumes that the prediction was stored in the same directory as the MCMC simulation, i.e. the \code{output.dir} argument of the \code{\link{e0.predict}} function was set to \code{NULL}. If it is not the case, the argument \code{mcmc.dir} should be used.}

\value{
	Function \code{has.e0.prediction} returns a logical indicating if a prediction exists for the given \code{mcmc}. 
	
	Function \code{get.e0.prediction} returns an
object of class \code{\link{bayesLife.prediction}}.}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesLife.prediction}}, \code{\link{e0.predict}}, \code{\link{summary.bayesLife.prediction}}
}
\examples{
sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
pred <- get.e0.prediction(sim.dir=sim.dir)
summary(pred, country='Canada')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

