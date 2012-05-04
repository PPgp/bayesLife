\name{e0.predict.extra}
\alias{e0.predict.extra}

\title{
Generating Posterior Trajectories of the Life Expectancy for Specific Countries or Regions
}
\description{
Using the posterior parameter samples the function generates posterior trajectories of the life expectancy for given countries or regions. It is intended to be used after running \code{\link{run.e0.mcmc.extra}}, but it can be also used for purposes of testing specific settings on one or a few countries. 
}
\usage{
e0.predict.extra(sim.dir = file.path(getwd(), 'bayesLife.output'), 
    prediction.dir = sim.dir, countries = NULL, 
    save.as.ascii = 1000, verbose = TRUE, \dots)
	}
\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results.}
  \item{prediction.dir}{Directory where the prediction object and the trajectories are stored.}
  \item{countries}{Vector of country codes for which the prediction should be made. If it is \code{NULL}, the prediction is run for all countries that are included in the MCMC object but for which no prediction was generated.}
  \item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ascii file, or \dQuote{all} in which case all trajectories are converted. It should be set to 0, if no converions is desired. Note that the convertion is done on all countries.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Additional arguments passed to a joint female-male prediction.}
  }
\details{
In order to use this function, a prediction object must exist, i.e. the function \code{\link{e0.predict}} must have been processed prior to using this function.

Trajectories for given countries or regions are generated and stored in binary format along with other countries (in \code{prediction.dir}). The existing prediction object is updated and stored in the same directory. If \code{save.as.ascii} is larger than zero, trajectories of ALL countries are converted to an ascii format.

If the prediction object contains joint male projections, these are also created for the given countries.
}
\value{
	Updated object of class \code{\link{bayesLife.prediction}}.
	}
	
\author{
Hana Sevcikova
}
\seealso{
\code{\link{e0.predict}}, \code{\link{run.e0.mcmc.extra}}
}
%\examples{

%}
\keyword{distribution}
\keyword{multivariate}
