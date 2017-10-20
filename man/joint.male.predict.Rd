\name{e0.jmale.predict}
\alias{e0.jmale.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction of the Joint Female-Male Model
}
\description{
Prediction of the joint female-male model of life expectancy, as described in Raftery et al. (2014, 2012) and Lalic (2011).
}
\usage{
e0.jmale.predict(e0.pred, estimates = NULL, gap.lim = c(0, 18), 
    max.e0.eq1.pred = 86, my.e0.file = NULL, my.locations.file = NULL,
    save.as.ascii = 1000, verbose = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}} containing female projections of life expectancy.}
  \item{estimates}{List of the same structure as returned by \code{\link{e0.jmale.estimate}}, containing the model estimation results. If it is not given, the \code{\link{e0.jmale.estimate}} function is invoked.}
  \item{gap.lim}{Vector of length two giving the minimum and maximum bounds for the female-male life expectancy gaps.}
 \item{max.e0.eq1.pred}{Maximum female life expectancy for male projections handled by the first equation (parameter \eqn{M} in Equation (1) in Raftery et al. (2012)).}
  \item{my.e0.file}{File name containing user-specified male time series for one or more countries. The function replaces the corresponding country data from the WPP dataset by values in this file. Only columns are replaced that match column names of the WPP dataset.}
  \item{my.locations.file}{File name containing user-specified locations if different from the default \code{\link[wpp2017]{UNlocations}} dataset. It should be the same file as passed to \code{\link{run.e0.mcmc}} for female life expectancy.}
  \item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ASCII file, or \dQuote{all} in which case all trajectories are converted. It should be set to 0, if no conversion is desired.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Further arguments passed to \code{\link{e0.jmale.estimate}}.}
}
\details{
If no \code{estimates} are given, the function invokes an estimation by calling \code{\link{e0.jmale.estimate}}. Using those estimates, the male life expectancy is projected forward in time (as a function of a female-male gap), using the female predictions from \code{e0.pred}. The initial male data point is extracted from the male WPP dataset (\code{\link[wpp2017]{e0M}}) and possibly partly replaced by any user-specified data given in \code{my.e0.file}.

The resulting trajectory files are stored in a subdirectory of the female prediction directory, called \sQuote{joint_male}. Furthermore, an object of class \code{\link{bayesLife.prediction}} is created and added to \code{e0.pred} as a component called \code{joint.male}. 

The predicted gaps can be viewed using the \code{\link{e0.gap.plot}} function.
}
\value{
	Updated \code{e0.pred} object where a new component was added, called \code{joint.male}. It is also an object of class \code{\link{bayesLife.prediction}} and it contains results from this prediction.
}
\references{
A. E. Raftery,  N. Lalic, P. Gerland (2014). Joint Probabilistic Projection of Female and Male Life Expectancy. Demographic Research, 30:795-822.
	
A. E. Raftery,  N. Li, H. Sevcikova , P. Gerland, G. K. Heilig (2012). Bayesian probabilistic population projections for all countries. Proceedings of the National Academy of Sciences 109:13915-13921.
	
Lalic, N. (2011). Master's thesis at the Department of Statistics, University of Washington. 
}
\author{
Hana Sevcikova
}

\seealso{
\code{\link{e0.jmale.estimate}}, \code{\link{get.e0.jmale.prediction}}, \code{\link{e0.predict}}, \code{\link{e0.gap.plot}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.e0.mcmc(sex='F', nr.chains=1, iter=30, thin=1, output.dir=sim.dir)
pred <- e0.predict(m, burnin=15, verbose=FALSE, save.as.ascii=0)
both.pred <- e0.jmale.predict(pred)
e0.trajectories.plot(both.pred, 'Guatemala') # Female
e0.trajectories.plot(get.e0.jmale.prediction(both.pred), 'Guatemala') # Male

# Marginal distribution of the sex-specific projections
e0.trajectories.plot(both.pred, 'Guatemala', both.sexes=TRUE,  pi=80)

# Plotting the gaps
e0.gap.plot(both.pred, 'Guatemala')

# Joint distribution of the sex-specific projections
e0.joint.plot(both.pred, 'Guatemala', pi=80, years=c(2013, 2043, 2093))

unlink(sim.dir, recursive=TRUE)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

