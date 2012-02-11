\name{joint.male.predict}
\alias{joint.male.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction of the Joint Female-Male Model
}
\description{
Prediction of the joint female-male model of life expectancy, as described in Lalic (2011).
}
\usage{
joint.male.predict(e0.pred, estimates = NULL, gap.lim = c(0, 18), 
    my.e0.file = NULL, verbose = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e0.pred}{Object of class \code{\link{bayesLife.prediction}} containing female projections of life expectancy.}
  \item{estimates}{List of the same structure as returned by \code{\link{joint.male.estimate}}, containing the model estimation results. If it is not given, the \code{\link{joint.male.estimate}} function is invoked.}
  \item{gap.lim}{Vector of length two giving the minimum and maximum bounds for the female-male life expectancy gaps.}
  \item{my.e0.file}{File name containing user-specified male time series for one or more countries. The function replaces the corresponding country data from the UN WPP file by values in this file. Only columns are replaced that match column names of the WPP file.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Further arguments passed to \code{\link{joint.male.estimate}}.}
}
\details{
If no \code{estimates} are given, the function invokes an estimation by calling \code{\link{joint.male.estimate}}. Using those estimates, the male life expectancy is projected forward in time (as a function of a female-male gap), using the female predictions from \code{e0.pred}. The initial male data point is extracted from the male UN WPP data (\link{UN2010e0}) and possibly partly replaced by any user-specified data given in \code{my.e0.file}.

The resulting trajectory files are stored in a subdirectory of the female prediction directory, called \sQuote{joint_male}. Furthermore, an object of class \code{\link{bayesLife.prediction}} is created and added to \code{e0.pred} as a component called \code{joint.male}. 

The predicted gaps can be viewed using the \code{\link{e0.gap.plot}} function.
}
\value{
	Updated \code{e0.pred} object where a new component was added, called \code{joint_male}. It is also an object of class \code{\link{bayesLife.prediction}} and it contains results from this prediction.
}
\references{
Lalic, N. (2011). Master's thesis at the Department of Statistics, University of Washington. 
}
\author{
Hana Sevcikova
}

\seealso{
\code{\link{joint.male.estimate}}, \code{\link{e0.predict}}, \code{\link{e0.gap.plot}}
}
\examples{
sim.dir <- tempfile()
m <- run.e0.mcmc(sex='F', nr.chains=1, iter=30, thin=1, output.dir=sim.dir)
pred <- e0.predict(m, burnin=15, verbose=FALSE, save.as.ascii=0)
new.pred <- joint.male.predict(pred)
e0.trajectories.plot(new.pred, 'Guatemala', main='Female - Guatemala')
e0.trajectories.plot(new.pred, 'Guatemala', joint.male=TRUE, 
        main='Male - Guatemala')

# Both projections in one plot
e0.trajectories.plot(new.pred, 'Guatemala', col=rep('pink', 5), 
        nr.traj=0, pi=c(80), show.legend=FALSE)
e0.trajectories.plot(new.pred, 'Guatemala', col=rep('blue', 5), 
        joint.male=TRUE, add=TRUE, nr.traj=0, pi=c(80), show.legend=FALSE)
legend('topleft', legend=c('female', 'male'), col=c('pink', 'blue'), 
        lty=1, bty='n')

# Plotting the gaps
e0.gap.plot(new.pred, 'Guatemala')

unlink(sim.dir, recursive=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

