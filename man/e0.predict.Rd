\name{e0.predict}
\alias{e0.predict}
\alias{bayesLife.prediction}

\title{
Generating Posterior Trajectories of the Life Expectancy
}
\description{
Using the posterior parameter samples simulated by \code{\link{run.e0.mcmc}} the function  generates posterior trajectories for the life expectancy for all countries of the world.
}
\usage{
e0.predict(mcmc.set = NULL, end.year = 2100, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    replace.output = FALSE, nr.traj = NULL, thin = NULL, 
    burnin = 2e+04, save.as.ascii = 1000, 
    output.dir = NULL, low.memory = TRUE, 
    seed = NULL, verbose = TRUE)
}

\arguments{
  \item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}}. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.
}
  \item{end.year}{End year of the prediction.}
  \item{sim.dir}{Directory with the MCMC simulation results. It should equal to the \code{output.dir} argument in \code{\link{run.e0.mcmc}}.}
  \item{replace.output}{Logical. If \code{TRUE}, existing predictions in \code{output.dir} will be replaced by results of this run.}
  \item{nr.traj}{Number of trajectories to be generated. If \code{NULL}, the argument \code{thin} is taken to determine the number of trajectories. If both are \code{NULL}, the number of trajectories corresponds to the minimum of the size of the parameter sample and 2000.}
  \item{thin}{Thinning interval used for determining the number of trajectories. Only relevant, if \code{nr.traj} is \code{NULL}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of the parameter traces.}
  \item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ASCII file, or \dQuote{all} in which case all trajectories are converted. It should be set to 0, if no conversion is desired.}
  \item{output.dir}{Directory into which the resulting prediction object and the trajectories are stored. If it is \code{NULL}, it is set to either \code{sim.dir}, or to \code{output.dir} of \code{mcmc.set$meta} if \code{mcmc.set} is given.}
  \item{low.memory}{Logical indicating if the prediction should run in a low-memory mode. If it is \code{FALSE}, the whole traces of all parameters, including the burnin, are loaded into memory. Otherwise, burnins are discarded and parameters are loaded as they are needed and are not kept in the memory.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible projections.}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
The trajectories are generated using the double logistic function (Chunn et al. 2010). Parameter samples  simulated via \code{\link{run.e0.mcmc}} are used from all chains, from which the given burnin was discarded. They are evenly thinned to match \code{nr.traj} or using the \code{thin} argument. Such thinned parameter traces, collapsed into one chain, if they do not already exist, are stored on disk into the sub-directory \file{\{thinned_mcmc_\emph{t}_\emph{b}} where \emph{t} is the value  of \code{thin} and \emph{b} the value of \code{burnin} (see \code{\link{create.thinned.e0.mcmc}}). 
	
%The projection is run for all missing values before the present year, if any. Medians over the trajectories are used as  imputed values and the trajectories are discarded. The process then continues by projecting the future values where all generated trajectories are kept.

The resulting prediction object is saved into \file{\{output.dir\}/predictions}. Trajectories for all countries are saved into the same directory in a binary format, one file per country. At the end of the projection, if \code{save.as.ascii} is larger than 0, the function converts the given number of trajectories into a CSV file of a UN-specific format. They are selected by equal spacing (see function \code{\link{convert.e0.trajectories}} for more details on the conversion). %In addition, two summary files are created: one in a user-friendly format, the other using a UN-specific coding of the variants and time (see \code{\link{write.projection.summary}} for more details).

}
\value{
	Object of class \code{bayesLife.prediction} which is a list containing components:
	\item{quantiles}{A \eqn{n \times q \times p}{n x q x p} array of quantile values computed on the trajectories. \eqn{n} is the number of countries,
		\eqn{q} is the number of quantile bounds and \eqn{p} is the number of projections.}
	\item{traj.mean.sd}{A \eqn{n \times 2 \times p}{n x 2 x p} array holding the mean of all trajectories in the first column and the standard deviation in the second column.
		\eqn{n} and \eqn{p} are the number of countries and number of projections, respectively.}
	\item{nr.traj}{Number of trajectories.}
	%\item{trf_matrix_reconstructed}{Matrix containing imputed TFR values on spots where the original TFR %matrix has missing values, i.e. between the last observed data point and the present year. }
	\item{output.directory}{Directory where trajectories corresponding to this prediction are stored.}
	\item{nr.projections}{Number of projections.}
	\item{burnin}{Burnin used for this prediction.}
	\item{end.year}{The end year of this prediction.}
	\item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}} used for this prediction, i.e. the burned, thinned, and collapsed MCMC chain.}
}
\references{
J. L. Chunn, A. E. Raftery, P. Gerland (2010): Bayesian Probabilistic Projections of Life Expectancy for All Countries.  Working Paper nr. 105, Center for Statistics and the Social Sciences, University of Washington. \url{http://www.csss.washington.edu/Papers}


}
\author{
Hana Sevcikova, using code from Jennifer Chunn
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.e0.mcmc}}, \code{\link{create.thinned.e0.mcmc}}, \code{\link{convert.e0.trajectories}}, 
\code{\link{get.e0.prediction}}, \code{\link{summary.bayesLife.prediction}}
}
\examples{
\dontrun{
m <- run.e0.mcmc(nr.chains=1, iter=50, thin=1, verbose=TRUE)
pred <- e0.predict(m, burnin=25, verbose=TRUE)
summary(pred, country='Portugal')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}