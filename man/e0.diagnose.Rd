\name{e0.diagnose}
\alias{e0.diagnose}
\alias{bayesLife.convergence}

\title{
Convergence Diagnostics for Markov Chain Monte Carlo of Life Expectancy
}
\description{
Function \code{e0.diagnose} runs convergence diagnostics of existing MCMCs, using the \code{raftery.diag} function from the \pkg{coda} package. 
}
\usage{
e0.diagnose(sim.dir, thin=225, burnin = 10000, express = FALSE, 
    country.sampling.prop = NULL, keep.thin.mcmc=FALSE, verbose = TRUE)
}
\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results.}
  \item{thin}{Thinning interval.}
  \item{burnin}{Number of iterations to be discarded from the beginning of the parameter traces.}
  \item{express}{Logical. If \code{TRUE}, the convergence diagnostics is run only on the country-independent parameters. If \code{FALSE}, the country-specific parameters are included in the diagnostics. The number of countries can be controlled by \code{country.sampling.prop}.}
  \item{country.sampling.prop}{Proportion of countries that are included in the diagnostics. If it is \code{NULL} and \code{express=FALSE}, all countries are included. Setting here a number between 0 and 1, one can limit the number of countries which are then randomly sampled. Note that for long MCMCs, this argument may significantly influence the run-time of this function.}
  \item{keep.thin.mcmc}{Logical. If \code{TRUE} the thinned traces used for computing the diagnostics are stored on disk (see \code{\link{create.thinned.e0.mcmc}}).}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
The function invokes the \code{\link{e0.raftery.diag}} function separately for country-independent parameters and for country-specific parameters. It results in two possible states: red, i.e. it did not converge, and green, i.e. it converged.
The resulting object is stored in \cr \file{\{sim.dir\}/diagnostics/bayesLife.convergence_\{thin\}_\{burnin\}.rda} and can be accessed using the function \code{\link{get.e0.convergence}}.

Function \code{\link[bayesTFR]{has.mcmc.converged}} from the \pkg{bayesTFR} package can be used to check if the existing diagnostics converged.
}
\value{
	\code{e0.diagnose} returns an object of class \code{bayesLife.convergence} with components:
	\item{result}{Table containing all not-converged parameters. Its columns include \sQuote{Total iterations needed} and \sQuote{Remaining iterations}.}
	\item{lresult.country.independent}{Number of rows in \code{result} that correspond to country-independent paramters. These rows are groupped at the beginning of the table.}
	\item{country.independent}{Result of \code{\link{e0.raftery.diag}} processed on country-independent parameters.}
	\item{country.specific}{Result of \code{\link{e0.raftery.diag}} processed on country-specific parameters.}
	\item{iter.needed}{Number of additional iterations suggested in order to achieve convergence.}
	\item{iter.total}{Total number of iterations of the original unthinned set of chains.}
	\item{use.nr.traj}{Suggestion for number of trajectories in generating predictions.}
	\item{burnin}{Burnin used.}
	\item{thin}{Thinning interval used.}
	\item{status}{Vector of character strings containing the result status. Possible values: \sQuote{green}, \sQuote{red}.}
	\item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}} that corresponds to the original set of MCMCs on which the diagnostics was run.}
	\item{thin.mcmc}{If \code{keep.thin.mcmc} is \code{TRUE}, it is an object of class \code{\link{bayesLife.mcmc.set}} that corresponds to the thinned mcmc set on which the diagnostics was run, otherwise \code{NULL}.}
	\item{express}{Value of the input argument \code{express}.}
	\item{nr.countries}{Vector with elements \code{used} - number of countries used in this diagnostics, and \code{total} - number of countries that this \code{mcmc.set} object was estimated on.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Adrian Raftery
}
%\note{
%%  ~~further notes~~
%}


\seealso{
\code{\link{e0.raftery.diag}}, \code{\link[coda]{raftery.diag}}, \code{\link{summary.bayesLife.convergence}}, \code{\link{get.e0.convergence}}, \code{\link{create.thinned.e0.mcmc}}
}

\keyword{ htest }

