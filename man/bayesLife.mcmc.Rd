\name{bayesLife.mcmc}
\alias{bayesLife.mcmc}

\title{
MCMC Simulation Object
}
\description{
MCMC simulation object \code{bayesLife.mcmc} containing information about one MCMC chain. A set of such objects belonging to the same simulation together with a \code{\link{bayesLife.mcmc.meta}} object constitute a \code{\link{bayesLife.mcmc.set}} object. 
}

\value{
A \code{bayesLife.mcmc} object contains parameters of the Bayesian hierarchical model, more specifically, their initial values (all names with the suffix \code{.ini}) and values from the last iteration.  These are: \cr
\code{Triangle/Triangle.ini, lambda/lambda.ini} - world parameters, containing four values each. They correspond to model parameters \eqn{\Delta_1, \dots, \Delta_4}{Delta_1, \dots, Delta_4} and \eqn{\lambda_1, \dots \lambda_4}{lambda_1, \dots lambda_4}, respectively.\cr
\code{k/k.ini, z/z.ini, omega/omega.ini, lambda.k/lambda.k.ini,}\cr \code{lambda.z/lambda.z.ini} -  world parameters, containing one value each. They correspond to model parameters \eqn{k}, \eqn{z}, \eqn{\omega}{omega}, \eqn{\lambda_k}{lambda_k}, and \eqn{\lambda_z}{lambda_z}, respectively.\cr
\code{Triangle.c} - country-specific parameter \eqn{\Delta^c_1, \dots, \Delta^c_4}{Delta^c_1, \dots, Delta^c_4} with four values for each country, i.e. an \eqn{4 \times C}{4 x C} matrix where \eqn{C} is the number of countries.\cr
\code{k.c, z.c} - country-specific parameters \eqn{k^c} and \eqn{z^c} (1d arrays of length \eqn{C}).\cr
Furthermore, the object contains components:
\item{iter}{Total number of iterations the simulation was started with.}
\item{finished.iter}{Number of iterations that were finished.  Results from the last finished iteration are stored in the parameters above.}
\item{length}{Length of the MCMC stored on disk. It differs from \code{finished.iter} only if \code{thin} is larger than one.}
\item{thin}{Thinning interval used when simulating the MCMCs.}
\item{id}{Identifier of this chain.}
\item{output.dir}{Subdirectory (relative to \code{output.dir} in the \code{\link{bayesLife.mcmc.meta}} object) where results of this chain are stored.}
\item{traces}{This is a placeholder for keeping whole parameter traces in the memory.  If the processing operates in a low memory mode, it will be 0.  It can be filled in using the function \code{\link{get.e0.mcmc}(\dots, low.memory=FALSE)}. In such a case, \code{traces} is a \eqn{I \times J}{I x J} array where \eqn{I} is the MCMC \code{length} and \eqn{J} is the number of parameters.}
\item{traces.burnin}{Burnin used to retrieve the traces, i.e. how many stored iterations are missing from the beginning in the \code{traces} array comparing to the \sQuote{raw} traces on the disk.}
\item{rng.state}{State of the random number generator at the end of the last finished interation.}
\item{meta}{Object of class \code{\link{bayesLife.mcmc.meta}} used for simulation of this chain.}
}
\details{
An object \code{bayesLife.mcmc} points to a place on disk (element \code{output.dir}) where MCMC results from all iterations are stored. They can be retrieved to the memory using \code{\link{get.e0.mcmc}(\dots)}. 

The object is in standard cases not to be manipulated by itself, but rather as part of a \code{\link{bayesLife.mcmc.set}} object.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{run.e0.mcmc}}, \code{\link{get.e0.mcmc}}, \code{\link{bayesLife.mcmc.set}}, \code{\link{bayesLife.mcmc.meta}}
}
\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
# loads traces from one chain
m <- get.e0.mcmc(sim.dir, low.memory=FALSE, burnin=35, chain.ids=1)
# should have 13 rows, since 60 iterations in total minus 35 burnin with thin=2
dim(e0.mcmc(m, 1)$traces)
summary(m)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

