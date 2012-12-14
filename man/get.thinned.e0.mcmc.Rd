\name{get.thinned.e0.mcmc}
\alias{get.thinned.e0.mcmc}
\alias{create.thinned.e0.mcmc}

\title{
Creating and Accessing Thinned MCMCs
}
\description{
The function \code{get.thinned.e0.mcmc} accesses 
a thinned and burned version of the given MCMC set. \code{create.thinned.e0.mcmc} creates such set.
}
\usage{
get.thinned.e0.mcmc(mcmc.set, thin = 1, burnin = 0)

create.thinned.e0.mcmc(mcmc.set, thin = 1, burnin = 0, 
    output.dir = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}}.}
  \item{thin, burnin}{Thinning interval and burnin used for creating or identifying the thinned object.}
  \item{output.dir}{Directory for storing the thinned object. By default it is stored into the same directory as \code{mcmc.set}.}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
The function \code{create.thinned.e0.mcmc} is called from \code{\link{e0.predict}} and thus, the resulting object contains exactly the same MCMCs used for generating projections. 

The thinning is done as follows: The given \code{burnin} is removed from the beginning of each chain in the original MCMC set. Then each chain is thinned by \code{thin} using equal spacing and all chains are collapsed into one single chain per parameter. They are stored in \code{output.dir} under the name \file{thinned_mcmc_\emph{t}_\emph{b}} where \emph{t} is the value  of \code{thin} and \emph{b} the value of \code{burnin}.   
}
\value{
Both functions return an object of class \code{\link{bayesLife.mcmc.set}}. \code{get.thinned.e0.mcmc} returns \code{NULL} if such object does not exist.
}

\author{
Hana Sevcikova
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bayesLife.mcmc.set}}, \code{\link{e0.predict}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.e0.mcmc(nr.chains=2, iter=60, thin=2, output.dir=sim.dir, verbose=TRUE)
e0.predict(m, burnin=30) # creates thinned MCMCs
mb <- get.thinned.e0.mcmc(m, thin=2, burnin=30)
summary(mb, meta.only=TRUE) # length 30 = 2chains x (60-30)iters./2thin
unlink(sim.dir, recursive=TRUE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
\keyword{ manip }% __ONLY ONE__ keyword per line
