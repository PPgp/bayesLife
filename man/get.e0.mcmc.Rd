\name{get.e0.mcmc}
\alias{get.e0.mcmc}
\alias{has.e0.mcmc}
\alias{e0.mcmc}
\alias{e0.mcmc.list}
\title{
Accessing MCMC Results
}
\description{
The function \code{get.e0.mcmc} retrieves results of an MCMC simulation and creates an object of class \code{\link{bayesLife.mcmc.set}}. Function \code{has.e0.mcmc} checks the existence of such results. Function \code{e0.mcmc} extracts a single chain and \code{e0.mcmc.list} extracts several or all chains from the simulation results. 

}
\usage{
get.e0.mcmc(sim.dir = file.path(getwd(), "bayesLife.output"), 
    chain.ids = NULL, low.memory = TRUE, burnin = 0, verbose = FALSE)
    
has.e0.mcmc(sim.dir)

e0.mcmc(mcmc.set, chain.id = 1)

e0.mcmc.list(mcmc.set, chain.ids=NULL)
}
\arguments{
  \item{sim.dir}{Directory where the simulation results are stored.}
  \item{chain.ids}{Chain identifiers in case only specific chains should be included in the resulting object. By default, all available chains are included.}
  \item{low.memory}{If \code{FALSE} full MCMC traces are loaded into memory.}
  \item{burnin}{Burnin used for loading traces. Only relevant, if \code{low.memory=FALSE}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{chain.id}{Chain identifier.}
  \item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{\code{get.e0.mcmc} returns an object of class \code{\link{bayesLife.mcmc.set}}. \code{has.e0.mcmc} returns a logical value.
	\code{e0.mcmc} returns an object of class \code{\link{bayesLife.mcmc}}, and \code{e0.mcmc.list} returns a list of \code{\link{bayesLife.mcmc}} objects. 
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesLife.mcmc.set}}
}

\examples{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
m <- get.e0.mcmc(sim.dir)
summary(m)

# summary of the world parameters for the single chains
for(mc in e0.mcmc.list(m)) print(summary(mc, par.names.cs=NULL))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }

