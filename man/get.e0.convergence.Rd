\name{get.e0.convergence}
\alias{get.e0.convergence}
\alias{get.e0.convergence.all}

\title{
Accessing a Convergence Object
}
\description{
The functions load objects of class \code{\link{bayesLife.convergence}} from disk that were created using the function \code{\link{e0.diagnose}}.
}
\usage{
get.e0.convergence(sim.dir = file.path(getwd(), "bayesLife.output"), 
    thin = 225, burnin = 10000)
	
get.e0.convergence.all(sim.dir = file.path(getwd(), "bayesLife.output"))
}
\arguments{
  \item{sim.dir}{Simulation directory used for computing the diagnostics.}
  \item{thin}{Thinning interval used with this diagnostics.}
  \item{burnin}{Burnin used for computing the diagnostics.}
}
\details{
Function \code{get.e0.convergence} loads an object of class \code{\link{bayesLife.convergence}} for the specific \code{thin} and \code{burnin}. Function \code{get.e0.convergence.all} loads all \code{\link{bayesLife.convergence}} objects available in \code{sim.dir}. 
}
\value{
\code{get.e0.convergence} returns an object of class \code{\link{bayesLife.convergence}}; \cr
\code{get.e0.convergence.all} returns a list of objects of class \code{\link{bayesLife.convergence}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{e0.diagnose}}, \code{\link{summary.bayesLife.convergence}}.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

