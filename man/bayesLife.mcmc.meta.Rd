\name{bayesLife.mcmc.meta}
\alias{bayesLife.mcmc.meta}

\title{
MCMC Simulation Meta Object
}
\description{
Simulation meta object \code{bayesLife.mcmc.meta} used by all chains of the same MCMC simulation. 
It contains information that is common to all chains. It is a part of a \code{\link{bayesLife.mcmc.set}} object.
}


\value{
	A \code{bayesLife.mcmc.meta} object contains components \code{sex, nr.chains, 
		start.year, present.year, wpp.year, buffer.size, my.e0.file, a, delta, tau,}\cr
		\code{Triangle.ini, k.ini, z.ini, omega.ini, lambda.ini, lambda.k.ini,}\cr 
		\code{lambda.z.ini, Triangle.ini.low, 
		Triangle.ini.up, k.ini.low, k.ini.up, z.ini.low, z.ini.up, lambda.ini.low, lambda.ini.up, 
		lambda.k.ini.low, lambda.k.ini.up, lambda.z.ini.low, lambda.z.ini.up, omega.ini.low, omega.ini.up,
		Triangle.c.ini.norm, k.c.ini.norm, z.c.ini.norm,}\cr
		\code{Triangle.c.width, k.c.width, z.c.width, nu, dl.p1, dl.p2}. 
		Their meaning and value are the same as the input arguments 
			to the \code{\link{run.e0.mcmc}} function. Furthermore, it contains components:
	\item{e0.matrix.all}{A \eqn{q \times n}{q x n} matrix with the United Nations life expectancy estimates. \eqn{q} is number of years, \eqn{n} is number of countries (see \code{nr.countries} below). The first \eqn{n_e} columns correspond to countries included in the MCMC estimation (see \code{nr.countries.estimation} below), where \eqn{n_e <= n}. The order of the countries corresponds to the order of countries in the element \code{regions}, see below.}
	\item{e0.matrix}{Like \code{e0.matrix.all}, but it has \code{NA} values for years where no historical data is available.}
	\item{d.ct}{A difference e0 matrix of size \eqn{(q-1) \times n}{(q-1) x n}. Each element \eqn{x_{t,i}} is a difference 
	\eqn{y_{t+1,i}-y_{t,i}} where \eqn{y} denotes elements of \code{e0.matrix}.}
	\item{loessSD}{The loess matrix of \code{d.ct}.}
	\item{nr.countries}{Number of countries included in the e0 matrices.}
	\item{nr.countries.estimation}{Number of countries included in the MCMC estimation. It must be smaller or equal to \code{nr.countries}.}
	\item{T.end.c}{For each country, number of years for which United Nations historical data are available.}
	\item{regions}{List of arrays of length \code{nr.countries}. These are: \cr
		\code{name} - Region name for each country.\cr
		\code{code} - Region code for each country. \cr
		\code{area_name} - Area name for each country. \cr
		\code{area_code} - Area code for each country. \cr
		\code{country_name} - Array of country names. \cr
		\code{country_code} - Array of country codes. \cr
		Any country indices in the \code{bayesLife.mcmc.meta} object are derived from this component.}
	\item{output.dir}{Directory for storing simulation output.}
}

\details{
The object is in standard cases not to be manipulated by itself, but rather as part of a \code{\link{bayesLife.mcmc.set}} object.
}


\author{
Hana Sevcikova
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.e0.mcmc}}, \code{\link{get.e0.mcmc}}
}
\examples{
sim.dir <- file.path(.find.package("bayesLife"), "ex-data", 'bayesLife.output')
m <- get.e0.mcmc(sim.dir)
summary(m, meta.only = TRUE)
names(m$meta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }

