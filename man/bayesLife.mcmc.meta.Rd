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
	A \code{bayesLife.mcmc.meta} object contains components \code{sex}, \code{nr.chains}, 
		\code{start.year}, \code{present.year}, \code{wpp.year}, \code{buffer.size}, \code{my.e0.file}, \code{a}, \code{delta}, \code{tau},
		\code{Triangle.ini}, \code{k.ini}, \code{z.ini}, \code{omega.ini}, \code{lambda.ini}, \code{lambda.k.ini},
		\code{lambda.z.ini}, \code{Triangle.ini.low}, 
		\code{Triangle.ini.up}, \code{k.ini.low}, \code{k.ini.up}, \code{z.ini.low}, \code{z.ini.up}, \code{lambda.ini.low}, \code{lambda.ini.up}, 
		\code{lambda.k.ini.low}, \code{lambda.k.ini.up}, \code{lambda.z.ini.low}, \code{lambda.z.ini.up}, \code{omega.ini.low}, \code{omega.ini.up},
		\code{Triangle.c.ini.norm}, \code{k.c.ini.norm}, \code{z.c.ini.norm}, \code{Triangle.c.prior.low}, \code{Triangle.c.prior.up},
		\code{k.c.prior.low}, \code{k.c.prior.up}, \code{z.c.prior.low}, \code{z.c.prior.up},
		\code{Triangle.c.width}, \code{k.c.width}, \code{z.c.width}, \code{nu}, \code{dl.p1}, \code{dl.p2}, \code{sumTriangle.lim}, \code{auto.conf}. 
		Their meaning and value are the same as the input arguments 
			to the \code{\link{run.e0.mcmc}} function. Furthermore, it contains components:
	\item{e0.matrix.all}{A \eqn{q \times n}{q x n} matrix with the United Nations life expectancy estimates. \eqn{q} is number of years, \eqn{n} is number of countries (see \code{nr.countries} below). The first \eqn{n_e} columns correspond to countries included in the MCMC estimation (see \code{nr.countries.estimation} below), where \eqn{n_e <= n}. The order of the countries corresponds to the order of countries in the element \code{regions}, see below.}
	\item{e0.matrix}{Like \code{e0.matrix.all}, but it has \code{NA} values for years where no historical data is available.}
	\item{d.ct}{A difference e0 matrix of size \eqn{(q-1) \times n}{(q-1) x n}. Each element \eqn{x_{t,i}} is a difference 
	\eqn{y_{t+1,i}-y_{t,i}} where \eqn{y} denotes elements of \code{e0.matrix}.}
	\item{loessSD}{The loess matrix of \code{d.ct}.}
	\item{nr.countries}{Number of countries included in the e0 matrices.}
	\item{nr.countries.estimation}{Number of countries included in the MCMC estimation. It must be smaller or equal to \code{nr.countries}.}
	\item{Tc.index}{A list with one element per country. For each country, it contains the index within \code{e0.matrix} where the observed historical data are not \code{NA}, i.e. it points to the data used in the estimation.}
	\item{regions}{List of arrays of length \code{nr.countries}. These are: \cr
		\code{name} - Region name for each country.\cr
		\code{code} - Region code for each country. \cr
		\code{area_name} - Area name for each country. \cr
		\code{area_code} - Area code for each country. \cr
		\code{country_name} - Array of country names. \cr
		\code{country_code} - Array of country codes. \cr
		Any country indices in the \code{bayesLife.mcmc.meta} object are derived from this component.}
	\item{output.dir}{Directory for storing simulation output.}
	\item{suppl.data}{If supplemental data were used in the simulation (i.e. start year was set prior to 1950), this is a list containing information about the additional data. It has the following components of the same form as described above, but related only to the additional data: \code{e0.matrix},  \code{regions}, \code{Tc.index}, \code{nr.countries}, \code{d.ct}, \code{loessSD}. In addition, it has the vectors\cr
		\code{index.from.all.countries} - a vector of size \code{nr.countries} (all countries) where each element is an index to the supplemental data, i.e. a match from regions$country_code to suppl.data$regions$country_code.\cr
		\code{index.to.all.countries} - a vector of size \code{suppl.data$nr.countries} (additional countries) where each element is an index to all data, i.e. a match from suppl.data$regions$country_code to regions$country_code.
}
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
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
m <- get.e0.mcmc(sim.dir)
summary(m, meta.only = TRUE)
names(m$meta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }

