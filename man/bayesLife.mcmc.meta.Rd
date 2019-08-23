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
	A \code{bayesLife.mcmc.meta} object stores values of the various input arguments 
			of the \code{\link{run.e0.mcmc}} function. These are \code{sex}, \code{nr.chains}, 
		\code{start.year}, \code{present.year}, \code{wpp.year}, \code{my.e0.file}, \code{compression.type}. \cr
		Furthermore, it contains components:
	\item{e0.matrix.all}{A \eqn{q \times n}{q x n} matrix with the United Nations life expectancy estimates. \eqn{q} is number of years, \eqn{n} is number of countries (see \code{nr.countries} below). The first \eqn{n_e} columns correspond to countries included in the MCMC estimation (see \code{nr.countries.estimation} below), where \eqn{n_e <= n}. The order of the countries corresponds to the order of countries in the element \code{regions} and \code{regionsDT}, see below.}
	\item{e0.matrix}{Like \code{e0.matrix.all}, but it has \code{NA} values for years that were excluded from the estimation due to setting the \code{start.year} and \code{present.year} argument, and the \sQuote{last.observed} column (see Details in \code{\link{run.e0.mcmc}}).}
	\item{d.ct}{A difference e0 matrix of size \eqn{(q-1) \times n}{(q-1) x n}. Each element \eqn{x_{t,i}} is a difference 
	\eqn{y_{t+1,i}-y_{t,i}} where \eqn{y} denotes elements of the \code{e0.matrix}. Outliers are set to \code{NA}.}
	\item{loessSD}{Matrix of the same dimension as \code{d.ct}. It contains values of the \code{loess_sd} dataset at the e0 levels of \code{e0.matrix}. See also \code{\link{compute.loess}}.}
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
	\item{regionsDT}{Like \code{regions} but for convenience it is a \code{data.table} object. A country in i-th row of this table has e0 data stored in the i-th column of \code{e0.matrix}.}
	\item{output.dir}{Directory for storing simulation output.}
	\item{mcmc.options}{List of various options used in the estimation. See \code{\link{e0mcmc.options}}.}
	\item{country.bounds}{List of country-specific bounds of the various parameter priors, as constructed from the \code{country.overwrites} option (see \code{\link{e0mcmc.options}}).}
	\item{suppl.data}{If supplemental data were used in the simulation (i.e. start year was set prior to 1950), this is a list containing information about the additional data. It has the following components of the same form as described above, but related only to the additional data: \code{e0.matrix},  \code{regions}, \code{regionsDT}, \code{Tc.index}, \code{nr.countries}, \code{d.ct}, \code{loessSD}. In addition, it has the vectors\cr
		\code{index.from.all.countries} - a vector of size \code{nr.countries} (all countries) where each element is an index to the supplemental data, i.e. a match from regions$country_code to suppl.data$regions$country_code.\cr
		\code{index.to.all.countries} - a vector of size \code{suppl.data$nr.countries} (additional countries) where each element is an index to all data, i.e. a match from \code{suppl.data$regions$country_code} to \code{regions$country_code}.
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
\code{\link{run.e0.mcmc}}, \code{\link{get.e0.mcmc}}, \code{\link{e0mcmc.options}}
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

