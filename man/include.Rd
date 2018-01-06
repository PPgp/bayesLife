\name{include}
\alias{include}
\alias{include_2017}
\alias{include_2015}
\alias{include_2010}
\alias{include_2012}

\docType{data}
\title{
Inclusion Codes
}
\description{
Data sets containing codes that determine which countries are to be included into a simulation or/and projections. 
}
\usage{
data(include_2017)
data(include_2015)
data(include_2012)
data(include_2010)
}

\format{
  Data frames containing one record per country or region. It has the following variables:
  \describe{
    \item{\code{country}}{Name of country or region. Not used.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{include_code}{Entries for which \code{include_code=2} are included in MCMC simulations (i.e. estimation of the model parameters). Entries for which \code{include_code} is 1 or 2 are included in the prediction.}}
}
\details{
	In a simulation, an \code{include_*} dataset is selected that corresponds to the given \code{wpp.year} passed to the function \code{\link{run.e0.mcmc}}. It is merged with a \code{\link[wpp2017]{e0}} dataset from the corresponding wpp package using the \code{country_code} column. Thus, the country entries in this dataset should correspond to entries in the \code{e0F} (\code{e0M}) dataset.

The package contains also a dataset called \file{my_e0_template} (in \file{extdata} directory) which is a template for user-specified e0 time series. It has the same structure as the \code{\link[wpp2017]{e0}} dataset, except that most of the columns are optional. The only required column is \code{country_code} (see description of the argument \code{my.e0.file} in \code{\link{run.e0.mcmc}}).
}
\note{
In all three datasets, countries affected by AIDS are not included in the estimation, i.e. the \code{include_code} is set to 3.  
}

\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(include_2017)
head(include_2017)
# select AIDS countries
subset(include_2017, include_code == 3)
}
\keyword{datasets}
