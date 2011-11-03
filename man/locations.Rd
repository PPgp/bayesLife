\name{WPP2008_LOCATIONS}
\alias{WPP2008_LOCATIONS}
\alias{WPP2010_LOCATIONS}
\docType{data}
\title{
United Nations Table of Locations
}
\description{
United Nations table of locations, including regions, as available in 2010 and 2008.  It is stored in a TAB-separated ASCII files, \dQuote{WPP2010_LOCATIONS.txt} and \dQuote{WPP2010_LOCATIONS.txt}. }
\usage{
	data(WPP2010_LOCATIONS)
	
	data(WPP2008_LOCATIONS)
	}
\format{
  A data frame with one observations per country or region. It contains the following 8 variables:
  \describe{
    \item{\code{Index}}{Numeric vector enumerating the observations. Not used by the package.}
    \item{\code{Major.area..region..country.or.area}}{Name of country or region (following ISO 3166 official short names in English - see \cr \url{http://www.iso.org/iso/country_codes/iso_3166_code_lists/english_country_names_and_code_elements.htm} and United Nations Multilingual Terminology Database - see \url{http://unterm.un.org}).}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{reg_code}}{Code of the regions.}
    \item{\code{reg_name}}{Name of the regions.}
    \item{\code{area_code}}{Area code.}
    \item{\code{area_name}}{Area names, such as \code{Africa}, \code{Asia}, \code{Europe} \code{Latin America and the Caribbean}, \code{Northern America}, \code{Oceania}, \code{World}.}
    \item{\code{location_type}}{Code giving the type of the observation (0=World, 2=Major Area, 3=Region, 4=Country/Area, 5=Development group, 12=Special groupings). Not used by the package.}
    \item{\code{include_code}}{Specific code. Only entries for which \code{include_code=2} are included in MCMC simulations (i.e., estimation of the model parameters). Entries for which \code{include_code} is 1 or 2 are included in the prediction. Values of 0 and 3 are ignored.}
  }
}
\note{
The dataset differs from the one in the \pkg{bayesTFR} package only in the values of the \code{include_code} column. Here, countries affected by AIDS are not included in the estimation, i.e. the \code{include_code} is set to 3.  

}
\source{
Data provided by the United Nations Population Division
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(WPP2010_LOCATIONS, package='bayesLife')
}
\keyword{datasets}
