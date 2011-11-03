\name{UN2008e0}
\alias{UN2008e0F}
\alias{UN2008e0M}
\alias{UN2008e0}
\alias{UN2010e0F}
\alias{UN2010e0M}
\alias{UN2010e0}
\alias{my_e0_template}
\docType{data}
\title{
United Nations Time Series of the Life Expectancy
}
\description{
Data sets containing the United Nations time series of the life expectancy for all countries of the world as available in 2010 and 2008. These are TAB-separated ASCII files, \dQuote{UN2010e0F.txt}, {UN2010e0M.txt}, {UN2008e0F.txt}, and \dQuote{UN2008e0M.txt}, corresponding to female and male life expectancy, respectively, in the particular years. \code{my_e0_template} is a template for user-specific data.
}
\usage{
data(UN2010e0F)

data(UN2010e0M)

data(UN2008e0F)

data(UN2008e0M)

data(my_e0_template)
}

\format{
  A data frame with 229 observations, one per country or region. It contains the following variables:
  \describe{
    \item{\code{Index}}{Numeric vector enumerating the observations. Not used by the package.}
    \item{\code{country}}{Name of country or region (following ISO 3166 official short names in English - see \url{http://www.iso.org/iso/country_codes/iso_3166_code_lists/english_country_names_and_code_elements.htm} and United Nations Multilingual Terminology Database - see \url{http://unterm.un.org}).}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{1950-1955}}{Life expectancy (e0) in time interval from 1950 to 1955.}
    \item{\code{1955-1960}}{e0 in time interval from 1955 to 1960.}
    \item{\code{1960-1965}}{e0 in time interval from 1960 to 1965.}
    \item{\code{1965-1970}}{e0 in the given time interval }
    \item{\code{1970-1975}}{e0 in the given time interval  }
    \item{\code{1975-1980}}{e0 in the given time interval  }
    \item{\code{1980-1985}}{e0 in the given time interval  }
    \item{\code{1985-1990}}{e0 in the given time interval  }
    \item{\code{1990-1995}}{e0 in the given time interval  }
    \item{\code{1995-2000}}{e0 in the given time interval  }
    \item{\code{2000-2005}}{e0 in the given time interval  }
    \item{\code{2005-2010}}{e0 in the given time interval  }
  }
}
\details{
All entries in these data sets must be included in the data set of locations, \code{\link{WPP2010_LOCATIONS}} and \code{\link{WPP2008_LOCATIONS}}, respectively. They are matched using the column \code{country_code}.

Regarding the time columns, only those are necessary that correspond to an earlier or equal time point than the value of the argument \code{present.year} passed to the function \code{\link{run.e0.mcmc}}. Thus, if \code{present.year=2003}, the column \code{2005-2010} is not used.

Optionally, the data sets can contain a column called \code{last.observed} containing the year of the last observation for each country. In such a case, the code would ignore any data after that time point. Furthermore, the function \code{\link{e0.predict}} fills in the missing values using the median of the BHM procedure (stored in \code{e0.matrix.reconstructed} of the \code{\link{bayesLife.prediction}} object). For \code{last.observed} values that are below a middle year of a time interval \eqn{[t_i, t_{i+1}]} (computed as \eqn{t_i+3}) the last valid data point is the time interval \eqn{[t_{i-1}, t_i]}, whereas for values larger equal a middle year, the data point in \eqn{[t_i, t_{i+1}]} is valid.

The package contains also a data set called \sQuote{my_e0_template} which is a template for user-specified e0 time series. It has the same structure, except that it can contain any subset of the set of the columns above. The only required column is \code{country_code} (see description of the argument \code{my.e0.file} in \code{\link{run.e0.mcmc}}).
}
\source{
Data provided by the United Nations Population Division; see also \url{http://esa.un.org/unpd/wpp}.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(UN2008e0M)
head(UN2008e0M)
}
\keyword{datasets}
