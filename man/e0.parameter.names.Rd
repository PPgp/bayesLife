\name{e0.parameter.names}
\alias{e0.parameter.names}
\alias{e0.parameter.names.cs}
\alias{e0.parameter.names.extended}
\alias{e0.parameter.names.cs.extended}

\title{
Accessing Parameter Names
}
\description{
Functions for accessing names of the MCMC parameters, either country-independent or country-specific.
}
\usage{
e0.parameter.names(...)
e0.parameter.names.cs()
e0.parameter.names.extended(...)
e0.parameter.names.cs.extended(country.code = NULL)
}
\arguments{
\item{country.code}{Country code. If it is given, the country-specific parameter names contain the suffix \sQuote{_c\eqn{x}} where \eqn{x} is the \code{country.code}.}
\item{...}{Here only the logical argument \code{is.hiv} is allowed. If \code{TRUE}, an additional parameter name (related to the HIV/AIDS model) is added.}
}

\value{
	\code{e0.parameter.names} returns names of the world parameters.\cr
	\code{e0.parameter.names.cs} returns names of the country-specific parameters.\cr
	\code{e0.parameter.names.extended} returns names of all world parameters in their extended format. I.e. parameters \sQuote{Triangle} and \sQuote{lambda} have the suffix \sQuote{_1}, \sQuote{_2}, \sQuote{_3}, and \sQuote{_4}. \cr
	\code{e0.parameter.names.cs.extended} returns names of all country-specific parameters in their extended format. I.e. parameters \sQuote{Triangle.c} and \sQuote{lambda.c} are in their extended format with the suffix \sQuote{_1}, \sQuote{_2} and \sQuote{_3}.
}

\author{
Hana Sevcikova
}

\examples{
e0.parameter.names()
e0.parameter.names.extended()
e0.parameter.names.cs()
e0.parameter.names.cs.extended()
}

\keyword{ manip }

