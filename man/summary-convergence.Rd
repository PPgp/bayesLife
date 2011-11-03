\name{summary.bayesLife.convergence}
\alias{summary.bayesLife.convergence}

\title{
Summary of a Life Expectancy Convergence Object
}
\description{
Summary of an object of class \code{\link{bayesLife.convergence}} created using the \code{\link{e0.diagnose}} function. It gives an overview about parameters that did not converge.
}
\usage{
\method{summary}{bayesLife.convergence}(object, expand = FALSE, \dots)
}

\arguments{
  \item{object}{Object of class \code{\link{bayesLife.convergence}}.}
  \item{expand}{By default, the function does not show parameters for each country for which there was no convergence, if the status is \sQuote{red}. This argument can switch that option on.}
  \item{\dots}{Not used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
\seealso{
\code{\link{e0.diagnose}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }

