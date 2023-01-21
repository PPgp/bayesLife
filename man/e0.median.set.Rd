\name{e0.median.set}
\alias{e0.median.set}
\alias{e0.median.shift}
\alias{e0.median.reset}
\alias{e0.median.adjust.jmale}
\alias{e0.shift.prediction.to.wpp}

\title{
Editing the Projection Medians
}
\description{
These functions are to be used by expert analysts. They allow to change the projection medians either to specific values, including the WPP values, or shift the medians by a given constant or a factor.
}
\usage{
e0.median.set(sim.dir, country, values, years = NULL, joint.male = FALSE)

e0.median.shift(sim.dir, country, reset = FALSE, shift = 0, 
    from = NULL, to = NULL, joint.male = FALSE)
    
e0.median.adjust.jmale(sim.dir, countries, factors = c(1.2, 1.1))

e0.median.reset(sim.dir, countries = NULL, joint.male = FALSE)

e0.shift.prediction.to.wpp(sim.dir, joint.male = FALSE, ...)
}
\arguments{
  \item{sim.dir}{Directory containing the prediction object.}
  \item{country}{Name or numerical code of a country.}
  \item{countries}{Vector of country names or codes. In the \code{e0.median.reset} function, if this argument is \code{NULL} (default), the reset is done for all countries.}
  \item{values}{Array of the new median values.}
  \item{years}{Numeric vector giving years which \code{values} correspond to. Ideally it should be of the same length as \code{values}. If it is \code{NULL},
  	\code{values} are set starting from the first prediction period. If \code{values} correspond to consecutive years, only the first year might be given here. A year \eqn{t} represents a prediction period \eqn{[t_i, t_{i+1}]} if \eqn{t_i < t \leq t_{i+1}}.}
  \item{joint.male}{Logical. If \code{TRUE}, the function is applied to a male prediction that was generated using the joint female-male model implemented in the function \code{\link{e0.jmale.predict}}.}
  \item{reset}{Logical. If \code{TRUE} medians in a range of \code{from} and \code{to} are reset to their original values.}
  \item{shift}{Constant by which the medians should be offset. It is not used if \code{reset} is \code{TRUE}.}
  \item{from}{Year from which the offset/reset should start. By default, it starts at the first prediction period.}
  \item{to}{Year until which the offset/reset should be done. By default, it is set to the last prediction period.}
  \item{factors}{It should be a vector where each element corresponds to one time period. The adjustment of male medians is done as \code{e0m(t) = e0f(t) - gap(t)*factor(t)}.}
  \item{\dots}{Additional arguments passed to the underlying adjustment function. It can be \code{verbose} to show/hide the progress of the adjustment and \code{wpp.year} to adjust it to if it differs from the wpp year of the simulation.}
}
\details{
The function \code{e0.median.set} can be used to set the medians of the given country to specific values. Function \code{e0.median.shift} can be used to offset the medians by a specific constant, or to reset the medians to their original BHM values. Function \code{e0.median.adjust.jmale} adjusts male medians using factors that can expand or shrink the female-male gap. 

Function\code{e0.shift.prediction.to.wpp} shifts the projected medians so that they correspond to the values found in the \code{e0Fproj} (\code{joint.male = FALSE}) or \code{e0Mproj} (\code{joint.male = TRUE}) datasets of the \pkg{wpp} package that either corresponds to the package used for the simulation itself or is given by the \code{wpp.year} argument. If using \pkg{wpp2022}, the dataset name is automatically adjusted depending if it is an annual or a 5-year simulation.

Function \code{e0.median.reset} resets medians of the given countries to the original values. By default it deletes adjustments for all countries.

In all cases, if a median is modified, the corresponding offset is stored in the prediction object (element \code{median.shift}). All functions write the updated prediction object back to disk. All functions in the package that use trajectories and trajectory statistics use the \code{median.shift} values to offset the results correspondingly, i.e. trajectories are shifted the same way as the medians.
}
\value{
All functions return an updated object of class \code{\link{bayesLife.prediction}}.
}

\author{Hana Sevcikova}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

