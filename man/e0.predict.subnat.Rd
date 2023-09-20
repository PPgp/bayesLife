\name{e0.predict.subnat}
\alias{e0.predict.subnat}
\alias{e0.jmale.predict.subnat}
\alias{subnat.gap.estimates}

\title{
Generating Posterior Trajectories of Subnational Life Expectancy at Birth
}
\description{
Generates posterior trajectories of the life expectancy at birth (e0) for subregions of given countries, for female and male. 
}
\usage{
e0.predict.subnat(countries, my.e0.file, 
    sim.dir = file.path(getwd(), "bayesLife.output"), 
    method = c("ar1", "shift", "scale"),
    predict.jmale = FALSE, my.e0M.file = NULL,
    end.year = 2100, start.year = NULL, output.dir = NULL, 
    annual = NULL, nr.traj = NULL, seed = NULL, 
    ar.pars = NULL, save.as.ascii = 0, verbose = TRUE, 
    jmale.estimates = NULL, \dots)
                
e0.jmale.predict.subnat(e0.pred, estimates = NULL, 	
    gap.lim = c(0,18), max.e0.eq1.pred = 86, my.e0.file = NULL, 
    save.as.ascii = 0, verbose = TRUE)
    
subnat.gap.estimates(annual = FALSE)
}

\arguments{
  \item{countries}{Vector of numerical country codes or country names.}
  \item{my.e0.file}{Tab-separated ASCII file containing the subnational e0 data. In \code{e0.predict.subnat}, if female and male e0 are projected jointly (i.e. \code{predict.jmale} is \code{TRUE}), this file should contain the female e0. In \code{e0.jmale.predict.subnat}, this file should contain the male e0. See Details for more information on its format.}
  \item{sim.dir}{Simulation directory with the national projections generated using \code{\link{e0.predict}}.}
  \item{method}{Method to use for the projections, see the reference paper.}
  \item{predict.jmale}{Logical determining if male projections should be generated as well. If \code{TRUE}, the argument \code{my.e0M.file} must be given.}
  \item{my.e0M.file}{Tab-separated ASCII file containing the subnational male e0 data.}
 \item{end.year}{End year of the projections.}
 \item{start.year}{Start year of the projections. By default, projections start at the same time point as the national projections.} 
   \item{output.dir}{Directory into which the resulting prediction objects and the trajectories are stored. See below for details.}
   \item{annual}{Logical indicating if the subnational projection should be on an annual scale or a 5-year scale. By default, 
the scale is matched to the national simulation given by \code{sim.dir}. If given, the scale must match to the scale of the subnational data provided in \code{my.e0.file}. 
If the subnational and national scales are not the same, 
the national trajectories are either interpolated (if \code{annual = TRUE} and the national simulation is not annual) or averaged to 5-year values 
(if \code{annual = FALSE} and the national simulation is annual).}
\item{nr.traj}{Number of trajectories to be generated. If \code{NULL}, the number of trajectories in the national projections is used.}
\item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible projections.}
\item{ar.pars}{Named vector containing the parameter estimates for the AR(1) method (i.e. if \code{method = "ar1"}, default). If given, it must have elements called \code{rho}, \code{U}, \code{a} and \code{b}. See the reference paper for details on the estimation. By default for a 5-year simulation, \code{c(rho = 0.95, U = 82.5, a = 0.0482, b = -0.0154)} is used. 
For an annual simulation these default parameters are scaled to \code{c(rho = 0.9898, U = 82.5, a = 0.01, b = -0.0032)}, see details below.}
\item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ASCII file, or \dQuote{all} in which case all trajectories are converted. By default no conversion is performed.}
\item{verbose}{Logical switching log messages on and off.}
\item{jmale.estimates, estimates}{List with estimates for the female-male gap model. The default values are retrieved using the function \code{subnat.gap.estimates()}.}
\item{\dots}{Additional arguments passed to \code{e0.jmale.predict.subnat}, which are \code{gap.lim} and \code{max.e0.eq1.pred}.}
\item{e0.pred}{Object of class \code{\link{bayesLife.prediction}}. It should be one element of the list returned by \code{e0.predict.subnat} corresponding to one country.}
\item{gap.lim, max.e0.eq1.pred}{The same meaning as in \code{\link{e0.jmale.predict}}.}
}

\details{
The \code{e0.predict.subnat} function implements the methodology described in Sevcikova and Raftery (2021). Given a set of national bayesLife projections, it applies one of the methods (AR(1), Shift or Scale) to each national trajectory and each subregion of given countries which yields subnational e0 projections. 

The file on subnational data passed into \code{my.e0.file} and \code{my.e0M.file} has to have a column \dQuote{country_code} with numerical values corresponding to countries given in the argument \code{countries}, and column \dQuote{reg_code} giving the numerical identifier of each subregion. Column \dQuote{name} should be used for subregion name, and column \dQuote{country_name} for country name. An optional column \dQuote{include_code} can be used to eliminate entries from processing. Entries with values of 1 or 2 will be included, all others will be ignored. Column \dQuote{last.observed} can be used to define which time period contains the last observed data point (given as integer, e.g. year in the middle of the time period). Remaining columns define the time periods, e.g. \dQuote{2000-2005}, \dQuote{2005-2010} for a 5-year simulation, or \dQuote{2020}, \dQuote{2021} for an annual simulation. The package contains an example of such dataset, see Example below.

The default AR(1) parameters for the \dQuote{ar1} method were designed for a 5-year simulation, see Sevcikova & Raftery (2021) for more detail. These are \eqn{\rho = 0.95, U = 82.5, a = 0.0482, b = -0.0154}. If an annual AR(1) process is desired, we use the following conversion for the autoregressive parameter \eqn{\rho} and the \eqn{a} and \eqn{b} parameters: 
\eqn{\rho^* = \rho^{(1/5)}, a^* = a * \sqrt{((1 - \rho^{(2/5)})/(1 - \rho^2))}, b^* = b * \sqrt{((1 - \rho^{(2/5)})/(1 - \rho^2))}}. The \eqn{U} parameter stays the same for both processes. Thus, the annual parameters are  \code{c(rho = 0.9898, U = 82.5, a = 0.01, b = -0.0032)}. Note that if the \code{ar.pars} argument is specified by the user, it is assumed that the parameters have been scaled appropriately and thus, no conversion takes place.

Argument \code{output.dir} gives a location on disk where results of the function should be stored. If it is \code{NULL} (default),  results are stored in the same directory as the national projections. In both cases a subdirectory called \dQuote{subnat_\code{method}} is created in which each country has its own subfolder with the country code in its name. Each such subfolder contains the same type of outputs as in the national case generated using \code{\link{e0.predict}}, most importantly a directory \dQuote{predictions} with trajectories for each region.

If the argument \code{predict.jmale} is \code{TRUE}, the \code{e0.predict.subnat} invokes  the \code{e0.jmale.predict.subnat} function for each country. However, one can call the \code{e0.jmale.predict.subnat} function explicitly. It applies the female-male gap model to regions of one country. See \code{\link{e0.jmale.predict}} for more detail on the model. The default  covariates of the model are not estimated on the fly. They were estimated externally using subnational data for about 30 countries and can be viewed using \code{subnat.gap.estimates()}, either for estimates derived from 5-year data (default) or annual data (\code{annual = TRUE}).
}

\value{
Function \code{e0.predict.subnat} returns a list of objects of class \code{\link{bayesLife.prediction}}. The name of each element includes its country code. Not all elements of the class \code{\link{bayesLife.prediction}} are available. For example, no \code{mcmc.set} is attached to these objects. Thus, not all functions that work with \code{\link{bayesLife.prediction}} can be applied to these results.

Function \code{e0.jmale.predict.subnat} returns an object of class \code{\link{bayesLife.prediction}} which updates the input \code{e0.pred} object by adding a new component called \code{joint.male}. This component is also an object of class \code{\link{bayesLife.prediction}} and it contains results of the male projections. 
}

\note{
Even though the resulting object contains subnational results, the names  of its elements are the same as in the national case. This allows to apply the same functions on both objects (subnational and national). However, it means that sometimes the meaning of the elements or function arguments does not match the subnational context. For example, various functions expect the argument \code{country}. When a subnational object is passed to such a function, \code{country} means a subregion. 
}

\references{
H. Sevcikova, A. E. Raftery (2021). Probabilistic Projection of Subnational Life Expectancy. Journal of Official Statistics, Vol. 37, no. 3, 591-610.
}

\author{
Hana Sevcikova
}


\seealso{
\code{\link{get.rege0.prediction}}, \code{\link{e0.predict}}, \code{\link{e0.jmale.predict}}
}

\examples{
# View the example data
my.sube0.file <- file.path(find.package("bayesLife"), 'extdata', 'subnational_e0_template.txt')
sube0 <- read.delim(my.sube0.file, check.names = FALSE)
head(sube0)

# Directory with national projections (contains 30 trajectories for each country)
nat.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")

# Subnational projections for Australia and Canada, 
# including the joint female-male gap model
subnat.dir <- tempfile()
preds <- e0.predict.subnat(c(36, 124), my.e0.file = my.sube0.file,
    sim.dir = nat.dir, output.dir = subnat.dir, start.year = 2018)
    
names(preds)
get.countries.table(preds[["36"]])
summary(preds[["36"]], "Queensland")
e0.trajectories.plot(preds[["36"]], "Queensland")

# plot subnational and national e0 in one plot
nat.pred <- get.e0.prediction(nat.dir)
e0.trajectories.plot(preds[["36"]], 4, pi = 80)
e0.trajectories.plot(nat.pred, "Australia", add = TRUE, 
      col = rep("darkgreen", 5), nr.traj = 0, show.legend = FALSE)
legend("top", c("regional e0", "national e0"), col = c("red", "darkgreen"), 
  lty = 1, bty = 'n')

# Add male projection to Canada, 
# using (wrongly) female data only for demonstration
predCan <- e0.jmale.predict.subnat(preds[["124"]], my.e0.file = my.sube0.file)

# retrieve male prediction object
predCanMale <- get.rege0.prediction(subnat.dir, 124, joint.male = TRUE)
# the same works using
predCanMale <- get.e0.jmale.prediction(predCan)

# Retrieve female and male trajectories
trajsF.Alberta <- get.e0.trajectories(predCan, "Alberta")
trajsM.Alberta <- get.e0.trajectories(predCanMale, "Alberta")
# summary of differences
summary(t(trajsF.Alberta - trajsM.Alberta))

# cleanup
unlink(subnat.dir)

# See more examples in ?get.rege0.prediction
}

\keyword{models}
\keyword{ts}

