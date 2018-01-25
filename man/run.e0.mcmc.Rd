\name{run.e0.mcmc}
\alias{run.e0.mcmc}
\alias{continue.e0.mcmc}
\alias{bayesLife.mcmc.set}

\title{
Running Markov Chain Monte Carlo for Parameters of Life Expectancy
}
\description{
Runs (or continues running) MCMCs for simulating the life expectancy for all countries of the world, using a Bayesian hierarchical model.
}
\usage{
run.e0.mcmc(sex = c("Female", "Male"), nr.chains = 3, iter = 160000, 
    output.dir = file.path(getwd(), "bayesLife.output"), 
    thin = 10, replace.output = FALSE, 
    start.year = 1873, present.year = 2015, wpp.year = 2017, hiv.model = FALSE,
    my.e0.file = NULL, my.locations.file = NULL, buffer.size = 100, 
    a = c(13.215, 41.070, 9.235, 17.605, 2.84, 0.385), 
    delta = c(3.844, 4.035, 11.538, 5.639, 0.901, 0.4), 
    tau = c(15.5976503, 23.650006, 14.5056919, 
            14.718598, 3.4514285, 0.5667531), 
    Triangle.ini = list(NULL, NULL, NULL, NULL), k.ini = NULL, 
    z.ini = NULL, lambda.ini = list(NULL, NULL, NULL, NULL), 
    lambda.k.ini = NULL, lambda.z.ini = NULL, omega.ini = NULL, 
    betanonART.ini = -0.77,
    Triangle.ini.low = c(10, 30, 0.1, 10), Triangle.ini.up = c(30, 50, 10, 30), 
    k.ini.low = 3, k.ini.up = 5, z.ini.low = 1e-04, z.ini.up = 0.653, 
    lambda.ini.low = c(0.01, 0.01, 0.01, 0.01), 
    lambda.ini.up = c(0.1, 0.1, 0.1, 0.1), 
    lambda.k.ini.low = 0.3, lambda.k.ini.up = 1, 
    lambda.z.ini.low = 1, lambda.z.ini.up = 40, 
    omega.ini.low = 0.1, omega.ini.up = 5, 
    Triangle.prior.low=c(0, 0, -20, 0), Triangle.prior.up=c(100, 100, 50, 100),
    k.prior.low = 0, k.prior.up = 10, z.prior.low = 0, z.prior.up = 0.653,
    Triangle.c.ini.norm = list(round(
    	Triangle.ini.low + (Triangle.ini.up - Triangle.ini.low)/2), c(2, 2, 2, 2)), 
    k.c.ini.norm = c(round(k.ini.low + (k.ini.up - k.ini.low)/2), 2), 
    z.c.ini.norm = c(round(z.ini.low + (z.ini.up - z.ini.low)/2, 2), 0.2), 
    Triangle.c.prior.low = c(0, 0, -20, 0), Triangle.c.prior.up = c(100, 100, 50, 100),
    k.c.prior.low = 0, k.c.prior.up = 10, z.c.prior.low = 0, z.c.prior.up = 0.653,
    country.overwrites = NULL, nu = 4, dl.p1 = 9, dl.p2 = 9, 
    sumTriangle.lim = c(30, 86), constant.variance = FALSE, 
    outliers=c(-5, 10), hiv.outliers = c(-15, 15), seed = NULL, 
    parallel = FALSE, nr.nodes = nr.chains, compression.type = 'None',
    auto.conf = list(max.loops=5, iter=160000, iter.incr=20000, 
        nr.chains=3, thin=225, burnin=10000),
    verbose = FALSE, verbose.iter = 100, \dots)
    
continue.e0.mcmc(iter, chain.ids = NULL, 
    output.dir = file.path(getwd(), "bayesLife.output"), 
    parallel = FALSE, nr.nodes = NULL, auto.conf = NULL, 
    verbose = FALSE, verbose.iter = 10, \dots)
}

\arguments{
  \item{sex}{Sex for which to run the simulation.}
  \item{nr.chains}{Number of MCMC chains to run.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} in which case the function runs for the number of iterations given in the \code{auto.conf} list (see below), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetition exceeded \code{auto.conf$max.loops}.}
  \item{output.dir}{Directory which the simulation output should be written into.}
  \item{thin}{Thinning interval between consecutive observations to be stored on disk.}
  \item{replace.output}{If \code{TRUE}, existing outputs in \code{output.dir} will be replaced by results of this simulation.}
  \item{start.year}{Start year for using historical data.}
  \item{present.year}{End year for using historical data.}
  \item{wpp.year}{Year for which WPP data is used. The functions loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses the \code{e0*} datasets.}
  \item{hiv.model}{Logical determining if HIV/AIDS countries should be included in the estimation and treated using the Godwin & Raftery (2017) model. If \code{FALSE}, HIV/AIDS countries are not included in the estimation.}
  \item{my.e0.file}{File name containing user-specified e0 time series for one or more countries. See Details below.}
  \item{my.locations.file}{File name containing user-specified locations. See Details below.}
  \item{buffer.size}{Buffer size (in number of [thinned] iterations) for keeping data in the memory. The smaller the \code{buffer.size} the more often will the process access the hard disk and thus, the slower the run. On the other hand, the smaller the \code{buffer.size} the less data will be lost in case of failure.}
  \item{a}{A vector of the \eqn{a_1}, \dots ,\eqn{a_6} parameters, which are the prior means of the world-level parameters 
  			(\eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}, \eqn{k}, \eqn{z}).}
  \item{delta}{A vector of the \eqn{\delta_1}{delta_1}, \dots ,\eqn{\delta_6}{delta_6} parameters, 
  	which are the prior standard deviations of the world-level parameters 
  			(\eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}, \eqn{k}, \eqn{z}).}
  \item{tau}{A vector of the \eqn{\tau_1}{tau_1}, \dots ,\eqn{\tau_6}{tau_6} parameters, which is the square root rate of the prior Gamma distribution of the world-level parameters 
  			(\eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}, \eqn{\lambda_k}{lambda_k}, \eqn{\lambda_z}{lambda_z}).}
  \item{Triangle.ini}{List (of length four) of initial values for \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}. Each list item should be of the length \code{nr.chains}. If a list item is \code{NULL}, the initial values are equally spaced between \code{Triangle.ini.low} and \code{Triangle.ini.up}. By default, if there is just one chain, the value is the middle point of the interval.}
  \item{k.ini, z.ini}{An array of length \code{nr.chains} of initial values for \eqn{k} (\eqn{z}). By default, the initial values are equally spaced between \code{k.ini.low} and \code{k.ini.up} (\code{z.ini.low} and \code{z.ini.up}). If \code{nr.chains=1} and \code{k.ini} (\code{z.ini})is \code{NULL}, the initial value is the middle point of the interval.}
  \item{lambda.ini}{List (of length four) of initial values for \eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}. Each list item should be of the length \code{nr.chains}. If a list item is \code{NULL}, the initial values are equally spaced between \code{lambda.ini.low} and \code{lambda.ini.up}. By default, if there is just one chain, the value is the middle point of the interval.}
  \item{lambda.k.ini, lambda.z.ini}{An array of length \code{nr.chains} of initial values for \eqn{\lambda_k}{lambda_k} (\eqn{\lambda_z}{lambda_z}). By default, the initial values are equally spaced between \code{lambda.k.ini.low} and \code{lambda.k.ini.up} (\code{lambda.z.ini.low} and \code{lambda.z.ini.up}). If \code{nr.chains=1} and \code{lambda.k.ini} (\code{lambda.z.ini}) is \code{NULL}, the initial value is the middle point of the interval.}
  \item{omega.ini}{An array of length \code{nr.chains} of initial values for \eqn{\omega}{omega}. By default, the initial values are equally spaced between \code{omega.ini.low} and \code{omega.ini.up}. If \code{nr.chains=1} and \code{omega.ini} is \code{NULL}, the initial value is the middle point of the interval.}
  \item{betanonART.ini}{Initital value for the regression coefficient of the HIV/AIDS model (1b). Only used if \code{hiv.model} is \code{TRUE}.}
  \item{Triangle.ini.low, Triangle.ini.up}{A vector of length four. They are the lower and upper bounds for initial values of \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}. An \eqn{i}-th item is only used if \code{Triangle.ini[[i]]} is \code{NULL}.}
  \item{k.ini.low, k.ini.up}{Single value giving the lower and upper bounds for initial values of \eqn{k}. It is only used if \code{k.ini} is \code{NULL}.}
  \item{z.ini.low, z.ini.up}{Single value giving the lower and upper bounds for initial values of \eqn{z}. It is only used if \code{z.ini} is   
  	\code{NULL}. Regarding defaults, see Note below.}
 \item{lambda.ini.low, lambda.ini.up}{A vector of length four. They are the lower and upper bounds for initial values of \eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}. An \eqn{i}-th item is only used if \code{lambda.ini[[i]]} is \code{NULL}.}
  \item{lambda.k.ini.low, lambda.k.ini.up}{Single value giving the lower and upper bounds for initial values of \eqn{\lambda_k}{lambda_k}. It is only used if \code{lambda.k.ini} is \code{NULL}.}
  \item{lambda.z.ini.low, lambda.z.ini.up}{Single value giving the lower and upper bounds for initial values of \eqn{\lambda_z}{lambda_z}. It is only used if \code{lambda.z.ini} is \code{NULL}.}
  \item{omega.ini.low, omega.ini.up}{Single value giving the lower and upper bounds for initial values of \eqn{\omega}{omega}. It is only used if \code{omega.ini} is \code{NULL}.}
  \item{Triangle.prior.low, Triangle.prior.up}{A vector of length four. They are the lower and upper bounds for the prior (truncated normal) distribution of \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}.}
    \item{k.prior.low, k.prior.up}{Single value giving the lower and upper bounds for the prior (truncated normal) distribution  of \eqn{k}. }
    \item{z.prior.low, z.prior.up}{Single value giving the lower and upper bounds for the prior (truncated normal) distribution  of \eqn{z}. Regarding defaults, see Note below.}
  \item{Triangle.c.ini.norm}{A list with two elements, each of them being a vector of size four. The first and second element in the list corresponds to the means and standard deviation, respectively, for the initial values of the country-specific parameters \eqn{\Delta_1^c}{Triangle_1_c}, \dots, \eqn{\Delta_4^c}{Triangle_4_c} which are drawn from a truncated normal distribution with bounds defined by \code{Triangle.c.prior.low} and \code{Triangle.c.prior.up}.}
  \item{k.c.ini.norm, z.c.ini.norm}{A vector of length two. The first and second element corresponds to the means and standard deviation, respectively, for the initial values of the country-specific parameters \eqn{k^c} (\eqn{z^c}) which are drawn from a normal distribution truncated between \code{k.c.prior.low} and \code{k.c.prior.up} (\code{z.c.prior.low} and \code{z.c.prior.up}).}
    \item{Triangle.c.prior.low, Triangle.c.prior.up}{A vector of length four. They are the lower and upper bounds for the prior (truncated normal) distribution of country-specific \eqn{\Delta_1^c}{Triangle_1_c}, \dots, \eqn{\Delta_4^c}{Triangle_4_c}.}
    \item{k.c.prior.low, k.c.prior.up}{Single value giving the lower and upper bounds for the prior (truncated normal) distribution  of country-specific \eqn{k^c}. }
    \item{z.c.prior.low, z.c.prior.up}{Single value giving the lower and upper bounds for the prior (truncated normal) distribution  of country-specific \eqn{z^c}. Regarding defaults, see Note below.}
  \item{country.overwrites}{This argument allows to overwrite some of the prior parameters for specific countries. If it is not \code{NULL} it should be a data frame with an obligatory column \sQuote{country_code}. Each row then corresponds to one country. Other columns can be \sQuote{k.c.prior.low}, \sQuote{k.c.prior.up}, \sQuote{z.c.prior.low}, \sQuote{z.c.prior.up}, \sQuote{Triangle_\eqn{x}.c.prior.low} and \sQuote{Triangle_\eqn{x}.c.prior.up} where \eqn{x} can be an integer from 1 to 4.}
  \item{nu}{The shape parameter of the Gamma distributions of all \eqn{\lambda}{lambda} parameters is \code{nu/2}.}
  \item{dl.p1, dl.p2}{Values of the parameters \eqn{p_1} and \eqn{p_2} of the double logistic function.}
  \item{sumTriangle.lim}{Lower and upper limits for the sum of the \eqn{\Delta_i}{Triangle_i} parameters. MCMC proposals that are outside of this limit are rejected. It is applied to both, the world parameters as well as the country specific parameters.}
  \item{constant.variance}{Logical indicating if the model should be estimated using constant variance.}
  \item{outliers}{Ranges for determining outliers in the historical data. If \code{outliers=c(x, y)} then any increase in life expectancy smaller than \code{x} or larger than \code{y} is considered as an outlier and removed from the estimation.}
  \item{hiv.outliers}{Ranges for outliers of the HIV/AIDS countries. Only used if \code{hiv.model} is \code{TRUE}.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible results.}
  \item{parallel}{Logical determining if the simulation should run multiple chains in parallel. If it is \code{TRUE}, the package \pkg{\link[snowFT]{snowFT}} is required.}
  \item{nr.nodes}{Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. By default it equals to the number of chains.}
  \item{compression.type}{One of \sQuote{None}, \sQuote{gz}, \sQuote{xz}, \sQuote{bz}, determining type of a compression of the MCMC files.}
    \item{auto.conf}{List containing a configuration for an \sQuote{automatic} run (see description of argument \code{iter}). Item \code{iter} gives the number of iterations in the first chunk of the MCMC simulation; item \code{iter.incr} gives the number of iterations in the following chunks; \code{nr.chains} gives the number of chains in all chunks of the MCMC simulation; items \code{thin} and \code{burnin} are used in the convergence diagnostics following each chunk; \code{max.loops} controls the maximum number of chunks. All items must be integer values. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{verbose.iter}{Integer determining how often (in number of iterations) log messages are outputted during the estimation.}
  \item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
  \item{chain.ids}{Array of chain identifiers that should be resumed. If it is \code{NULL}, all existing chains in \code{output.dir} are resumed.}
}
\details{
	The function \code{run.e0.mcmc} creates an object of class \code{\link{bayesLife.mcmc.meta}} and stores it in \code{output.dir}.  It launches \code{nr.chains} MCMCs, either sequentially or in parallel.  Parameter traces of each chain are stored as (possibly compressed) ASCII files in a subdirectory of \code{output.dir}, called \code{mc}\emph{x} where \emph{x} is the identifier of that chain.  There is one file per parameter, named after the parameter with the suffix \dQuote{.txt}, possibly followed by a compression suffix if \code{compression.type} is given.  Country-specific parameters have the suffix \code{_country}\emph{c} where \emph{c} is the country code.  In addition to the trace files, each \code{mc}\emph{x} directory contains the object \code{\link{bayesLife.mcmc}} in binary format.  All chain-specific files  are written into disk after the first, last and each \code{buffer.size}-th (thinned) iteration.
	
Using the function \code{continue.e0.mcmc} one can continue simulating an existing MCMCs by \code{iter} iterations for either all or selected chains.

The function loads observed data (further denoted as WPP dataset), depending on the specified sex,  from the \code{\link[wpp2017]{e0F}} (\code{\link[wpp2017]{e0M}}) and \code{\link[wpp2017]{e0F_supplemental}} (\code{\link[wpp2017]{e0M_supplemental}}) datasets in a \pkg{wpp}\eqn{x} package where \eqn{x} is the \code{wpp.year}. It is then merged with the  \code{\link{include}} dataset that corresponds to the same \code{wpp.year}. The argument \code{my.e0.file} can be used to overwrite those default data. Such a file can include a subset of countries contained in the WPP dataset, as well as a set of new countries. In the former case,
  	the function replaces the corresponding country data from the WPP dataset with values in this file. Only columns are replaced that match column names of the WPP dataset, and in addition, columns \sQuote{last.observed} and \sQuote{include_code} are used, if present. Countries are merged with WPP using the column \sQuote{country_code}. In addition, in order the countries to be included in the simulation, in both cases (whether they are included in the WPP dataset or not), they must be contained in the table of locations (\code{\link[wpp2017]{UNlocations}}). In addition, their corresponding \sQuote{include_code} must be set to 2. If the column \sQuote{include_code} is present in \code{my.e0.file}, its value overwrites the default include code, unless is -1.

The default UN table of locations mentioned above can be overwritten/extended by using a file passed as the \code{my.locations.file} argument. Such a file must have the same structure as the \code{\link[wpp2017]{UNlocations}} dataset. Entries in this file will overwrite corresponding entries in \code{\link[wpp2017]{UNlocations}} matched by the column \sQuote{country_code}. If there is no such entry in the default dataset, it will be appended. This option of appending new locations is especially useful in cases when \code{my.e0.file} contains new countries/regions that are not included in \code{\link[wpp2017]{UNlocations}}. In such a case, one must provide a \code{my.locations.file} with a definition of those countries/regions.


For simulation of the hyperparameters of the Bayesian hierarchical model, all countries are used that are included in the WPP dataset, possibly complemented by the \code{my.e0.file}, that have \code{include_code} equal to 2. The hyperparameters are used to simulate country-specific parameters, which is done for all countries with \code{include_code} equal 1 or 2.  The following values of \code{include_code} in \code{my.e0.file} are recognized: -1 (do not overwrite the default include code), 0 (ignore), 1 (include in prediction but not estimation), 2 (include in both, estimation and prediction). Thus, the set of countries included in the estimation and prediction can be fully specified by the user.

Optionally, \code{my.e0.file} can contain a column called \code{last.observed} containing the year of the last observation for each country. In such a case, the code would ignore any data after that time point. Furthermore, the function \code{\link{e0.predict}} fills in the missing values using the median of the BHM procedure (stored in \code{e0.matrix.reconstructed} of the \code{\link{bayesLife.prediction}} object). For \code{last.observed} values that are below a middle year of a time interval \eqn{[t_i, t_{i+1}]} (computed as \eqn{t_i+3}) the last valid data point is the time interval \eqn{[t_{i-1}, t_i]}, whereas for values larger equal a middle year, the data point in \eqn{[t_i, t_{i+1}]} is valid.

The package contains a dataset called \file{my_e0_template} (in \file{extdata} directory) which is a template for user-specified \code{my.e0.file}.
}
\note{
	Parameter \eqn{z} determines the asymptote in gains in life expectancy. The following text gives an explanation for the choice of  upper limits on \eqn{z}-related parameters: 
	
	The pace of improvement and the asymptotic limit in future gains in female life expectancy vary for each projected trajectory, but ultimately is informed and constrained by the finding that the rate of increase of maximum female life expectancy over the past 150 year has been highly linear (2a, 2b) (i.e., about 2.4 years per decade), albeit at slightly lower pace once the leading countries started to exceed 75 years of female life expectancy at birth in the 1960s (3) (about 2.26 years of gains per decade). By assuming that the asymptotic average rate of increase in life expectancy is nonnegative, life expectancy is assumed to continually increase (on average), and no limit is imposed to life expectancy in the foreseeable future. The increase in maximum female life span among countries with highest life expectancy and reliable data on very old age provide further guidance on future rate of progress which has also been increasingly linear at least since the 1970s (4a-4c) (about 1.25 years per decade for countries like Sweden and Norway), and is used to inform the asymptotic average rate of increase in female life expectancy used in the 2012 WPP Revision. To set the posterior median to an annual gain of 0.125 year (or 5-year gain of 0.625 in this context) the upper bound value of 0.653 is used for the world prior (\eqn{z}) and country-specific prior (\eqn{z_c}) as default values in the estimation of the double-logistic parameters.
}
\value{
	An object of class \code{bayesLife.mcmc.set} which is a list with two components:
	\item{meta}{An object of class \code{\link{bayesLife.mcmc.meta}}.}
	\item{mcmc.list}{A list of objects of class \code{\link{bayesLife.mcmc}}, one for each MCMC.}
}
\references{
(1a) Raftery, A.E., Chunn, J.L., Gerland, P., and Sevcikova, H. (2013). Bayesian probabilistic projections of life expectancy for all countries. Demography 50(3): 777–801. \url{https://dx.doi.org/10.1007/s13524-012-0193-x}

(1b) J. Godwin and A. E. Raftery (2017): Bayesian projection of life expectancy accounting for the HIV/AIDS epidemic. Demographic Research 37(48):1549–1610. \url{http://www.demographic-research.org/Volumes/Vol37/48}

(2a)  Oeppen J, and J.W. Vaupel (2002) Broken limits to life expectancy. Science 296:1029-1031. 

(2b) Vaupel, J.W. and K.G.V. Kistowski. 2005. Broken Limits to Life Expectancy. Ageing Horizons (3):6-13.

(3) Vallin, J., and F. Mesle (2009). The Segmented Trend Line of Highest Life Expectancies. Population and Development Review, 35(1), 159-187. doi:10.1111/j.1728-4457.2009.00264.x

(4a) Wilmoth, J. R., L. J. Deegan, H. Lundstrom, and S. Horiuchi (2000). Increase of maximum life-span in Sweden, 1861-1999. Science, 289(5488), 2366-2368. 

(4b) Wilmoth, J. R. and J-M. Robine. (2003). The world trend in maximum life span, in: J. R. Carey and S. Tuljapurkar (eds.), Life Span: Evolutionary, Ecological, and Demographic Perspectives, supplement to vol. 29, Population and Development Review, pp. 239-257. 

(4c) Wilmoth, J. R. and N. Ouellette (2012). Maximum human lifespan: Will the records be unbroken?, Paper presented at the European Population Conference, Stockholm, Sweden, 13-16 June.
}
\author{
Hana Sevcikova, Patrick Gerland contributed to the documentation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get.e0.mcmc}}, \code{\link{summary.bayesLife.mcmc.set}}.
}
\examples{
\dontrun{
m <- run.e0.mcmc(nr.chains=1, iter=5, thin=1, verbose=TRUE)
summary(m)
m <- continue.e0.mcmc(iter=5, verbose=TRUE)
summary(m)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
