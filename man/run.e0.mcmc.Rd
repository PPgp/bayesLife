\name{run.e0.mcmc}
\alias{run.e0.mcmc}
\alias{continue.e0.mcmc}
\alias{bayesLife.mcmc.set}

\title{
Running Bayesian Hierarchical Model for Life Expectancy via Markov Chain Monte Carlo
}
\description{
Runs (or continues running) MCMCs for simulating the life expectancy for all countries of the world, using a Bayesian hierarchical model.
}
\usage{
run.e0.mcmc(sex = c("Female", "Male"), nr.chains = 3, iter = 160000, 
    output.dir = file.path(getwd(), "bayesLife.output"), 
    thin = 10, replace.output = FALSE, 
    start.year = 1873, present.year = 2020, wpp.year = 2019, 
    my.e0.file = NULL, my.locations.file = NULL, constant.variance = FALSE, 
    seed = NULL, parallel = FALSE, nr.nodes = nr.chains, compression.type = 'None',
    verbose = FALSE, verbose.iter = 100, mcmc.options = NULL, \dots)
    
continue.e0.mcmc(iter, chain.ids = NULL, 
    output.dir = file.path(getwd(), "bayesLife.output"), 
    parallel = FALSE, nr.nodes = NULL, auto.conf = NULL, 
    verbose = FALSE, verbose.iter = 10, \dots)
}

\arguments{
  \item{sex}{Sex for which to run the simulation.}
  \item{nr.chains}{Number of MCMC chains to run.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} for an automatic assessment of the convergence. In such a case, the function runs for the number of iterations given in the global option \code{auto.conf} list (see \code{\link{e0mcmc.options}}), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetition exceeded \code{auto.conf$max.loops}.}
  \item{output.dir}{Directory which the simulation output should be written into.}
  \item{thin}{Thinning interval between consecutive observations to be stored on disk.}
  \item{replace.output}{If \code{TRUE}, existing outputs in \code{output.dir} will be replaced by results of this simulation.}
  \item{start.year}{Start year for using historical data.}
  \item{present.year}{End year for using historical data.}
  \item{wpp.year}{Year for which WPP data is used. The functions loads a package called \pkg{wpp}\eqn{x} where \eqn{x} is the \code{wpp.year} and uses the \code{\link[wpp2019]{e0}*} datasets.}
  \item{my.e0.file}{File name containing user-specified e0 time series for one or more countries. See Details below.}
  \item{my.locations.file}{File name containing user-specified locations. See Details below.}
  \item{constant.variance}{Logical indicating if the model should be estimated using constant variance. It should only be used if the standard deviation lowess is to be analysed, see \code{\link{compute.loess}}.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible results.}
  \item{parallel}{Logical determining if the simulation should run multiple chains in parallel. If it is \code{TRUE}, the package \pkg{\link[snowFT]{snowFT}} is required.}
  \item{nr.nodes}{Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. By default it equals to the number of chains.}
  \item{compression.type}{One of \sQuote{None}, \sQuote{gz}, \sQuote{xz}, \sQuote{bz}, determining type of a compression of the MCMC files.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{verbose.iter}{Integer determining how often (in number of iterations) log messages are outputted during the estimation.}
  \item{mcmc.options}{List of options that overwrites global MCMC options as defined in \code{\link{e0mcmc.options}}. Type \code{e0mcmc.options()} to view default values.}
  \item{auto.conf}{In \code{continue.e0.mcmc}, one can overwrite the global \code{auto.conf} option, see \code{\link{e0mcmc.options}} for its definition. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
  \item{chain.ids}{Array of chain identifiers that should be resumed. If it is \code{NULL}, all existing chains in \code{output.dir} are resumed.}
}
\details{
The function \code{run.e0.mcmc} uses a set of global options (for priors, initial values etc.), possibly modified by the \code{mcmc.options} argument. One can also modify these options using \code{\link{e0mcmc.options}}. Call \code{e0mcmc.options()} for the full set of options. Function \code{continue.e0.mcmc} inherits its set of options from the corresponding \code{run.e0.mcmc} call.

The function \code{run.e0.mcmc} creates an object of class \code{\link{bayesLife.mcmc.meta}} and stores it in \code{output.dir}.  It launches \code{nr.chains} MCMCs, either sequentially or in parallel.  Parameter traces of each chain are stored as (possibly compressed) ASCII files in a subdirectory of \code{output.dir}, called \code{mc}\emph{x} where \emph{x} is the identifier of that chain.  There is one file per parameter, named after the parameter with the suffix \dQuote{.txt}, possibly followed by a compression suffix if \code{compression.type} is given.  Country-specific parameters have the suffix \code{_country}\emph{c} where \emph{c} is the country code.  In addition to the trace files, each \code{mc}\emph{x} directory contains the object \code{\link{bayesLife.mcmc}} in binary format.  All chain-specific files  are written into disk after the first, last and each \eqn{i}-th (thinned) iteration, where \eqn{i} is given by the global option \code{buffer.size}.
	
Using the function \code{continue.e0.mcmc} one can continue simulating an existing MCMCs by \code{iter} iterations for either all or selected chains. The global options used for generating the existing MCMCs will be used. Only the \code{auto.conf} option can be overwritten by passing the new value as an argument.

The function loads observed data (further denoted as WPP dataset), depending on the specified sex,  from the \code{\link[wpp2019]{e0F}} (\code{\link[wpp2019]{e0M}}) and \code{\link[wpp2019]{e0F_supplemental}} (\code{\link[wpp2019]{e0M_supplemental}}) datasets in a \pkg{wpp}\eqn{x} package where \eqn{x} is the \code{wpp.year}. It is then merged with the  \code{\link{include}} dataset that corresponds to the same \code{wpp.year}. The argument \code{my.e0.file} can be used to overwrite those default data. Such a file can include a subset of countries contained in the WPP dataset, as well as a set of new countries. In the former case,
  	the function replaces the corresponding country data from the WPP dataset with values in this file. Only columns are replaced that match column names of the WPP dataset, and in addition, columns \sQuote{last.observed} and \sQuote{include_code} are used, if present. Countries are merged with WPP using the column \sQuote{country_code}. In addition, in order the countries to be included in the simulation, in both cases (whether they are included in the WPP dataset or not), they must be contained in the table of locations (\code{\link[wpp2019]{UNlocations}}). In addition, their corresponding \sQuote{include_code} must be set to 2. If the column \sQuote{include_code} is present in \code{my.e0.file}, its value overwrites the default include code, unless is -1.

The default UN table of locations mentioned above can be overwritten/extended by using a file passed as the \code{my.locations.file} argument. Such a file must have the same structure as the \code{\link[wpp2019]{UNlocations}} dataset. Entries in this file will overwrite corresponding entries in \code{\link[wpp2019]{UNlocations}} matched by the column \sQuote{country_code}. If there is no such entry in the default dataset, it will be appended. This option of appending new locations is especially useful in cases when \code{my.e0.file} contains new countries/regions that are not included in \code{\link[wpp2019]{UNlocations}}. In such a case, one must provide a \code{my.locations.file} with a definition of those countries/regions.


For simulation of the hyperparameters of the Bayesian hierarchical model, all countries are used that are included in the WPP dataset, possibly complemented by the \code{my.e0.file}, that have \code{include_code} equal to 2. The hyperparameters are used to simulate country-specific parameters, which is done for all countries with \code{include_code} equal 1 or 2.  The following values of \code{include_code} in \code{my.e0.file} are recognized: -1 (do not overwrite the default include code), 0 (ignore), 1 (include in prediction but not estimation), 2 (include in both, estimation and prediction). Thus, the set of countries included in the estimation and prediction can be fully specified by the user.

Optionally, \code{my.e0.file} can contain a column called \code{last.observed} containing the year of the last observation for each country. In such a case, the code would ignore any data after that time point. Furthermore, the function \code{\link{e0.predict}} fills in the missing values using the median of the BHM procedure (stored in \code{e0.matrix.reconstructed} of the \code{\link{bayesLife.prediction}} object). For \code{last.observed} values that are below a middle year of a time interval \eqn{[t_i, t_{i+1}]} (computed as \eqn{t_i+3}) the last valid data point is the time interval \eqn{[t_{i-1}, t_i]}, whereas for values larger equal a middle year, the data point in \eqn{[t_i, t_{i+1}]} is valid.

The package contains a dataset called \file{my_e0_template} (in \file{extdata} directory) which is a template for user-specified \code{my.e0.file}.
}

\value{
	An object of class \code{bayesLife.mcmc.set} which is a list with two components:
	\item{meta}{An object of class \code{\link{bayesLife.mcmc.meta}}.}
	\item{mcmc.list}{A list of objects of class \code{\link{bayesLife.mcmc}}, one for each MCMC.}
}
\references{
J. L. Chunn, A. E. Raftery, P. Gerland, H. Sevcikova (2013): Bayesian Probabilistic Projections of Life Expectancy for All Countries. Demography 50(3):777-801. <doi:10.1007/s13524-012-0193-x>
}
\author{
Hana Sevcikova, Patrick Gerland contributed to the documentation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get.e0.mcmc}}, \code{\link{summary.bayesLife.mcmc.set}}, \code{\link{e0mcmc.options}}, \code{\link{e0.predict}}.
}
\examples{
\dontrun{
m <- run.e0.mcmc(nr.chains = 1, iter = 5, thin = 1, verbose = TRUE)
summary(m)
m <- continue.e0.mcmc(iter = 5, verbose = TRUE)
summary(m)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
