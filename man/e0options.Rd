\name{e0options}
\alias{e0options}
\alias{e0mcmc.options}
\alias{e0pred.options}
\alias{using.bayesLife}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Global options
}
\description{
Setting and retrieving global options.
}
\usage{
using.bayesLife()
e0options()
e0mcmc.options(...)
e0pred.options(...)
}

\arguments{
  \item{...}{Arguments in \code{tag = value} form, or names of options to retrieve.
}
}
\details{
Function \code{using.bayesLife} sets all global options to their default values. Function \code{e0options} is used to get all options as a named list. 

The global options are divided into two main categories, namely options used for MCMC estimations and options used for predictions. To set or retrieve options of the first category, use \code{e0mcmc.options} (see section MCMC Options below), while the second category is controlled by \code{e0pred.options} (see section Prediction Options below). 

Many options are in form of a list and it is possible to overwrite only single elements of the list. However, if an option is a vector, all elements of the vector have to be defined when updating (see Example).
}

\value{
\code{e0options} returns a list of all global options.\cr
\code{e0mcmc.options}, when called with no argument, it returns a list of options related to the MCMC estimation.\cr
\code{e0pred.options}, when called with no argument, it returns a list of options related to the prediction. \cr
For both, \code{e0mcmc.options} and \code{e0pred.options}, when a specific option is queried, it returns the value of that option. When an option is set, a list of the previous values of all MCMC/prediction options is returned invisibly. 
}

\section{MCMC Options}{
\describe{
  \item{a}{vector of the \eqn{a_1}, \dots ,\eqn{a_6} parameters, which are the prior means of the world-level parameters (\eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}, \eqn{k}, \eqn{z}).}
  \item{delta}{vector of the \eqn{\delta_1}{delta_1}, \dots ,\eqn{\delta_6}{delta_6} parameters, 
  	which are the prior standard deviations of the world-level parameters (\eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}, \eqn{k}, \eqn{z}).}
  \item{tau}{vector of the \eqn{\tau_1}{tau_1}, \dots ,\eqn{\tau_6}{tau_6} parameters, which is the square root rate of the prior Gamma distribution of the world-level parameters (\eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}, \eqn{\lambda_k}{lambda_k}, \eqn{\lambda_z}{lambda_z}).}
  \item{Triangle}{list with elements:
    \describe{
        \item{ini}{list with elements:
            \describe{
                \item{T1, T2, T3, T4}{initial values for \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}. If not \code{NULL}, then each element should be of the same length as the number of MCMC chains. If it is \code{NULL}, the initial values are equally spaced between \code{ini.low} and \code{ini.up} for the respective parameter. By default in the estimation, if there is just one chain, the initial value is the middle point of the interval.}
            }}
        \item{ini.low, ini.up}{vectors of length four. They are the lower and upper bounds for initial values of \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}. An \eqn{i}-th item is only used if \code{ini$Ti} is \code{NULL}.}
       \item{prior.low, prior.up}{vectors of length four. They are the lower and upper bounds for the prior (truncated normal) distribution of \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}.}   
        \item{slice.width}{vector of length four defining the slice width for MCMC slice sampling for the four parameters, \eqn{\Delta_1}{Triangle_1}, \dots, \eqn{\Delta_4}{Triangle_4}.}
    }}
    \item{k, z}{lists with elements:
        \describe{
        \item{ini}{vector of initial values for \eqn{k} (\eqn{z}). Its length (if not \code{NULL}) should correspond to the number of MCMC chains. By default, the initial values are equally spaced between \code{ini.low} and \code{ini.up}. In case of one chain, the initial value is by default the middle point of the interval.}
        \item{ini.low, ini.up}{single value giving the lower and upper bounds for initial values of \eqn{k} (\eqn{z}). It is only used if \code{ini} is \code{NULL}. Regarding defaults for the \eqn{z} parameter, see Note below.}
         \item{prior.low, prior.up}{single value giving the lower and upper bounds for the prior (truncated normal) distribution  of \eqn{k} (\eqn{z}). Regarding defaults for the \eqn{z} parameter, see Note below.}
         \item{slice.width}{single value giving the slice width for MCMC slice sampling of the \eqn{z} parameter (not available for \eqn{k}).}
        }}
    \item{lambda}{list with elements:
        \describe{
        \item{ini}{list with elements:
            \describe{
                \item{T1, T2, T3, T4}{initial values for \eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}. Each element should be of the same length as the number of MCMC chains. If it is \code{NULL}, the initial values are equally spaced between \code{ini.low} and \code{ini.up} of the respective parameter. By default, if there is just one chain, the value is the middle point of the interval.}
            }}
            \item{ini.low, ini.up}{vectors of length four. They are the lower and upper bounds for initial values of \eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}. An \eqn{i}-th item is only used if \code{ini$Ti} is \code{NULL}.}
            \item{slice.width}{vector of length four defining the slice width for MCMC slice sampling for the four parameters, \eqn{\lambda_1}{lambda_1}, \dots, \eqn{\lambda_4}{lambda_4}.}
            }}
    \item{lambda.k, lambda.z}{lists with elements:
        \describe{     
        \item{ini}{vector of initial values for \eqn{\lambda_k}{lambda_k} (\eqn{\lambda_z}{lambda_z}). Its length (if not \code{NULL}) should correspond to the number of MCMC chains. By default, the initial values are equally spaced between \code{ini.low} and \code{ini.up}. In case of one chain, the initial value is by default the middle point of the interval.}
       \item{ini.low, ini.up}{single value giving the lower and upper bounds for initial values of \eqn{\lambda_k}{lambda_k} (\eqn{\lambda_z}{lambda_z}). It is only used if \code{ini} is \code{NULL}.} 
        \item{slice.width}{single value giving the slice width for MCMC slice sampling of the \eqn{\lambda_z}{lambda_z} parameter (not available for \eqn{\lambda_k}{lambda_k}).}
        }}
    \item{omega}{list with elements:
            \describe{ 
                \item{ini}{vector of initial values for \eqn{\omega}{omega}. Its length (if not \code{NULL}) should correspond to the number of MCMC chains. By default, the initial values are equally spaced between \code{ini.low} and \code{ini.up}. In case of one chain, the initial value is by default the middle point of the interval.}
                \item{ini.low, ini.up}{single value giving the lower and upper bounds for initial values of \eqn{\omega}{omega}. It is only used if \code{ini} is \code{NULL}.}
    }}
    \item{Triangle.c}{list with elements:
        \describe{  
            \item{ini.norm}{list with elements:
                \describe{
                    \item{mean, sd}{vectors of size four. They correspond to the means and standard deviations, respectively, for the initial values of the country-specific parameters \eqn{\Delta_1^c}{Triangle_1_c}, \dots, \eqn{\Delta_4^c}{Triangle_4_c} which are drawn from a truncated normal distribution with bounds defined by \code{prior.low} and \code{prior.up}.}
                    }}
            \item{prior.low, prior.up}{vectors of length four. They are the lower and upper bounds for the prior (truncated normal) distribution of country-specific \eqn{\Delta_1^c}{Triangle_1_c}, \dots, \eqn{\Delta_4^c}{Triangle_4_c}.} 
            \item{slice.width}{vector of length four defining the slice width for MCMC slice sampling of the country-specific \eqn{\Delta_1^c}{Triangle_1_c}, \dots, \eqn{\Delta_4^c}{Triangle_4_c}.}
        }}
    \item{k.c, z.c}{list with elements:
        \describe{  
            \item{ini.norm}{named vector of length two, called \dQuote{mean} and \dQuote{sd}. The elements correspond to the means and standard deviations, respectively, for the initial values of the country-specific parameters \eqn{k^c} (\eqn{z^c}) which are drawn from a normal distribution truncated between \code{prior.low} and \code{prior.up}.}
            \item{prior.low, prior.up}{single values giving the lower and upper bounds for the prior (truncated normal) distribution  of country-specific \eqn{k^c} (\eqn{z^c}). Regarding defaults for \eqn{z^c}, see Note below.}
            \item{slice.width}{single value giving the slice width for MCMC slice sampling of the \eqn{k^c} (\eqn{z^c}) parameter.}
    }}
    \item{nu}{the shape parameter of the Gamma distributions of all \eqn{\lambda}{lambda} parameters is \code{nu/2}.}
    \item{dl.p1, dl.p2}{values of the parameters \eqn{p_1} and \eqn{p_2} of the double logistic function.}
    \item{sumTriangle.lim}{lower and upper limits for the sum of the \eqn{\Delta_i}{Triangle_i} parameters. MCMC proposals that are outside of this limit are rejected. It is applied to both, the world parameters as well as the country specific parameters.}
    \item{world.parameters}{named vector where names are the world parameters and values are the number of sub-parameters. For example, \eqn{\Delta}{Triangle} has 4 sub-parameters, while \eqn{k} and \eqn{z} are both just one parameter.}
    \item{country.parameters}{named vector where names are the country-specific parameters and values are the number of sub-parameters.}
    \item{outliers}{ranges for determining outliers in the historical data. If \code{outliers=c(x, y)} then any increase in life expectancy smaller than \code{x} or larger than \code{y} is considered as an outlier and removed from the estimation.}
    \item{buffer.size}{buffer size (in number of [thinned] iterations) for keeping data in the memory. The smaller the \code{buffer.size} the more often will the process access the hard disk and thus, the slower the run. On the other hand, the smaller the \code{buffer.size} the less data will be lost in case of failure.}
    \item{auto.conf}{list containing a configuration for an \sQuote{automatic} run. All items in this list must be integer values. The option is only used if the argument \code{iter} in \code{\link{run.e0.mcmc}} is set to \sQuote{auto} (see description of argument \code{iter} in \code{\link{run.e0.mcmc}}). The list contains the following elements:
    \describe{  
            \item{iter}{gives the number of iterations in the first chunk of the MCMC simulation.}
            \item{iter.incr}{gives the number of iterations in the following chunks.}
            \item{nr.chains}{gives the number of chains in all chunks of the MCMC simulation.}
            \item{thin, burnin}{used in the convergence diagnostics following each chunk.}
            \item{max.loops}{controls the maximum number of chunks.}
            }}
  \item{country.overwrites}{This option allows to overwrite some of the prior parameters for specific countries. If it is not \code{NULL} it should be a data frame with an obligatory column \sQuote{country_code}. Each row then corresponds to one country. Other columns can be \sQuote{k.c.prior.low}, \sQuote{k.c.prior.up}, \sQuote{z.c.prior.low}, \sQuote{z.c.prior.up}, \sQuote{Triangle_\eqn{x}.c.prior.low} and \sQuote{Triangle_\eqn{x}.c.prior.up} where \eqn{x} can be an integer from 1 to 4.}
}}

%\section{Prediction Options}{blabla}



\note{
	Parameter \eqn{z} determines the asymptote in gains in life expectancy. The following text gives an explanation for the choice of  upper limits on \eqn{z}-related parameters: 
	
	The pace of improvement and the asymptotic limit in future gains in female life expectancy vary for each projected trajectory, but ultimately is informed and constrained by the finding that the rate of increase of maximum female life expectancy over the past 150 year has been highly linear (2a, 2b) (i.e., about 2.4 years per decade), albeit at slightly lower pace once the leading countries started to exceed 75 years of female life expectancy at birth in the 1960s (3) (about 2.26 years of gains per decade). By assuming that the asymptotic average rate of increase in life expectancy is nonnegative, life expectancy is assumed to continually increase (on average), and no limit is imposed to life expectancy in the foreseeable future. The increase in maximum female life span among countries with highest life expectancy and reliable data on very old age provide further guidance on future rate of progress which has also been increasingly linear at least since the 1970s (4a-4c) (about 1.25 years per decade for countries like Sweden and Norway), and is used to inform the asymptotic average rate of increase in female life expectancy used in the 2012 WPP Revision. To set the posterior median to an annual gain of 0.125 year (or 5-year gain of 0.625 in this context) the upper bound value of 0.653 is used for the world prior (\eqn{z}) and country-specific prior (\eqn{z_c}) as default values in the estimation of the double-logistic parameters.
}


\author{
Hana Sevcikova, Patrick Gerland contributed to the documentation.
}

\references{
(1) J. L. Chunn, A. E. Raftery, P. Gerland, H. Sevcikova (2013): Bayesian Probabilistic Projections of Life Expectancy for All Countries. Demography 50(3):777-801. <doi:10.1007/s13524-012-0193-x>

(2a)  Oeppen J, and J.W. Vaupel (2002) Broken limits to life expectancy. Science 296:1029-1031. 

(2b) Vaupel, J.W. and K.G.V. Kistowski. 2005. Broken Limits to Life Expectancy. Ageing Horizons (3):6-13.

(3) Vallin, J., and F. Mesle (2009). The Segmented Trend Line of Highest Life Expectancies. Population and Development Review, 35(1), 159-187. doi:10.1111/j.1728-4457.2009.00264.x

(4a) Wilmoth, J. R., L. J. Deegan, H. Lundstrom, and S. Horiuchi (2000). Increase of maximum life-span in Sweden, 1861-1999. Science, 289(5488), 2366-2368. 

(4b) Wilmoth, J. R. and J-M. Robine. (2003). The world trend in maximum life span, in: J. R. Carey and S. Tuljapurkar (eds.), Life Span: Evolutionary, Ecological, and Demographic Perspectives, supplement to vol. 29, Population and Development Review, pp. 239-257. 

(4c) Wilmoth, J. R. and N. Ouellette (2012). Maximum human lifespan: Will the records be unbroken?, Paper presented at the European Population Conference, Stockholm, Sweden, 13-16 June.

}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.e0.mcmc}}, \code{\link{e0.predict}}
}
\examples{
e0mcmc.options("z", "Triangle")
# Set new z$ini.up and Triangle$prior.up
# Modifying single elements of the z-list and Triangle-list.
# However, Triangle$prior.up is a vector and needs all four values.
e0mcmc.options(z = list(ini.up = 0.8), Triangle = list(prior.up = rep(120, 4)))
e0mcmc.options("z", "Triangle")

# revert to defaults
using.bayesLife()
e0mcmc.options("z", "Triangle")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~manip }% use one of  RShowDoc("KEYWORDS")

