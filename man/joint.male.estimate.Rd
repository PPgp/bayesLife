\name{e0.jmale.estimate}
\alias{e0.jmale.estimate}

\title{
Estimation of the Joint Female-Male Model
}
\description{
The function estimates the joint female-male model of life expectancy, as described in Raftery et al. (2014, 2012) and Lalic (2011). It consist of two equations with t-distributed errors, see Details below.
}
\usage{
e0.jmale.estimate(mcmc.set, countries.index = NULL, 
    estDof.eq1 = TRUE, start.eq1 = list(dof = 2), max.e0.eq1 = 83, 
    estDof.eq2 = TRUE, start.eq2 = list(dof = 2),
    constant.gap.eq2 = TRUE, my.e0.file = NULL, 
    my.locations.file = NULL, verbose = FALSE)
}

\arguments{
  \item{mcmc.set}{Object of class \code{\link{bayesLife.mcmc.set}} containing estimation results of female life expectancy.}
  \item{countries.index}{Index of countries (within the mcmc.set object) to be included in the estimation. By default, all countries included in the estimation of \code{mcmc.set} are used.}
  \item{estDof.eq1, estDof.eq2}{Logical, controlling whether the degrees of freedom of the first and second  equation, respectively, should be estimated. If it is \code{FALSE} the degrees of freedom are set using the arguments \code{start.eq1} and \code{start.eq2}, respectively.}
  \item{start.eq1, start.eq2}{Argument \code{start} of the \code{\link[hett]{tlm}} function of the \pkg{hett} package, which is applied to the first and second equation of the model, respectively. It is a list of possibly four named components, ("beta", "lambda", "dof", "omega"), for the location, scale, degrees of freedom parameters and random scale effects respectively. If \code{estDof.eq1} (\code{estDof.eq2}) is \code{FALSE}, the "dof" component must be given.
}
 \item{max.e0.eq1}{Maximum female life expectancy of records included in the estimation of the first equation (parameter \eqn{M} in Details below).}
 \item{constant.gap.eq2}{Logical. If \code{TRUE} the coefficient of the second equation (\eqn{\gamma_1}) is set to one and the standard deviation is estimated under the assumption of normal distribution.}
 \item{my.e0.file}{File name containing user-specified male time series for one or more countries. The function replaces the corresponding country data from the  WPP dataset by values in this file. Only columns are replaced that match column names of the WPP dataset.}
 \item{my.locations.file}{File name containing user-specified locations if different from the default \code{\link[wpp2017]{UNlocations}} dataset. It should be the same file as passed to \code{\link{run.e0.mcmc}} for female life expectancy.}
 \item{verbose}{Logical switching log messages on and off. If \code{TRUE} summary results from the \code{\link[hett]{tlm}} function of both equations are shown.}
}
\details{
The joint female-male life expectancy model is a model for estimating gaps \eqn{G} between female and male life expectancy. It consists of two parts, see Equation (1) in Raftery et  al. (2012): \cr
1. If \eqn{l_{c,t} \leq M}{l_(c,t) <= M}, then 
\deqn{
	G_{c,t} = \beta_0 + \beta_1 l_{c,1953} + \beta_2 G_{c,t-1} + \beta_3 l_{c,t} + \beta_4 (l_{c,t}-75)_+ + \epsilon_{c,t}
}{
	G_(c,t) = beta_0 + beta_1 * l_(c,1953) + beta_2 * G_(c,t-1) + beta_3 * l_(c,t) + beta_4 * [l_(c,t)-75]_+ + epsilon_(c,t)
}
where \eqn{\epsilon_{c,t}}{epsilon_(c,t)} is iid \eqn{t(\mu=0, \sigma_1^2, \nu_1)}{t(mu=0, sigma_1^2, nu_1)}.

2. If \eqn{l_{c,t} > M}{l_(c,t) > M}, then 
\deqn{
	G_{c,t} = \gamma_1 G_{c,t-1} + \epsilon_{c,t}
}{
	G_(c,t) = gamma_1 * G_(c,t-1) + epsilon_(c,t)
}
where \eqn{\epsilon_{c,t}}{epsilon_(c,t)} is iid \eqn{t(\mu=0, \sigma_2^2, \nu_2)}{t(mu=0, sigma_2^2, nu_2)}.

Here, \eqn{t} is the time and \eqn{c} is the country index. \eqn{G_{c,t}}{G_(c,t)} is the gap for country \eqn{c} at time \eqn{t} and \eqn{l_{c,t}}{l_(c,t)} is the female life expectancy for country \eqn{c} at time \eqn{t}. \eqn{M} can be set in the \code{max.e0.eq1} argument.

Using the \code{\link[hett]{tlm}} function of the \pkg{hett} package, the function  estimates the coefficients \eqn{\beta_i}{beta_i} (\eqn{i=1,\dots,4}) and \eqn{\gamma_1}{gamma_1}, as well as paramteres \eqn{\sigma_j}{sigma_j} (\eqn{j=1,2}) and optionally the degrees of freedom \eqn{\nu_j}{nu_j} (\eqn{j=1,2}). If \code{constant.gap.eq2} is \code{TRUE}, \eqn{\gamma_1}{gamma_1} is set to 1 and \eqn{\epsilon_{c,t}}{epsilon_(c,t)} is iid \eqn{N(\mu=0, \sigma_2^2)}{N(mu=0, sigma_2^2)}.

The \code{mcmc.set} object should be a \code{\link{bayesLife.mcmc.set}} object obtained from a simulation of a female life expectancy. Note that since only the observed data and no MCMC results are used in this estimation, the \code{mcmc.set} object can be obtained from a toy simulation such as in the example below. The function extracts observed data from this object and treats them as \eqn{l_{c,t}}. For the male historical time series, the function takes the male WPP dataset  from the same \pkg{wpp} package as the female data (\code{\link[wpp2017]{e0M}}) and possibly partly replaces it by any user-specified data given in \code{my.e0.file}. 
}
\value{
List with the components, \code{eq1} and \code{eq2}, each containing estimation results from the first and second equation, respectively. These are:
    \item{coefficients}{Estimated coefficients \eqn{\beta_i}{beta_i}.}
    \item{sigma}{Parameter \eqn{\sigma_j}{sigma_j}.}
    \item{dof}{Degrees of freedom \eqn{\nu_j}{nu_j}. If \code{estDof.eq1} (\code{estDof.eq2}) is \code{TRUE} this parameter is estimated, otherwise it is set to the value of \code{start.eq1$dof} (\code{start.eq2$dof}).}
}
\references{
A. E. Raftery,  N. Lalic, P. Gerland (2014). Joint Probabilistic Projection of Female and Male Life Expectancy. Demographic Research, 30:795-822.
	
A. E. Raftery,  N. Li, H. Sevcikova , P. Gerland, G. K. Heilig (2012). Bayesian probabilistic population projections for all countries. Proceedings of the National Academy of Sciences 109:13915-13921.	
	
Lalic, N. (2011). Master's thesis at the Department of Statistics, University of Washington. 
}
\author{
Hana Sevcikova
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{e0.jmale.predict}}
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
m <- get.e0.mcmc(sim.dir)
fit <- e0.jmale.estimate(m, verbose=TRUE)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

