\name{e0.map}
\alias{e0.map}
\alias{e0.map.all}
\alias{get.e0.map.parameters}
\alias{e0.map.gvis}
\title{
World Map of the Life Expectancy 
}
\description{
Generates a world map of the life expectancy for given quantile and projection or estimation period.
}
\usage{
e0.map(pred, \dots)
	
e0.map.all(pred, output.dir, output.type = "png", 
    e0.range = NULL, nr.cats = 50, same.scale = TRUE, 
    quantile = 0.5, file.prefix = "e0wrldmap_", \dots)
			
get.e0.map.parameters(pred, e0.range = NULL, 
    nr.cats = 50, same.scale = TRUE, quantile = 0.5, \dots)
    
e0.map.gvis(pred, \dots)
}

\arguments{
  \item{pred}{Object of class \code{\link{bayesLife.prediction}}.}
 \item{output.dir}{Directory into which resulting maps are stored.}
 \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
 \item{e0.range}{Range of the life expectancy to be displayed. It is of the form \code{c(}\var{e0.min}, \var{e0.max}\code{)}.
 By default, the whole range is considered. Note that countries with values outside of the given range will appear white.}
  \item{nr.cats}{Number of color categories.}
  \item{same.scale}{Logical controlling if maps for all years of this prediction object should be on the same color scale.}
   \item{quantile}{Quantile for which the map should be generated. It must be equal to one of the values in \code{dimnames(pred$quantiles[[2]])}, i.e. 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 1. Value 0.5 corresponds to the median.}
  \item{file.prefix}{Prefix for file names.}
  \item{\dots}{In \code{e0.map}, \dots are all arguments that can be passed to \code{\link[bayesTFR]{tfr.map}}, such as \code{quantile}, \code{year}, \code{projection.index}, \code{par.name}, \code{adjusted}, \code{device}, \code{main}, \code{device.args}, and \code{data.args}. In \code{e0.map.gvis}, \dots are all arguments that can be passed to \code{\link[bayesTFR]{tfr.map.gvis}}. In addition, the first three functions accept
arguments passed to the \code{\link[rworldmap]{mapCountryData}} function of the \pkg{rworldmap} package.
}
}
\details{
\code{e0.map} creates a single map for the given time period and quantile. \code{e0.map.all} generates a sequence of maps, namely one for each projection period. If the package \pkg{fields} is installed, a color bar legend at the botom of the map is created.

Function \code{get.e0.map.parameters} can be used in combination with \code{e0.map}. (Note that \code{get.e0.map.parameters} is called from inside of \code{e0.map.all}.) It sets breakpoints for the color scheme using quantiles of a fitted gamma distribution.

Function \code{e0.map.gvis} creates an interactive map using the \pkg{googleVis} package and opens it in an internet browser. It also generates a table of the mapped values that can be sorted by columns interactively in the browser. 

By default, both \code{e0.map} and \code{e0.map.gvis} produce maps of life expectancy. Alternatively, the functions can be used to plot country-specific MCMC parameters into a world map. They are given by the argument \code{par.name}. One can pass any value from \code{\link{e0.parameter.names.cs.extended}()}.
}
\value{
\code{get.e0.map.parameters} returns a list with elements:
  \item{pred}{The object of class \code{\link{bayesLife.prediction}} used in the function.}
  \item{quantile}{Value of the argument \code{quantile}.}
  \item{catMethod}{If the argument \code{same.scale} is \code{TRUE}, this element contains breakpoints for categorization. It is generated from a fitted gamma distribution. Otherwise, it is \code{NULL}.}
  \item{numCats}{Number of categories.}
  \item{coulourPalette}{Subset of the rainbow palette, starting from dark blue and ending at red.}
  \item{\dots}{Additional arguments passed to the function.}
}

\author{
Hana Sevcikova, Adrian Raftery
}

\seealso{\code{\link[bayesTFR]{tfr.map}}}

\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
pred <- get.e0.prediction(sim.dir=sim.dir)
# Uses heat colors and seven categories by default
e0.map(pred)
# Uses more colors with more suitable categorization
params <- get.e0.map.parameters(pred)
do.call('e0.map', params)
# Another projection year on the same scale
do.call('e0.map', c(list(year=2043), params))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
