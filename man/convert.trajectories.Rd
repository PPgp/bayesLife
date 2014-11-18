\name{convert.e0.trajectories}
\alias{convert.e0.trajectories}

\title{
Converting Trajectories of Life Expectancy into ACSII Files}
\description{
Converts trajectories of the life expectancy stored in a binary format into two CSV files of a UN-specific format.
}
\usage{
convert.e0.trajectories(dir = file.path(getwd(), "bayesLife.output"), 
    n = 1000, output.dir = NULL, verbose = FALSE)
}
\arguments{
  \item{dir}{Directory containing the prediction object. It should correspond to the \code{output.dir} argument of the \code{\link{e0.predict}} function.}
  \item{n}{Number of trajectories to be stored. It can be either a single number or the word \dQuote{all} in which case all available trajectories are converted.}
  \item{output.dir}{Directory in which the resulting files will be stored. If \code{NULL} the same directory is used as for the prediction. Otherwise, if the directory contains joint predictions for both sexes, the ouptuts are stored into subdirectories \file{F} and \file{M}.}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
The function creates two files per sex. One is called \dQuote{ascii_trajectories.csv}, it is a comma-separated table with the following columns: 
\itemize{\item{\dQuote{LocID}: }{country code} 
		\item{\dQuote{Period}: }{prediction interval, e.g. 2015-2020} 
		\item{\dQuote{Year}: }{middle year of the prediction interval}
		\item{\dQuote{Trajectory}: }{identifier of the trajectory}
		\item{\dQuote{e0}: }{life expectancy}
}

The second file is called \dQuote{ascii_trajectories_wide.csv}, it is also a comma-separated table and it contains the same information as above but in a \sQuote{transposed} format. I.e. the data for one country are ordered in columns, thus, there is one column per country. The country columns are ordered alphabetically. 

If \code{n} is smaller than the total number of trajectories, the trajectories are selected using equal spacing.}

\note{This function is automatically called from the \code{\link{e0.predict}} function, therefore in standard cases it will not be needed to call it directly. However, it can be useful for example, if different number of trajectories are to be converted, without having to re-run the prediction.}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{write.e0.projection.summary}}, \code{\link{e0.predict}}
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")
pred.dir <- file.path(getwd(), "exampleLEpred")

# stores 10 trajectories out of 26 (2x(60-35)/2 + 1) into 
# exampleLEpred/predictions/ascii_trajectories.csv
e0.predict(sim.dir=sim.dir, output.dir=pred.dir, 
            burnin=35, save.as.ascii=10, verbose=TRUE)
            
# stores all 26 trajectories into the current directory
convert.e0.trajectories(dir=pred.dir, n="all", output.dir=".", verbose=TRUE)

# Note: If the output.dir argument in e0.predict is omitted, 
# call convert.e0.trajectories with dir=sim.dir 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ IO }

