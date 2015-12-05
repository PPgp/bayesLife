\name{write.e0.projection.summary}
\alias{write.e0.projection.summary}
\title{
Writing Projection Summary Files
}
\description{
The function creates two files containing projection summaries, such as the median, the lower and upper bound of the 80 and 90\% probability intervals, respectively, and the constant variant. One file is in a user-friendly format, whereas the other is in a UN-specific format with internal coding of the time and the variants.
}
\usage{
write.e0.projection.summary(dir = file.path(getwd(), "bayesLife.output"), 
    output.dir = NULL, revision = NULL, adjusted = FALSE)
}
\arguments{
  \item{dir}{Directory containing the prediction object. It should correspond to the \code{output.dir} argument of the \code{\link{e0.predict}} function.}
  \item{output.dir}{Directory in which the resulting file will be stored. If \code{NULL} the same directory is used as for the prediction.}
  \item{revision}{UN revision number. If \code{NULL} it is determined from the corresponding WPP year: WPP 2008 corresponds to revision 13, every subsequent WPP increases the revision number by one. Used as a constant in the second file only.}
   \item{adjusted}{Logical. By default the function writes summary using the original BHM projections. If the projection medians are adjusted (using e.g. \code{\link{e0.median.set}}), setting this argument to \code{TRUE} causes writing the adjusted projections.}
}
\details{
The first file that the function creates is called \file{projection_summary_user_friendly.csv}, it is a comma-separated table with the following columns: 
\itemize{\item{\dQuote{country_name}: }{country name}
		 \item{\dQuote{country_code}: }{country code}
		 \item{\dQuote{variant}: }{name of the variant, such as \dQuote{median}, \dQuote{lower 80}, \dQuote{upper 80}, \dQuote{lower 95}, \dQuote{upper 95}, \dQuote{constant}}
		\item{period1: }{e.g. \dQuote{2010-2015}: life expectancy for the first time period}
		\item{period2: }{e.g. \dQuote{2015-2020}: life expectancy for the second time period}
		\item{... }{further columns with life expectancy projections}
	}

The second file, called \file{projection_summary.csv}, also comma-separated table, contains the same information as above in a UN-specific format:
\itemize{\item{\dQuote{RevID}: }{revision number, passed to the function as an argument;}
		  \item{\dQuote{VarID}: }{variant identifier, extracted from the \code{\link[bayesTFR]{UN_variants}} dataset in the \pkg{bayesTFR} package;}
		  \item{\dQuote{LocID}: }{country code;}
		  \item{\dQuote{TimeID}: }{time identifier, extracted from the \code{\link[bayesTFR]{UN_time}} dataset in the \pkg{bayesTFR} package;}
		  \item{\dQuote{e0}: }{the life expectancy for this variant, location and time period.}
}
If the simulation directory contains joint male predictions, summary files for those are created as well. In such a case, if \code{output.dir} is given, separate subdirectories for female and male are created.
}

\author{
Hana Sevcikova
}
\note{
This function is automatically called from the \code{\link{e0.predict}} and \code{\link{e0.jmale.predict}} functions, therefore in standard cases it will not be needed to call it directly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{convert.e0.trajectories}}, \code{\link{e0.predict}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
