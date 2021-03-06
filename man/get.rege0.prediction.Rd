\name{get.rege0.prediction}
\alias{get.rege0.prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessing Subnational Prediction Objects
}
\description{
Retrieve subnational (regional) prediction results produced by \code{\link{e0.predict.subnat}}, either for one country or for all available countries.
}

\usage{
get.rege0.prediction(sim.dir, country = NULL, method = "ar1", joint.male = FALSE)
}

\arguments{
  \item{sim.dir}{Simulation directory of the subnational predictions. It corresponds to the argument \code{output.dir} in \code{\link{e0.predict.subnat}}.
}
  \item{country}{Numerical country code. If it is not given, all available subnational predictions are retrieved.}
  \item{method}{Method used for generating the projections. It corresponds to the \code{method} argument in \code{\link{e0.predict.subnat}}.}
  \item{joint.male}{Logical. If \code{TRUE}, the function is applied to a male prediction that was generated using the joint female-male gap model implemented in the function \code{\link{e0.jmale.predict.subnat}}.}
}
\details{
Predictions for country \eqn{x} are assumed to be stored in \dQuote{\code{sim.dir}/subnat_\code{method}/c\eqn{x}}.
}

\value{
If argument \code{country} is given, the function returns an object of class \code{\link{bayesLife.prediction}}. If it is \code{NULL}, it returns a list of such objects. Names of the list items are the country codes.
}

\seealso{
\code{\link{e0.predict.subnat}}
}
\examples{
# Subnational example data
my.sube0.file <- file.path(find.package("bayesLife"), 'extdata', 'subnational_e0_template.txt')
sube0 <- read.delim(my.sube0.file, check.names = FALSE)
countries <- unique(sube0[, c("country_code", "country_name")])

# Directory with national projections (contains 30 trajectories for each country)
nat.dir <- file.path(find.package("bayesLife"), "ex-data", "bayesLife.output")

# Subnational projections for all three countries,
# including male projections where female 
# data is used only for demonstration purposes
# (my.e0M.file should contain male e0).
subnat.dir <- tempfile()
e0.predict.subnat(countries$country_code, my.e0.file = my.sube0.file,
    sim.dir = nat.dir, output.dir = subnat.dir, start.year = 2013,
    predict.jmale = TRUE, my.e0M.file = my.sube0.file)
    
# Retrieve results for all countries
preds <- get.rege0.prediction(subnat.dir)
names(preds)

# View tables of subregions for each country
for(i in 1:nrow(countries)) {
  cat("\n\n", countries$country_name[i], "\n")
  print(get.countries.table(preds[[as.character(countries$country_code[i])]]))
}
# Quantiles for individual subregions for female
e0.trajectories.table(preds[["36"]], "Victoria")

# Retrieve results for one country (Canada)
pred <- get.rege0.prediction(subnat.dir, 124)
e0.trajectories.plot(pred, "Quebec", both.sexes = TRUE)

# Retrieve only male results
predM <- get.rege0.prediction(subnat.dir, 124, joint.male = TRUE)
e0.trajectories.table(predM, "Quebec")

# cleanup
unlink(subnat.dir)

# See more examples in ?e0.predict.subnat
}

\keyword{ manip }% use one of  RShowDoc("KEYWORDS")

