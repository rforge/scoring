\name{WorldEvents}
\alias{WorldEvents}
\docType{data}
\title{
Forecasts of world events
}
\description{
Probabilistic forecasts of three world events, provided by seven MTurkers.
}
\usage{data("WorldEvents")}
\format{
  A data frame with forecasts of three world events provided by seven Mechanical Turk users.
  \describe{
    \item{\code{forecaster}}{Forecaster ID}
    \item{\code{item}}{Item ID (see details)}
    \item{\code{answer}}{Item resolution (0/1)}
    \item{\code{forecast}}{Forecast associated with outcome 1}
  }
}
\details{
The three forecasted items were:

1. The UK will leave the European Union before the end of 2012.

2. Before Jan 1, 2013, Apple will announce it has sold more than 10 million iPad minis.

3. Japan's nuclear plant in Tsuruga will remain idle between June 1 and December 31, 2012.

For each item, \code{outcome=1} implies that the item text did occur and \code{outcome=0} implies that the item text did not occur.  Forecasts were provided on Dec 20, 2012.
}
\source{
Unpublished data provided by Ed Merkle.
}
\examples{
data("WorldEvents")

## Average forecast for each item
with(WorldEvents, tapply(forecast, item, mean))

## Brier scores
bs <- calcscore(answer ~ forecast, data = WorldEvents, bounds=c(0,1))
}
\keyword{datasets}
