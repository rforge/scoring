\name{WeatherProbs}
\alias{WeatherProbs}
\docType{data}
\title{
Three-category weather forecasts
}
\description{
Probabilistic forecasts from the U.S. National Oceanic and Atmospheric
Administration, concerning below/near/above average temperatures and
below/near/above median precipitation.
}
\usage{data("WeatherProbs")}
\format{
  A data frame with 8976 observations on the following 11 variables.
  \describe{
    \item{\code{stn}}{Station World Meteorological Organization (WMO) number}
    \item{\code{made}}{Forecast issuance date}
    \item{\code{valid}}{Center of forecast valid period}
    \item{\code{tblw}}{Probability of below normal temperatures}
    \item{\code{tnrm}}{Probability of near normal temperatures}
    \item{\code{tabv}}{Probability of above normal temperatures}
    \item{\code{tcat}}{Realized temperature category (1=below, 2=near, 3=above)}
    \item{\code{pblw}}{Probability of below median precipitation}
    \item{\code{pnrm}}{Probability of near median precipitation}
    \item{\code{pabv}}{Probability of above median precipitation}
    \item{\code{pcat}}{Realized precipitation category (1=below, 2=near, 3=above)}
  }
}
\details{
The forecasts are valid for a period of 6 to 10 days from the date
that the forecast was made.  The forecasts were supplied every
weekday during April, 2009, and they specifically predict the average
temperature or total precipitation for the entire valid period.
}
\source{
Data were obtained from
http://www.cpc.ncep.noaa.gov/products/archives/short_range/ (see URL in references).
}
\references{
See
\url{http://www.cpc.ncep.noaa.gov/products/archives/short_range/README.6-10day.txt}
for more details on the data.

For an application of similar data (different dates, same source), see:

Wilks, D. S. (in press).  The calibration simplex: A generalization of
the reliability diagram for 3-category probability forecasts.
\emph{Weather and Forecasting}.
}
\examples{
data("WeatherProbs")

## Brier score for temperature forecasts
## (Warning arises because some forecast rows don't sum to 1.)
res <- calcscore(tcat ~ tblw + tnrm + tabv, data=WeatherProbs,
                 bounds=c(0,1))

## Ordered Brier score for temperature forecasts
res2 <- calcscore(tcat ~ tblw + tnrm + tabv, data=WeatherProbs,
                  bounds=c(0,1), ordered=TRUE)

## Spherical score for temperature forecasts
res3 <- calcscore(tcat ~ tblw + tnrm + tabv, data=WeatherProbs,
                  fam="sph", bounds=c(0,1))

## Average scores by station
avgbrier <- with(WeatherProbs, tapply(res, stn, mean))
avgobrier <- with(WeatherProbs, tapply(res2, stn, mean))
avgsph <- with(WeatherProbs, tapply(res3, stn, mean))

## Conclusions vary across Brier and ordinal Brier scores
plot(avgbrier, avgobrier, pch=20, xlab="Brier", ylab="Ordinal Brier")
}
\keyword{datasets}
