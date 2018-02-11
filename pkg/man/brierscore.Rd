\name{brierscore}
\alias{brierscore}

\title{
Calculate Brier Scores And Decompositions
}
\description{
Calculate Brier scores, average Brier scores by a grouping variable, and
Brier score decompositions for two-alternative forecasts.
}
\usage{
brierscore(object, data, group = NULL, decomp = FALSE, bounds = NULL,
           reverse = FALSE, wt = NULL, decompControl = list())
}

\arguments{
  \item{object}{an object of class "formula", of the form
    \code{outcome ~ forecast}.  See \code{calcscore()} documentation for further details.}
  \item{data}{an optional data frame or list containing the
          variables in the formula.  If not found in \code{data}, the
          variables are taken from the environment from which
	  \code{calcscore} is called.}
  \item{group}{the name of a grouping variable within \code{data}, which
  is used to calculate average Brier score by group.}
  \item{decomp}{if \code{TRUE}, Brier score decompositions are
  calculated.}
  \item{bounds}{a vector of length 2 corresponding to the desired
    minimum and maximum Brier score, respectively.}
  \item{reverse}{if \code{FALSE} (default), smaller scores imply better
    forecasts.  If \code{TRUE}, larger scores imply better forecasts.}
  \item{wt}{a vector of weights for computing a weighted Brier score. If \code{NULL}, the weights are set to 1/n, where n is the number of forecasts (this corresponds to a simple average Brier score).}
  \item{decompControl}{a list of additional settings for the Brier
  decomposition. See options below.}
}
\details{
If \code{decomp=TRUE} or \code{group} is supplied, the function returns
  a list (see value section).  Otherwise, the function returns a numeric
  vector containing the Brier score associated with each forecast.

Some \code{decompControl} arguments are specifically designed for forecasting tournaments and may not be useful in other situations. Possible arguments for \code{decompControl} include:
\describe{
  \item{wt}{A vector of weights, for performing a weighted Brier
    decomposition (could also use the simple \code{wt} argument).}
  \item{qid}{A vector of question ids, for use with the \code{qtype} argument.}
  \item{bin}{If \code{TRUE} (default), forecasts are binned prior to decomposition. If \code{FALSE}, the original forecasts are maintained.}
  \item{qtype}{A data frame with columns \code{qid}, \code{ord}, \code{squo}. For each unique question id in the \code{qid} argument above, this describes whether or not the question is ordinal (1=yes,0=no) and whether or not the question has a "status quo" interpretation (1=yes,0=no).}
  \item{scale}{Should Brier components be rescaled, such that 1 is
    always best and 0 is always worst? Defaults to \code{FALSE}.}
  \item{roundto}{To what value should forecasts be rounded (necessary
    for Murphy decomposition)? Defaults to
    .1, meaning that forecasts are rounded to the nearest .1.}
  \item{binstyle}{Method for ensuring that each forecast sums to 1. If
    equal to 1 (default), the smallest forecast is one minus the sum of
    the other forecasts. If equal to 2, the forecast furthest from its
    rounded value is one minus the sum of other forecasts.}
  \item{resamples}{Desired number of Brier resamples (useful for questions with inconsistent alternatives). Defaults to 0; see Merkle & Hartman reference for more detail.}
}
}
\value{
Depending on input arguments, \code{brierscore} may return an object of
  class \code{numeric} containing raw Brier scores.  It may also return
  a list containing the objects below.
  
\item{rawscores}{an object of class \code{numeric} containing raw Brier
  scores for each forecast.}
\item{brieravg}{an object of class \code{numeric} containing average Brier
  scores for each unique value of \code{group}. If \code{wt} was
  supplied, this is a weighted sum. Otherwise, it is a simple average
  (equal weights summing to 1).}
\item{decomp}{an object of class \code{matrix} containing Brier score
  decompositions and mean Brier scores for each unique value of \code{group}.}
}
\references{
Brier, G. W. (1950).  Verification of forecasts expressed in terms of
probability.  \emph{Monthly Weather Review, 78},  1-3.

Merkle, E. C. & Hartman, R. (2018). Weighted Brier score decompositions
for topically heterogenous forecasting tournaments. \emph{Working paper}.

Murphy, A. H. (1973).  A new vector partition of the probability score.
\emph{Journal of Applied Meteorology, 12}, 595-600.

Yates, J. F. (1982).  External correspondence: Decompositions of the
mean probability score.  \emph{Organizational Behavior and Human
  Performance, 30}, 132-156.

Young, R. M. B. (2010). Decomposition of the Brier score for weighted
forecast-verification pairs. \emph{Quarterly Journal of the Royal
  Meteorological Society, 136}, 1364-1370.
}
\author{
Ed Merkle
}


\seealso{
\code{\link{calcscore}}
}
\examples{
data("WorldEvents")
## Raw Brier scores
brier1 <- brierscore(answer ~ forecast, data=WorldEvents)
## Raw Brier scores plus group means and decompositions
brier2 <- brierscore(answer ~ forecast, data=WorldEvents,
                     group="forecaster", decomp=TRUE)
## Obtain Brier scores via calcscore
brier3 <- calcscore(answer ~ forecast, data=WorldEvents,
                    param=2, fam="pow")
all.equal(brier1, brier3)
}
