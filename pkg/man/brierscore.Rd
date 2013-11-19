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
brierscore(object, data, group = NULL, decomp = FALSE, bounds = NULL, reverse = FALSE)
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
}
\details{
If \code{decomp=TRUE} or \code{group} is supplied, the function returns
  a list (see value section).  Otherwise, the function returns a numeric
  vector containing the Brier score associated with each forecast.
}
\value{
Depending on input arguments, \code{brierscore} may return an object of
  class \code{numeric} containing raw Brier scores.  It may also return
  a list containing the objects below.
  
\item{rawscores}{an object of class \code{numeric} containing raw Brier
  scores for each forecast.}
\item{mnbrier}{an object of class \code{numeric} containing mean Brier
  scores for each unique value of \code{group}.}
\item{decomp}{an object of class \code{matrix} containing Brier score
  decompositions and mean Brier scores for each unique value of \code{group}.}
}
\references{
Brier, G. W. (1950).  Verification of forecasts expressed in terms of
probability.  \emph{Monthly Weather Review, 78},  1-3.

Murphy, A. H. (1973).  A new vector partition of the probability score.
\emph{Journal of Applied Meteorology, 12}, 595-600.

Yates, J. F. (1982).  External correspondence: Decompositions of the
mean probability score.  \emph{Organizational Behavior and Human
  Performance, 30}, 132-156.
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
