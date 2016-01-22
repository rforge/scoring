\name{logscore}
\alias{logscore}

\title{
Calculate Logarithmic Scores
}
\description{
Calculate logarithmic scores and average logarithmic scores by a grouping variable.
}
\usage{
logscore(object, data, group = NULL, reverse = FALSE)
}

\arguments{
  \item{object}{an object of class "formula", of the form
    \code{outcome ~ forecast}.  See \code{calcscore()} documentation for further details.}
  \item{data}{an optional data frame or list containing the
          variables in the formula.  If not found in \code{data}, the
          variables are taken from the environment from which
	  \code{calcscore} is called.}
  \item{group}{the name of a grouping variable within \code{data}, which
    is used to calculate average log score by group.}
  \item{reverse}{if \code{FALSE} (default), smaller scores imply better
    forecasts.  If \code{TRUE}, larger scores imply better forecasts.}
}
\details{
If \code{group} is supplied, the function returns
  a list (see value section).  Otherwise, the function returns a numeric
  vector containing the log score associated with each forecast.

The argument \code{bounds} is not available because
the upper bound of the logarithmic score is infinity.  If one wants a
bounded rule that approximates the logarithmic rule, try using
\code{calcscore()} with \code{fam="pow"} and \code{param=1.001}.
}
\value{
Depending on input arguments, \code{logscore} may return an object of
  class \code{numeric} containing raw logarithmic scores.  It may also return
  a list containing the objects below.
  
\item{rawscores}{an object of class \code{numeric} containing raw log
  scores for each forecast.}
\item{mnlog}{an object of class \code{numeric} containing mean log
  scores for each unique value of \code{group}.}
}
\references{
Toda, M. (1963).  Measurement of subjective probability distributions. 
  ESD-TDR-63-407, Decision Sciences Laboratory, L. G. Hanscom Field,
  Bedford, Mass.

Shuford, E. H., Albert, A., & Massengill, H. E. (1966).  Admissible
probability measurement procedures.  \emph{Psychometrika, 31}, 125-145.
}
\author{
Ed Merkle
}


\seealso{
\code{\link{calcscore}}
}
\examples{
data("WorldEvents")
## Raw log scores (0 best, infinity worst)
log1 <- logscore(answer ~ forecast, data=WorldEvents)
## Raw log scores (0 best, -infinity worst)
log1 <- logscore(answer ~ forecast, data=WorldEvents,
                 reverse = TRUE)
## Raw log scores plus group means
log2 <- logscore(answer ~ forecast, data=WorldEvents,
                 group="forecaster")
}
