\name{sphscore}
\alias{sphscore}

\title{
Calculate Spherical Scores
}
\description{
Calculate spherical scores and average spherical scores by a grouping variable.
}
\usage{
sphscore(object, data, group = NULL, bounds = NULL, reverse = FALSE)
}

\arguments{
  \item{object}{an object of class "formula", of the form
    \code{outcome ~ forecast}.  See \code{calcscore()} documentation for further details.}
  \item{data}{an optional data frame or list containing the
          variables in the formula.  If not found in \code{data}, the
          variables are taken from the environment from which
	  \code{calcscore} is called.}
  \item{group}{the name of a grouping variable within \code{data}, which
    is used to calculate average spherical score by group.}
  \item{bounds}{a vector of length 2 corresponding to the desired
    minimum and maximum spherical score, respectively.}
  \item{reverse}{if \code{FALSE} (default), smaller scores imply better
    forecasts.  If \code{TRUE}, larger scores imply better forecasts.}
}
\details{
If \code{group} is supplied, the function returns
  a list (see value section).  Otherwise, the function returns a numeric
  vector containing the spherical score associated with each forecast.
}
\value{
Depending on input arguments, \code{sphscore} may return an object of
  class \code{numeric} containing raw spherical scores.  It may also return
  a list containing the objects below.
  
\item{rawscores}{an object of class \code{numeric} containing raw spherical
  scores for each forecast.}
\item{mnsph}{an object of class \code{numeric} containing mean spherical
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
## Raw spherical scores
sph1 <- sphscore(answer ~ forecast, data=WorldEvents)
## Raw spherical scores plus group means
sph2 <- sphscore(answer ~ forecast, data=WorldEvents,
                 group="forecaster")
}
