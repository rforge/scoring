\name{calcscore}
\alias{calcscore}
\alias{calcscore.formula}
\alias{calcscore.default}

\title{
Calculate Scores Under A Specific Rule
}

\description{
Given parameters of a scoring rule family, calculate scores for probabilistic forecasts and associated outcomes.
}
\usage{
\method{calcscore}{formula}(object, fam="pow", param, data, bounds=NULL,
                  reverse=FALSE, ...)

\method{calcscore}{default}(object, outcome, fam="pow",
                  param=c(2,rep(1/max(2,NCOL(forecast)),max(2,NCOL(forecast)))),
                  bounds=NULL, reverse=FALSE, ...)
}
\arguments{
  \item{object}{an object of class  "formula", of the form
    \code{outcome ~ forecast} (see details).  Alternatively, a matrix of forecasts,
          with observations in rows and
          forecast alternatives in columns.  For two-alternative forecasts, this
          can be a vector reflecting forecasts for one alternative.}
  \item{outcome}{a vector of outcomes.  For each row of the
          \code{forecast} matrix, \code{outcome} should contain an entry
          reflecting the column number associated with the event that occurred.}	
  \item{fam}{scoring rule family. \code{pow} (default) is the power
    family, \code{beta} is the beta family, \code{sph} is the
    pseudospherical family.}
  \item{param}{for family \code{beta}, a numeric vector of length 2
    containing the scoring rule family parameters.
    For other families, a numeric vector of length
    \code{NCOL(forecast)+1}, containing first the family parameter gamma
    and then the baseline parameters (see details).}
  \item{data}{an optional data frame or list containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from the environment from which
	  \code{calcscore} is called.}
  \item{bounds}{a vector of length 2 corresponding to the desired
    minimum value and maximum value of the scoring rule, respectively.  Entries of \code{NA}
  imply that the minimum and/or maximum bound will not be modified from
                  the natural, family-implied bounds.}
  \item{reverse}{if \code{FALSE} (default), smaller scores imply better
    forecasts.  If \code{TRUE}, larger scores imply better forecasts.}
  \item{...}{Additional arguments.}
}
\details{
  The formula is of the form \code{outcome ~ forecast}, where \code{forecast} describes the column(s) containing forecasts associated with the possible outcomes.  Multiple columns are separated by \code{+}.  \code{outcome} is always a vector describing the outcome associated with each forecast.  These should be coded 1, 2, ..., reflecting the column associated with the outcome (see examples).

For events with only two alternatives, one can take a shortcut and supply only forecasts and baseline parameters (for families \code{pow} and \code{sph}) associated with a single outcome.  In this case, the \code{outcome} vector should contain zeros and ones, where `one' means that the forecasted alternative occurred.

When the \code{bounds} argument is supplied, the code attempts to scale the scores so that the maximum score is \code{bounds[2]} and the minimum score is \code{bounds[1]}.  This scaling cannot be accomplished when the scoring rule allows scores of infinity (the log score is the most common case here).  If \code{reverse=TRUE}, the bounds are applied after the reversal (so that the supplied lower bound reflects the worst score and upper bound reflects the best score).
}

\note{  The beta family was originally proposed by Buja et al. (2005); the power and pseudospherical families with baseline are described by Jose et al. (2009).  A discussion of choosing specific rules from these families is provided by Merkle and Steyvers (2013).

Some notable special cases of these families are:

Beta family: Log score when parameters are (0,0); Brier score when parameters are (1,1).

Power family with baseline parameters all equal (to 1/(number of alternatives)): The family approaches the log score as gamma goes to 1 (but the family is undefined for gamma=1).  One obtains the Brier score for gamma=2.

Pseudospherical family with baseline parameters all equal: The family approaches the log score as gamma goes to 1 (but the family is undefined for gamma=1).  One obtains the spherical score for gamma=2.

}
\value{
\code{calcscore} returns a numeric vector that has length equal to 
\code{length(outcome)}, containing scores under the selected scoring rule.
}
\references{
Buja, A., Stuetzle, W., & Shen, Y. (2005). Loss functions for binary
class probability estimation and classification: Structure and
applications. (Obtained from \url{http://stat.wharton.upenn.edu/~buja/PAPERS/})

Jose, V. R. R., Nau, R. F., & Winkler, R. L. (2008). Scoring rules, generalized entropy, and utility maximization. \emph{Operations Research, 56}, 1146--1157.

Jose, V. R. R., Nau, R. F., & Winkler, R. L. (2009). Sensitivity to distance and baseline distributions in forecast evaluation. \emph{Management Science, 55}, 582--590.

Merkle, E. C. & Steyvers, M. (in press). Choosing a strictly proper scoring rule. \emph{Decision Analysis}.
}
\author{
Ed Merkle
}

\seealso{\code{\link{plotscore}}}
\examples{
## Brier scores for two alternatives, with bounds of 0 and 1
data("WorldEvents")
scores <- calcscore(answer ~ forecast, fam="beta",
                    param=c(1,1), data=WorldEvents,
                    bounds=c(0,1))

## Calculate Brier scores manually
scores.man <- with(WorldEvents, (forecast - answer)^2)

## Comparison
all.equal(scores, scores.man)

## Average Brier score for each forecaster
with(WorldEvents, tapply(scores, forecaster, mean))

## Brier scores for 3 alternatives, with bounds of 0 and 1
data("WeatherProbs")
scores2 <- calcscore(tcat ~ tblw + tnrm + tabv, fam="pow",
                     param=c(2,rep(1/3,3)), data=WeatherProbs,
                     bounds=c(0,1))

## Spherical scores for 3 alternatives, reversed so 0 is worst and
## 1 is best
scores3 <- calcscore(tcat ~ tblw + tnrm + tabv, fam="sph",
                     param=c(2,rep(1/3,3)), data=WeatherProbs,
                     bounds=c(0,1), reverse=TRUE)
}
