\name{plotscore}
\alias{plotscore}
\title{
Plot a Two-Alternative Scoring Rule
}
\description{
Given parameters for a two-alternative scoring rule, plot the hypothetical scores that would be obtained for each forecast/outcome combination.
}
\usage{
plotscore(param = c(2, 0.5), fam = "pow", bounds, reverse = FALSE,
          legend = TRUE, ...)
}
\arguments{
  \item{param}{Numeric vector of length 2, containing the parameters for
    \code{fam}.  For family \code{beta}, these are the parameters
    commonly denoted alpha and beta.  For families \code{pow} and
    \code{sph}, these correspond to the family parameter gamma and the
    baseline parameter associated with the focal outcome, respectively.}
  \item{fam}{scoring rule family. \code{pow} (default) is the power
    family, \code{beta} is the beta family, \code{sph} is the
    pseudospherical family.}
  \item{bounds}{Lower and upper bounds supplied to \code{calcscore}.}
  \item{reverse}{\code{reverse} argument supplied to \code{calcscore}.}
  \item{legend}{Should a legend be displayed?  Defaults to \code{TRUE}}
  \item{...}{Other arguments to \code{plot()}}
}
\details{
For more information on the scoring rule families and the \code{bounds} and \code{reverse} arguments, see the details of \code{calcscore()}.
}
\value{
Returns the result of a \code{plot()} call that graphs the scoring rule.
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

\seealso{
\code{\link{calcscore}}
}
\examples{
## Plot Brier score from power family with natural bounds
plotscore(c(2,.5), fam="pow")

## Plot Brier score from beta family with bounds of 0 and 1
plotscore(c(1,1), fam="beta", bounds=c(0,1))

## Plot log score
plotscore(c(0,0), fam="beta")

## Score from pseudospherical family with
## baseline of .3 and (0,1) bounds
plotscore(c(3, .3), fam="sph", bounds=c(0,1))
}
