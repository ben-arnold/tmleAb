% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ab_cvRF.R
\name{ab_cvRF}
\alias{ab_cvRF}
\title{Select optimal tree depth for SL.randomForest() using cross-validation}
\usage{
ab_cvRF(Y, X, id = NULL, family = gaussian(), SL.library,
  cvControl = list(), print = FALSE, RFnodesize = seq(15, 40, by = 5))
}
\arguments{
\item{Y}{The outcome. Must be a numeric vector.}

\item{X}{A matrix of features that predict Y, usually a data.frame.}

\item{id}{An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.}

\item{family}{Model family (gaussian for continuous outcomes, binomial for binary outcomes)}

\item{SL.library}{SuperLearner library}

\item{cvControl}{Optional list to control cross-valiation (see \code{\link[SuperLearner]{SuperLearner}} for details).}

\item{print}{logical. print messages? Defaults to FALSE}

\item{RFnodesize}{a sequence of nodes used by the random forest algorithm. Defaults to a sequence from 15 to 40 by every 5 nodes}
}
\value{
returns a list with updated SuperLearner library, the optimal node size, and cvRisks (cross-validated risks for each nodesize evaluated)
}
\description{
ab_cvRF is an internal tuning function called by \code{agecurveAb} and \code{tmleAb} that selects optimal tree depth by tuning the nodesize parameter
}
\details{
\code{ab_cvRF} is an internal function called by \code{\link[tmleAb]{agecurveAb}} or \code{\link[tmleAb]{tmleAb}} if SL.randomForest() is included in the algorithm library. It performs an addition pre-screen step of selecting the optimal node depth for random forest using cross validation. The default range of node sizes evaluated is 15, 20, ..., 40.  In the context of Age-antibody curves, without this tuning step random forest will fit extremely jagged curves that are clear overfits. This additional selection step prevents overfitting. Cross-validated risks are estimated using \code{\link[SuperLearner]{SuperLearner}}.
}
\examples{
# TBD
}
\keyword{internal}
