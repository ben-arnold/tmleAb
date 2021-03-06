% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL.Yman2016.R
\name{SL.Yman2016}
\alias{SL.Yman2016}
\title{SuperLearner wrapper for antibody acquisition models}
\usage{
SL.Yman2016(Y, X, newX = NULL, ...)
}
\arguments{
\item{Y}{log quantitative antibody level}

\item{X}{age}

\item{newX}{(optional) new values of age over which to predict}

\item{...}{other arguments passed to the function (currently ignored)}
}
\value{
\code{pred} Predicted outcome values

\code{fit} Maximum likelihood fit object returned from optim
}
\description{
An internal SuperLearner wrapper for the antibody acquisition and loss model in Yman et al. 2016
}
\details{
The \code{SL.Yman2016} wrapper enables the inclusion of this model in a \code{\link[SuperLearner]{SuperLearner}} ensemble. The model assumes that rates of antibody acquistion and loss are constant. It can only use information about age and not additional covariates. If \code{X} includes more than 1 column, then the function will assume that the first column of \code{X} is age and will throw a warning.
}
\seealso{
\code{\link[tmleAb]{Yman2016}}

\code{\link[SuperLearner]{SuperLearner}}
}
\keyword{internal}
