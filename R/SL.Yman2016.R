
#' SuperLearner wrapper for antibody acquisition models
#'
#' An internal SuperLearner wrapper for the antibody acquisition and loss model in Yman et al. 2016
#'
#' @param Y log quantitative antibody level
#' @param X age
#' @param newX (optional) new values of age over which to predict
#' @param ... other arguments passed to the function (currently ignored)
#'
#' @details The \code{SL.Yman2016} wrapper enables the inclusion of this model in a \code{\link[SuperLearner]{SuperLearner}} ensemble. The model assumes that rates of antibody acquistion and loss are constant. It can only use information about age and not additional covariates. If \code{X} includes more than 1 column, then the function will assume that the first column of \code{X} is age and will throw a warning.
#' @return \code{pred} Predicted outcome values
#' @return \code{fit} Maximum likelihood fit object returned from optim
#'
#'
#' @seealso \code{\link[tmleAb]{Yman2016}}
#' @seealso \code{\link[SuperLearner]{SuperLearner}}
#'
#' @keywords internal
#' @export

SL.Yman2016 <- function(Y,X,newX=NULL,...) {

  if(is.null(newX)) {
    newX <- X
  }

  # Stop if X does not include Age
  if(length(grep("Age",names(X)))<=0 ) {
    stop("For the Yman2016 model, X must include a variable named 'Age'.")
  }
  X <- X[,"Age"]

  if(length(grep("Age",names(newX)))<=0 ) {
    stop("For the Yman2016 model, newX must include a variable named 'Age'.")
  }
  newX <- newX[,"Age"]

  # ML fit of L(theta | X, Y)
  mlfit <- Yman2016(logAb=Y,Age=X,print=FALSE)
  fit <- list(object = mlfit)
  class(fit) <- "SL.Yman2016"

  # predicted values
  pred <- predict(fit,newdata=newX)

  # return results
  out <- list(pred=pred, fit=fit)
  return(out)
}


#' Internal prediction method for SL.Yman2016
#'
#' @param object An object of class \code{Yman2016}
#' @param newdata New data used for predicting the outcome (must be univariate and equal to age, or if multivariate, it must include a variable 'Age')
#' @param ... Other arguments passed to the function (currently ignored)
#'
#' @return pred Predicted values given \code{newdata}
#'
#' @keywords internal
#' @export
#'
predict.SL.Yman2016 <- function(object,newdata,...) {

  if( is.null(dim(newdata)) ){
    age <- as.numeric(newdata)
  }
  else {
    if(length(grep("Age",names(newdata)))<=0 ) {
      error("For the Yman2016 model, newX must include a variable named 'Age'.")
    }
    age <- as.numeric(newdata[,"Age"])
  }

  pred <-  log( (object$object$par[1]/object$object$par[2])*(1-exp(-object$object$par[2]*age)) )
  pred
}

