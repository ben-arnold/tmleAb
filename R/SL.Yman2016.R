
#' Antibody acquisition model
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
#' @export
#' @references  Yman V, White MT, Rono J, Arcà B, Osier FH, Troye-Blomberg M, et al. Antibody acquisition models: A new tool for serological surveillance of malaria transmission intensity. Sci Rep. ; 2016;6: 19472. (\url{http://www.nature.com/articles/srep19472})
#'
#' @examples
#' # simulate data assuming the Yman 2016 model is the truth
#' set.seed(1234)
#' age <- runif(1000,min=1,max=20)
#' alpha <- 0.5
#' r <- 0.07
#' sigma <- log(1.7)
#' abtrue <- (alpha/r)*(1-exp(-r*age))
#' abobs <- exp( log(abtrue) +rnorm(length(age),mean=0,sd=sigma) )
#'
#' # Fit the model and obtain predictions
#' X <- data.frame(age=age)
#' Y <- log(abobs)
#' yman2016fit <- SL.Yman2016(Y=Y,X=X)
#' pY <- predict(yman2016fit$fit,newdata=X)
#'
#' # compare against truth
#' cbind(c(alpha,r,sigma),yman2016fit$fit$object$par)
#' # figure (not run)
#' # plot(age,log(abobs),cex=0.25,col="gray40",pch=16)
#' # lines(age[order(age)],log(abtrue[order(age)]),lwd=3,lty=2)
#' # lines(age[order(age)],pY[order(age)],lwd=1)
#'
#' # super learner fit
#' # (note that it selects Yman 2016 under
#' # conditions where it models the truth)
#' SL.library <- c("SL.glm","SL.loess","SL.Yman2016")
#' SLfit <- SuperLearner(Y=log(abobs),X=data.frame(age=age),SL.library=SL.library)
#' SLfit
#'
#' # The Yman 2016 model can only use age, not additional covariates,
#' # so it will throw a warning if X includes more than 1 covariate
#' # and it will assume that the first column of X is age
#' SLfit2 <- SuperLearner(Y=log(abobs),X=data.frame(age=age,age2=age^2),SL.library=SL.library)
#' SLfit2
#'
SL.Yman2016 <- function(Y,X,newX=NULL,...) {

  # objective function
  LL <- function(theta, Age, logAb) {
    # theta : parameters for optimization, length 3, including:
    # alpha : antibody acquisition rate
    # r     : antibody loss rate
    # sigma : sd of the log-normal distribution
    # X     : age of individual i
    # Y     : log antibody level of individual i

    alpha <- theta[1]
    r     <- theta[2]
    sigma <- theta[3]

    # get age-specific predicted geomean antibody level given alpha, r
    Aa <- (alpha/r)*(1-exp(-r*Age))

    # likelihood function for log Ab levels
    L <- dnorm(logAb,mean=log(Aa),sd=sigma)

    # negative log likelihood function (to minimize with optim)
    -sum(log(L))
  }

  if(is.null(newX)) {
    newX <- X
  }

  # Stop if X does not include Age
  if(length(grep("Age",names(X)))<=0 ) {
    error("For the SL.Yman2016 model, X must include a variable named 'Age'.")
  }
  X <- X[,"Age"]

  if(length(grep("Age",names(newX)))<=0 ) {
    error("For the SL.Yman2016 model, newX must include a variable named 'Age'.")
  }
  newX <- newX[,"Age"]

  # ML fit of L(theta | X, Y)
  mlfit <- optim(c(0.1,0.01,1),fn=LL,Age=X,logAb=Y)
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
#' @param object An object of class \code{SL.Yman2016}
#' @param newdata New data used for predicting the outcome (must be univariate and equal to age, or if multivariate, it must include a variable 'Age')
#' @param ... Other arguments passed to the function (currently ignored)
#'
#' @return pred Predicted values given \code{newdata}
#' @export
#'
#' @examples
predict.SL.Yman2016 <- function(object,newdata,...) {

  if( is.null(dim(newdata)) ){
    age <- as.numeric(newdata)
  }
  else {
    if(length(grep("Age",names(newdata)))<=0 ) {
      error("For the SL.Yman2016 model, newX must include a variable named 'Age'.")
    }
    age <- as.numeric(newdata[,"Age"])
  }

  pred <-  log( (object$object$par[1]/object$object$par[2])*(1-exp(-object$object$par[2]*age)) )
  pred
}


