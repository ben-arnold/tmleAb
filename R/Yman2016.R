



#' Antibody acquisition models
#'
#' Fit a parametric model of age-dependent antibody acquisition (and loss)
#'
#' @param logAb Log antibody level (numeric)
#' @param Age Age of the individual at the time of measurement (numeric)
#' @param type Specify type of model to fit -- either \code{constant}, \code{changepoint}, or \code{linear}. At this time only constant is supported and this argument is currently ignored.
#' @param print Logical. Print parameter estimates?
#'
#' @return Returns an object of class \code{Yman2016} that is just the maximum likelihood fit results passed from optim.
#' @details \code{Yman2016} estimates rates of antibody acquisision (alpha) and loss (r) using maximum likelihood under a parametric model that assumes constant rates of acquisition and loss, and that antibody levels follow a log-normal distribution. The function also estimates nuissance parameter (sigma, the standard deviation of the log-normal distribution). Future updates may include the changepoint and linear decrease models from the original paper
#' @references Yman V, White MT, Rono J, Arcà B, Osier FH, Troye-Blomberg M, et al. Antibody acquisition models: A new tool for serological surveillance of malaria transmission intensity. Sci Rep. ; 2016;6: 19472. (\url{http://www.nature.com/articles/srep19472})
#' @export
#'
#' @examples
#' \dontrun{
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
#' yman2016fit <- Yman2016(logAb=Y,Age=X$Age)
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
#' }
Yman2016 <- function(logAb,Age,type="constant",print=TRUE) {

  # objective function
  LL <- function(theta, Age, logAb) {
    # theta : parameters for optimization, length 3, including:
      # alpha : antibody acquisition rate
      # r     : antibody loss rate
      # sigma : sd of the log-normal distribution
    # Age     : age of individual i
    # logAb   : log antibody level of individual i

    if(is.numeric(logAb)==FALSE){
      stop("logAb must be numeric")
    }
    if(is.numeric(Age)==FALSE){
      stop("Age must be numeric")
    }

    alpha <- theta[1]
    r     <- theta[2]
    sigma <- theta[3]

    # get age-specific predicted geomean antibody level given alpha, r
    Aa <- as.numeric((alpha/r)*(1-exp(-r*Age)))

    # likelihood function for log Ab levels
    L <- dnorm(logAb,mean=log(Aa),sd=sigma)

    # negative log likelihood function (to minimize with optim)
    -sum(log(L))
  }

  # ML fit of L(theta | X, Y)
  mlfit <- optim(c(0.1,0.01,1),fn=LL,Age=Age,logAb=logAb)
  class(mlfit) <- "Yman2016"

  if(print==TRUE) {
    p_res <- matrix(mlfit$par,nrow=3,ncol=1)
    rownames(p_res) <- c("Ab acquisition rate","Ab loss rate","Sigma")
    colnames(p_res) <- c("ML Estimate")
    cat("\n-------------------------------------------\nYman et al. 2016 constant rate model\nMaximum likelihood parameter estimates:\n\n")
    print(p_res)
    cat("-------------------------------------------\n")
  }


  return(mlfit)

}


#' Internal prediction method for Yman2016
#'
#' @param object An object of class \code{Yman2016}
#' @param newdata New data used for predicting the outcome (must be univariate and equal to age, or if multivariate, it must include a variable 'Age')
#' @param ... Other arguments passed to the function (currently ignored)
#'
#' @return pred Predicted values given \code{newdata}
#' @keywords internal
#' @export
#'
predict.Yman2016 <- function(object,newdata,...) {

  if( is.null(dim(newdata)) ){
    age <- as.numeric(newdata)
  }
  else {
    if(length(grep("Age",names(newdata)))<=0 ) {
      error("For the Yman2016 model, newdata must include a variable named 'Age'.")
    }
    age <- as.numeric(newdata[,"Age"])
  }

  pred <-  log( (object$par[1]/object$par[2])*(1-exp(-object$par[2]*age)) )
  pred
}

