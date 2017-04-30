#' Select optimal degrees of freedom for SL.gam() using cross-validation
#'
#' ab_cvGAM is an internal tuning function called by \code{agecurveAb} and \code{tmleAb} that selects degrees of freedom for natural splines in a GAM model using cross-validation
#'
#' @param Y The outcome. Must be a numeric vector.
#' @param X A matrix of features that predict Y, usually a data.frame.
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Model family (gaussian for continuous outcomes, binomial for binary outcomes)
#' @param SL.library SuperLearner library
#' @param cvControl Optional list to control cross-valiation (see \code{\link[SuperLearner]{SuperLearner}} for details).
#' @param print logical. print messages? Defaults to FALSE
#' @param df a sequence of degrees of freedom to control the smoothness of natural splines in the GAM model. Defaults to 2:6
#'
#' @return returns a list with updated SuperLearner library, the optimal node size, and cvRisks
#' @details \code{ab_cvGAM} is an internal function called by \code{\link[tmleAb]{agecurveAb}} or \code{\link[tmleAb]{tmleAb}} if SL.gam() is included in the algorithm library. It performs an addition pre-screen step of selecting the optimal spline degress of freedom using cross validation. The default is to search over degrees 2,3,...10, which is usually pretty good. This additional selection step enables you to tune the smoothing parameter. Cross-validated risks are estimated using \code{\link[SuperLearner]{SuperLearner}}.
#' @examples
#' # TBD
#' @keywords internal
#' @export

ab_cvGAM <- function(Y,X,id=NULL,family=gaussian(),SL.library,cvControl=list(),print=FALSE, df=2:10) {
  if(is.null(df)) {
    df <- 2:10
  }
  if(print==TRUE) {
    cat("\nThe ensemble library includes SL.gam.")
    cat("\nThe default R implementation of gam() may over- or under-smooth the data")
    cat("\nTuning the fit by selecting the optimal df for the smoothing splines")
    cat("\nfrom ",df[1]," to ",df[length(df)]," using V-fold cross-validation.")
  }

  if(is.null(id)) id <- 1:length(Y)


  create.SL.gam <- function(tune = list(df = df)) {
    for(mm in seq(length(tune$df))) {
      eval(parse(file = "", text = paste("SL.gam.df", tune$df[mm], "<- function(...,deg.gam = ", tune$df[mm], ") SL.gam(..., deg.gam = deg.gam)", sep = "")), envir = .GlobalEnv)
    }
    invisible(TRUE)
  }
  create.SL.gam()

  # estimate cross-validated Risks for different node sizes
  cvRisks <- rep(NA,length(df))
  for(nn in seq(length(df))) {
    fit <- SuperLearner(Y=Y,X=X,id=id,family=family,SL.library=paste("SL.gam.df",df[nn],sep=""),cvControl=cvControl)
    cvRisks[nn] <- fit$cvRisk
  }
  # identify the lowest risk
  df_opt <- df[order(cvRisks)][1]
  cvr_tab <- cbind(df,cvRisks)
  colnames(cvr_tab) <- c("df","CVRisk")

  # update the library
  SLlib2 <- gsub("SL.gam",paste("SL.gam.df",df_opt,sep=""),SL.library)

  if(print==TRUE) {
    cat("\n-----------------------------------")
    cat("\nGeneralized additive model with natural splines\n")
    cat("\nOptimal smoothing df: ",df_opt,"\n")
    print(cvr_tab)
    cat("-----------------------------------\n")
  }

  # return results
  return(list(SL.library=SLlib2,df_opt=df_opt,cvRisks=cvr_tab))

}

