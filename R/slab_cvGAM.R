#' Selects optimal degrees of freedom for SL.gam() using cross-validation by tuning the df parameter
#'
#' @param Y The outcome. Must be a numeric vector.
#' @param X A matrix of features that predict Y, usually a data.frame.
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Model family (gaussian for continuous outcomes, binomial for binary outcomes)
#' @param SL.library SuperLearner library
#' @param print logical. print messages? Defaults to FALSE
#' @param df a sequence of degrees of freedom to control the smoothness of natural splines in the GAM model. Defaults to a sequence from 2 to 10
#'
#' @return returns a list with updated SuperLearner library, the optimal node size, and cvRisks
#' @details \code{slab_cvGAM} is an internal function called by \code{\link[SLAb]{slab_curve}} or \code{slab_tmle} if SL.gam() is included in the algorithm library. It performs an addition pre-screen step of selecting the optimal node depth for random forest using cross validation. The default range of degrees of freedom is 2, 3, ...10, which is usually sufficient.  In the context of Age-antibody curves, without this tuning step the default parameters for GAM will usually under-smooth. This additional selection step tunes the smoothing parameter. Cross-validated risks are estimated using \code{\link[SuperLearner]{SuperLearner}}.
#' @examples TBD
#' @export

slab_cvGAM <- function(Y,X,id=NULL,family=gaussian(),SL.library,print=FALSE, df=2:10) {
  if(print==TRUE) {
    cat("\nThe ensemble library includes SL.gam.")
    cat("\nThe default R implementation of gam() may over- or under-smooth the data")
    cat("\nTuning the fit by selecting the optimal df for the smoothing splines")
    cat("\nfrom ",df[1]," to ",df[length(df)]," using V-fold cross-validation.")
  }

  if(is.null(id)) id <- 1:length(Y)


  create.SL.gam <- function(tune = list(df = df)) {
    for(mm in seq(length(tune$df))) {
      eval(parse(file = "", text = paste("SL.gam.df", tune$df[mm], "<- function(...,df = ", tune$df[mm], ") SL.gam(..., df = df)", sep = "")), envir = .GlobalEnv)
    }
    invisible(TRUE)
  }
  create.SL.gam()

  # estimate cross-validated Risks for different node sizes
  cvRisks <- rep(NA,length(df))
  for(nn in seq(length(df))) {
    fit <- SuperLearner(Y=Y,X=X,id=id,family=family,SL.library=paste("SL.gam.df",df[nn],sep=""))
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
    cat("\nOptimal smoothing df: ",df_opt,"\n")
    print(cvr_tab)
    cat("-----------------------------------\n")
  }

  # return results
  return(list(SL.library=SLlib2,df_opt=df_opt,cvRisks=cvr_tab))

}

