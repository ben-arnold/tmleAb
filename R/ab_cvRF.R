#' Select optimal tree depth for SL.randomForest() using cross-validation
#'
#' ab_cvRF is an internal tuning function called by \code{agecurveAb} and \code{tmleAb} that selects optimal tree depth by tuning the nodesize parameter
#'
#' @param Y The outcome. Must be a numeric vector.
#' @param X A matrix of features that predict Y, usually a data.frame.
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Model family (gaussian for continuous outcomes, binomial for binary outcomes)
#' @param SL.library SuperLearner library
#' @param print logical. print messages? Defaults to FALSE
#' @param RFnodesize a sequence of nodes used by the random forest algorithm. Defaults to a sequence from 15 to 40 by every 5 nodes
#'
#' @return returns a list with updated SuperLearner library, the optimal node size, and cvRisks (cross-validated risks for each nodesize evaluated)
#' @details \code{ab_cvRF} is an internal function called by \code{\link[tmleAb]{agecurveAb}} or \code{\link[tmleAb]{tmleAb}} if SL.randomForest() is included in the algorithm library. It performs an addition pre-screen step of selecting the optimal node depth for random forest using cross validation. The default range of node sizes evaluated is 15, 20, ..., 40.  In the context of Age-antibody curves, without this tuning step random forest will fit extremely jagged curves that are clear overfits. This additional selection step prevents overfitting. Cross-validated risks are estimated using \code{\link[SuperLearner]{SuperLearner}}.
#' @examples
#' # TBD
#' @keywords internal
#' @export

ab_cvRF <- function(Y,X,id=NULL,family=gaussian(),SL.library,print=FALSE, RFnodesize=seq(15,40,by=5)) {
  if(is.null(RFnodesize)){
    RFnodesize <- seq(15,40,by=5)
  }
  if(print==TRUE) {
    cat("\nThe ensemble library includes SL.randomForest.")
    cat("\nThe default R implementation of randomForest tends to overfit the data")
    cat("\n\tby growing trees that are too deep.")
    cat("\nPrevent overfitting by properly tuning the node size.")
    cat("\nSelecting the optimal node size (tree depth)")
    cat("\n",RFnodesize[1],"to",RFnodesize[length(RFnodesize)],"using V-fold cross-validation.")
    cat("\nThis could take a few minutes, depending on the size of your dataset...\n")
  }

  if(is.null(id)) id <- 1:length(Y)


  create.SL.randomForest <- function(tune = list(nodesize = RFnodesize)) {
    for(mm in seq(length(tune$nodesize))) {
      eval(parse(file = "", text = paste("SL.randomForest.ns", tune$nodesize[mm], "<- function(...,nodesize = ", tune$nodesize[mm], ") SL.randomForest(..., nodesize = nodesize)", sep = "")), envir = .GlobalEnv)
    }
    invisible(TRUE)
  }
  create.SL.randomForest()

  # estimate cross-validated Risks for different node sizes
  cvRisks <- rep(NA,length(RFnodesize))
  for(nn in seq(length(RFnodesize))) {
    fit <- SuperLearner(Y=Y,X=X,id=id,family=family,SL.library=paste("SL.randomForest.ns",RFnodesize[nn],sep=""))
    cvRisks[nn] <- fit$cvRisk
  }
  # identify the lowest risk
  ns_opt <- RFnodesize[order(cvRisks)][1]
  cvr_tab <- cbind(RFnodesize,cvRisks)
  colnames(cvr_tab) <- c("nodesize","CVRisk")

  # update the library
  SLlib2 <- gsub("SL.randomForest",paste("SL.randomForest.ns",ns_opt,sep=""),SL.library)

  if(print==TRUE) {
    cat("\n-----------------------------------")
    cat("\nRandom Forest\n")
    cat("\nOptimal minimum node size: ",ns_opt,"\n")
    print(cvr_tab)
    cat("-----------------------------------\n")
  }

  # return results
  return(list(SL.library=SLlib2,ns_opt=ns_opt,cvRisks=cvr_tab))

}

