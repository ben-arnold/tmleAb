#' Selects optimal tree depth for SL.randomForest() using cross-validation by tuning the nodesize parameter
#'
#' @param Y The outcome. Must be a numeric vector.
#' @param X A matrix of features that predict Y, usually a data.frame.
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param SL.library SuperLearner library
#' @param print logical. print messages? Defaults to FALSE
#'
#' @return returns a list with updated SuperLearner library, the optimal node size, and cvRisks
#' @details \code{slab_cvRF} is an internal function called by \code{\link[SLAb]{slab_curve}} or \code{slab_tmle} if SL.randomForest() is included in the algorithm library. It performs an addition pre-screen step of selecting the optimal node depth for random forest using cross validation. The default range of node sizes evaluated is 15, 20, ..., 40.  In the context of Age-antibody curves, without this tuning step random forest will fit extremely jagged curves that are clear overfits. This additional selection step prevents overfitting. Cross-validated risks are estimated using \code{\link[SuperLearner]{SuperLearner}}.
#' @examples TBD


slab_cvRF <- function(Y,X,id=NULL,SL.library,print=FALSE) {
  if(print==TRUE) {
    cat("\nThe ensemble library includes SL.randomForest.")
    cat("\nThe default R implementation of randomForest tends to overfit the data")
    cat("\n  by growing trees that are too deep.")
    cat("\nPrevent overfitting by properly tuning the node size.")
    cat("\nSelecting the optimal node size (tree depth)")
    cat("\n  from 15,20,...,40 using V-fold cross-validation.")
    cat("\nThis could take a few minutes, depending on the size of your dataset...\n")
  }

  if(is.null(id)) id <- 1:length(Y)

  nodesizes <- seq(15,40,by=5)
  create.SL.randomForest <- function(tune = list(nodesize = nodesizes)) {
    for(mm in seq(length(tune$nodesize))) {
      eval(parse(file = "", text = paste("SL.randomForest.ns", tune$nodesize[mm], "<- function(...,nodesize = ", tune$nodesize[mm], ") SL.randomForest(..., nodesize = nodesize)", sep = "")), envir = .GlobalEnv)
    }
    invisible(TRUE)
  }
  create.SL.randomForest()

  # estimate cross-validated Risks for different node sizes
  cvRisks <- rep(NA,length(nodesizes))
  for(nn in seq(length(nodesizes))) {
    fit <- SuperLearner(Y=Y,X=X,id=id,SL.library=paste("SL.randomForest.ns",nodesizes[nn],sep=""))
    cvRisks[nn] <- fit$cvRisk
  }
  # identify the lowest risk
  ns.opt <- nodesizes[order(cvRisks)][1]
  cvr.tab <- cbind(nodesizes,cvRisks)
  colnames(cvr.tab) <- c("nodesize","CVRisk")

  # update the library
  SLlib2 <- gsub("SL.randomForest",paste("SL.randomForest.ns",ns.opt,sep=""),SL.library)

  if(print==TRUE) {
    cat("\n-----------------------------------")
    cat("\nOptimal node size: ",ns.opt,"\n")
    print(cvr.tab)
    cat("-----------------------------------\n")
  }

  # return results
  return(list(SL.library=SLlib2,ns_opt=ns.opt,cvRisks=cvr.tab))

}

