#' Selects optimal tree depth for randomForest() using corss-validation by tuning the nodesize parameter
#'
#' @param Y outcome
#' @param X matrix of features that predict Y
#' @param id id variable to identify independent observations for V-splits
#' @param SLlib SuperLearner library
#' @param print logical. print messages? Defaults to FALSE
#'
#' @return returns a list with updated SuperLearner library, the optimal node size, and cvRisks
#' @examples
#'
#Additional Context can be found at bottom of file

SLAb.cvRF <- function(Y,X,id,SLlib,print=FALSE) {
  if(print==TRUE) {
    cat("\nThe ensemble library includes SL.randomForest.")
    cat("\nThe default R implementation of randomForest tends to overfit the data")
    cat("\n  by growing trees that are too deep.")
    cat("\nPrevent overfitting by properly tuning the node size.")
    cat("\nSelecting the optimal node size (tree depth)")
    cat("\n  from 15,20,...,40 using V-fold cross-validation.")
    cat("\nThis could take a few minutes, depending on the size of your dataset...\n")
  }

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
  SLlib2 <- gsub("SL.randomForest",paste("SL.randomForest.ns",ns.opt,sep=""),SLlib)

  if(print==TRUE) {
    cat("\n-----------------------------------")
    cat("\nOptimal node size: ",ns.opt,"\n")
    print(cvr.tab)
    cat("-----------------------------------\n")
  }

  # return results
  return(list(SLlib=SLlib2,ns.opt=ns.opt,cvRisks=cvr.tab))

}


#---------------------------------------
# SLAb-cvRF.R
#
# Developed by Ben Arnold (benarnold@berkeley.edu)
#
# select optimal tree depth for
# randomForest() using cross-validation
# by tuning the nodesize parameter
#
# if you don't do this, RF will usually
# grow the trees that are too deep
# (e.g., to nodes with 1 obs)
# and this will overfit the data
#
# In the context of Age-antibody curves,
# it will make extremely jagged curves
# that are clear overfits. This adds
# the one-step optimization of that parameter
# that should probably be included in the
# default implementation of the
# randomForest package, but isn't
#
# this function is designed to be
# called by SLAb.curve or SLAb.tmle
# if they include SL.randomForest in the
# ensemble library
# it leverages SuperLearner() to estimate
# cross-validated Risks
#---------------------------------------
