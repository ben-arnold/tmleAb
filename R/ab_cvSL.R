
#' Wrapper function to estimate cross-validated risk for the super learner ensemble using antibody data
#'
#' A convenience wrapper for \code{\link[SuperLearner]{CV.SuperLearner}}
#'
#' @param Y Antibody measurement. Must be a numeric vector.
#' @param Age Age of the individual at the time of measurement. Must be a numeric vector.
#' @param W An optional vector, matrix, or data.frame of covariates for each individual used to marginally adjust the curve
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Model family (gaussian for continuous outcomes, binomial for binary outcomes)
#' @param SL.library Library of algorithms to include in the ensemble (see the \code{\link[SuperLearner]{SuperLearner}} package for details).
#' @param RFnodesize Optional argument to specify a range of minimum node sizes for the random Forest algorithm. If \code{SL.library} includes \code{SL.randomForest}, then the default is to search over node sizes of 15,20,...40. Specifying this option will override the default.
#' @param gamdf Optional argument to specify a range of degrees of freedom for natural smoothing splines in a generalized additive model. If \code{SL.library} includes \code{SL.gam}, then the default is to search over degrees of freedom 2,3,...10. Specifying this option will override the default.
#'
#' @details TBD
#'
#' @return
#' @export
#'
#' @examples TBD
#'
ab_cvSL <-function(Y,Age,W=NULL,id=1:length(Y),family=gaussian(),V=10,SL.library= c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.randomForest"),RFnodesize=NULL,gamdf=NULL) {
	# wrapper function for the CV.SuperLearner algorithm
  # to compute cross-validated risk for the super learner and its constituent algorithms
	#
	# Y   : outcome_i (log10 antibody)
	# Age : age_i
  # W    : matrix of covariates (optional)
	# id  : ID variable for individuals in the dataset (for repeated observations)
  # family : gaussian or binomial
  # V : number of splits for V-fold cross validation
	# SL.Library : SuperLearner algorithm library (if different from the default, above)
	#
	# the function returns a data frame with id, Y, Age, W, and pY (SL predictions of marginally averaged Y at A=a)

	require(SuperLearner)


  # convert W into a design matrix (SuperLearner does not acommodate factor variables)
  if (is.null(W)) {
    fitd <- data.frame(id,Y,Age)
	} else{
	  Wdesign <- design_matrix(W)
	  fitd <- data.frame(id,Y,Age,Wdesign)
	}

	# restrict dataset to non-missing observations
	fitd <- fitd[complete.cases(fitd),]
	X <- subset(fitd,select=-c(1:2) )

	# If SL.randomForest is included in the library,
	# select optimal node size (tree depth) using cross-validated risk
	# and then update the ensemble library to include the optimal node size
	if (length(grep("SL.randomForest",SL.library))>0) {
	  if(is.null(RFnodesize)) RFnodesize <- seq(15,40,by=5)
	  cvRF <- ab_cvRF(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,RFnodesize=RFnodesize)
	  SL.library <- cvRF$SL.library
	}

	# If SL.gam is included in the library,
	# select the optimal degrees of freedom for the smoothing splines using cross-validated risk
	# and then updated the ensemble library to include the optimal df
	if (length(grep("SL.gam",SL.library))>0) {
	  if(is.null(gamdf)) gamdf <- 2:10
	  cvGAM <- ab_cvGAM(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,df=gamdf)
	  SL.library <- cvGAM$SL.library
	}


	# Fit CV.SuperLearner
	cvSL.fit <- CV.SuperLearner(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,V=V,family=family,control=list(saveFitLibrary=TRUE))
	return(cvSL.fit)
}




