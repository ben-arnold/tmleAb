#'
#' Estimate cross-validated risk for the super learner fit to antibody data
#'
#' A convenience wrapper for \code{\link[SuperLearner]{CV.SuperLearner}} for antibody measurements.
#'
#' @param Y Antibody measurement. Must be a numeric vector.
#' @param X A vector, matrix, or data.frame of covariates for each individual used to predict antibody levels
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param V Number of folds to use in the cross validation (default is 10)
#' @param SL.library Library of algorithms to include in the ensemble (see the \code{\link[SuperLearner]{SuperLearner}} package for details).
#' @param RFnodesize Optional argument to specify a range of minimum node sizes for the random forest algorithm. If \code{SL.library} includes \code{SL.randomForest}, then the default is to search over node sizes of 15,20,...40. Specifying this option will override the default.
#' @param gamdf Optional argument to specify a range of degrees of freedom for natural smoothing splines in a generalized additive model. If \code{SL.library} includes \code{SL.gam}, then the default is to search over degrees of freedom 2,3,...10. Specifying this option will override the default.
#'
#' @details The SuperLearner function builds a estimator, but does not contain an estimate on the performance of the estimator. Various methods exist for estimator performance evaluation. If you are familiar with the super learner algorithm, it should be no surprise we recommend using cross-validation to evaluate the honest performance of the super learner estimator. The function \code{cvSLAb} provides a convenient wrapper for the \code{CV.SuperLearner} routine to compute the V-fold cross-validated risk estimate for the super learner (and all algorithms in \code{SL.library} for comparison). The wrapper adds convenience by restricting the dataset to complete cases, transforming the covariate matrix (\code{W}) into a data.frame, and allowing the user to tune parameters in the Random Forest and GAM libraries if they are included in \code{SL.library}. It assumes a continuous outcome (\code{family=gaussian()}), but can be run on binary outcomes without problems.
#'
#' @return This function returns an object of class \code{CV.SuperLearner} (see the \code{\link[SuperLearner]{SuperLearner}} package for details)
#'
#' @seealso
#' \code{\link{tmleAb}}, \code{\link[SuperLearner]{SuperLearner}}
#'
#' @examples
#' \dontrun{
#' # load the Garki project serology data, subset to round 5 intervention
#' data("garki_sero")
#' garki_sero$village <- factor(garki_sero$village)
#' garki_sero$sex <- factor(garki_sero$sex)
#' garki_sero$tr01 <- ifelse(garki_sero$tr=="Intervention",1,0)
#' d <- subset(garki_sero,serosvy==5 & tr=="Intervention")
#'
#' # fit the cross-validated super learner
#' # with just Age as the predictor
#' set.seed(62522)
#' CVfit <- cvSLAb(Y=log10(d$ifatpftitre+1),X=data.frame(Age=d$ageyrs),id=d$id)
#'
#' # plot cross-validated MSE ("Risk")
#' plot(CVfit)
#'}
#' @export
#'
cvSLAb <-function(Y,X,id=1:length(Y),V=10,SL.library= c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.randomForest"),RFnodesize=NULL,gamdf=NULL) {

  # convert X into a design matrix (SuperLearner does not acommodate factor variables)
  if (is.null(X)) {
    stop("You must include at least one covariate in X, such as age.")
	} else{
	  Xdesign <- design_matrix(X)
	  fitd <- data.frame(id,Y,Xdesign)
	}

	# restrict dataset to non-missing observations
	fitd <- fitd[complete.cases(fitd),]
	fitX <- subset(fitd,select=-c(1:2) )

	# If SL.randomForest is included in the library,
	# select optimal node size (tree depth) using cross-validated risk
	# and then update the ensemble library to include the optimal node size
	if (length(grep("SL.randomForest",SL.library))>0) {
	  cvRF <- ab_cvRF(Y=fitd$Y,X=fitX,id=fitd$id,SL.library=SL.library,RFnodesize=RFnodesize)
	  SL.library <- cvRF$SL.library
	}

	# If SL.gam is included in the library,
	# select the optimal degrees of freedom for the smoothing splines using cross-validated risk
	# and then updated the ensemble library to include the optimal df
	if (length(grep("SL.gam",SL.library))>0) {
	  cvGAM <- ab_cvGAM(Y=fitd$Y,X=fitX,id=fitd$id,SL.library=SL.library,df=gamdf)
	  SL.library <- cvGAM$SL.library
	}


	# Fit CV.SuperLearner
	cvSL.fit <- SuperLearner::CV.SuperLearner(Y=fitd$Y,X=fitX,id=fitd$id,SL.library=SL.library,V=V,family="gaussian",control=list(saveFitLibrary=TRUE))
	return(cvSL.fit)
}




