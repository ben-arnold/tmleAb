
#' slab_tmle
#'
#' Targeted maximum liklihood estimation (TMLE) for mean antibody levels or differences in antibody levels
#'
#' @param Y antibody measurement
#' @param Age individual's age
#' @param X comparison group  (must be binary, 0/1, for tmle() ).
#' @param W matrix of covariates (optional)
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
slab_tmle <- function(Y,Age,X=NULL,W=NULL,id=NULL,family=gaussian(),SL.library=c("SL.mean","SL.glm","SL.loess","SL.gam","SL.randomForest"),diff=FALSE,RFnodesize=NULL,gamdf=NULL) {

	# if X is null, create a row of 1s
	if (is.null(X)) {
		X <- rep(1,length(Y))
	}

  # if id is null, fill it in
  if (is.null(id)) {
    id <- 1:length(Y)
  }

  # restrict dataset to non-missing observations
  if (is.null(W)) {
    nullW <-TRUE
    fulld <- data.frame(id,Y,X,Age)
  } else{
    nullW <-FALSE
    # convert W into a design matrix (SuperLearner does not acommodate factor variables)
    Wdesign <- design_matrix(W)
    fulld <- data.frame(id,Y,X,Age,Wdesign)
  }
  fitd <- fulld[complete.cases(fulld),]
  fitW <- fitd[,-c(1:3)]

	# If SL.randomForest is included in the library,
	# select optimal node size (tree depth) using cross-validated risk
	# and then update the ensemble library to include the optimal node size
	if (length(grep("SL.randomForest",SL.library))>0) {
	  if(is.null(RFnodesize)) RFnodesize <- seq(15,40,by=5)
	  cvRF <- slab_cvRF(Y=fitd$Y,X=X,id=fitd$id,family=family,SL.library=SL.library,RFnodesize=RFnodesize)
	  SL.library <- cvRF$SL.library
	}

	# If SL.gam is included in the library,
	# select the optimal degrees of freedom for the smoothing splines using cross-validated risk
	# and then updated the ensemble library to include the optimal df
	if (length(grep("SL.gam",SL.library))>0) {
	  if(is.null(gamdf)) gamdf <- 2:10
	  cvGAM <- slab_cvGAM(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,df=gamdf)
	  SL.library <- cvGAM$SL.library
	}

	# estimate either the difference (A=fitd$X) or the adjusted mean (A=NULL)
	if (diff==TRUE) {
		tmle.fit <- tmle(Y=fitd$Y,
			A=fitd$X,
			W=fitW,
			id=fitd$id,
			Q.SL.library=SL.library
		)
		print(tmle.fit)
		psi  <- tmle.fit$estimates$ATE$psi
		se   <- sqrt(tmle.fit$estimates$ATE$var.psi)
		ci   <- tmle.fit$estimates$ATE$CI
		p    <- tmle.fit$estimates$ATE$pvalue
	} else {
		tmle.fit <- tmle(Y=fitd$Y,
			A=NULL,
			W=data.frame(fitd$Age,fitd$W),
			id=fitd$id,
			Q.SL.library=SLlib
		)
		print(tmle.fit)
		psi  <- tmle.fit$estimates$EY1$psi
		se   <- sqrt(tmle.fit$estimates$EY1$var.psi)
		ci   <- tmle.fit$estimates$EY1$CI
		p    <- tmle.fit$estimates$EY1$pvalue
	}
	# return estimate, SE, 95% CI, and P-value
	list(psi=psi,se=se,lb=ci[1],ub=ci[2],p=p)
}






