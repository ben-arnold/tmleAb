
#' Targeted maximum liklihood estimation for antibody measurements
#'
#' Targeted maximum liklihood estimation (TMLE) for mean antibody levels or differences in antibody levels
#'
#' @param Y antibody measurement
#' @param X comparison group  (must be binary, 0/1). If \code{X=NULL}, then the function returns the mean (rather than the difference between levels of \code{X}).
#' @param W matrix of covariates -- should probably at minimum include the individual's age (if available).
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Model family (gaussian for continuous outcomes, binomial for binary outcomes)
#' @param SL.library Library of algorithms to include in the ensemble (see the \code{\link[SuperLearner]{SuperLearner}} package for details).
#' @param RFnodesize Optional argument to specify a range of minimum node sizes for the random Forest algorithm. If \code{SL.library} includes \code{SL.randomForest}, then the default is to search over node sizes of 15,20,...40. Specifying this option will override the default.
#' @param gamdf Optional argument to specify a range of degrees of freedom for natural smoothing splines in a generalized additive model. If \code{SL.library} includes \code{SL.gam}, then the default is to search over a range of df=2-10. Specifying this option will override the default.
#'
#' @details TBD
#'
#' @return \code{psi} Mean (if \code{X=NULL}) or difference
#' @return \code{se} Standard error of \code{psi}, estimated from the influence curve
#' @return \code{lb} Lower bound of the 95 percent confidence interval of \code{psi}
#' @return \code{ub} Upper bound of the 95 percent confidence interval of \code{psi}
#' @return \code{p} P-value for a test that \code{psi=0}
#' @return \code{tmle_fit} The original \code{tmle()} fit (see the \code{\link[tmle]{tmle}} package for details).
#' @export
#'
#' @examples TBD
tmleAb <- function(Y,X=NULL,W=NULL,id=NULL,family="gaussian",SL.library=c("SL.mean","SL.glm","SL.gam","SL.loess"),RFnodesize=NULL,gamdf=NULL) {

  # ensure SuperLeaner and tmle packages are loaded
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("You need the SuperLearner package for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("tmle", quietly = TRUE)) {
    stop("You need the tmle package for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.null(id)) {
    id <- 1:length(Y)
  }

  # restrict dataset to non-missing observations
  if (is.null(W)) {
    nullW <-TRUE
    fulld <- data.frame(id,Y)
  } else{
    nullW <-FALSE
    # convert W into a design matrix (SuperLearner does not acommodate factor variables)
    Wdesign <- design_matrix(W)
    fulld <- data.frame(id,Y,Wdesign)
  }

  if (is.null(X)) {
    nullX <- TRUE
  } else{
    nullX <- FALSE
    if(min(X,na.rm=T)!=0 & max(X,na.rm=T)!=1){
      error("X must be a binary variable 0/1")
    }
    fulld$X <- X
  }

  fitd <- fulld[complete.cases(fulld),]
  fitW <- subset(fitd,select=c(-id,-Y))

	# If SL.randomForest is included in the library,
	# select optimal node size (tree depth) using cross-validated risk
	# and then update the ensemble library to include the optimal node size
	if (length(grep("SL.randomForest",SL.library))>0) {
	  cvRF <- ab_cvRF(Y=fitd$Y,X=fitW,id=fitd$id,family=family,SL.library=SL.library,RFnodesize=RFnodesize)
	  SL.library <- cvRF$SL.library
	}

	# If SL.gam is included in the library,
	# select the optimal degrees of freedom for the smoothing splines using cross-validated risk
	# and then updated the ensemble library to include the optimal df
	if (length(grep("SL.gam",SL.library))>0) {
	  cvGAM <- ab_cvGAM(Y=fitd$Y,X=fitW,id=fitd$id,SL.library=SL.library,df=gamdf)
	  SL.library <- cvGAM$SL.library
	}

	# estimate either the difference (A=fitd$X) or the adjusted mean (A=NULL)
	if (nullX==FALSE) {
		tmle_fit <- tmle::tmle(Y=fitd$Y,
			A=fitd$X,
			W=fitW,
			id=fitd$id,
			Q.SL.library=SL.library,
			family=family,
			fluctuation = "logistic"
		)
		tmle::print.tmle(tmle_fit)
		psi  <- tmle_fit$estimates$ATE$psi
		se   <- sqrt(tmle_fit$estimates$ATE$var.psi)
		ci   <- tmle_fit$estimates$ATE$CI
		p    <- tmle_fit$estimates$ATE$pvalue
	} else {
		tmle_fit <- tmle::tmle(Y=fitd$Y,
			A=NULL,
			W=fitW,
			id=fitd$id,
			Q.SL.library=SL.library,
			family=family,
			fluctuation = "logistic"
		)
		tmle::print.tmle(tmle_fit)
		psi  <- tmle_fit$estimates$EY1$psi
		se   <- sqrt(tmle_fit$estimates$EY1$var.psi)
		ci   <- tmle_fit$estimates$EY1$CI
		p    <- tmle_fit$estimates$EY1$pvalue
	}
	# return estimate, SE, 95% CI, and P-value, along with the tmle() object
	list(psi=psi,se=se,lb=ci[1],ub=ci[2],p=p,tmle_fit=tmle_fit)
}






