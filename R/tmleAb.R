
#' Targeted maximum liklihood estimation for antibody measurements
#'
#' Targeted maximum liklihood estimation (TMLE) for mean antibody measurements or difference in antibody measurements between groups
#'
#' @param Y antibody measurement
#' @param X comparison group  (must be binary, 0/1). If \code{X=NULL}, then the function returns the mean (rather than the difference between levels of \code{X}).
#' @param W matrix of covariates -- should probably at minimum include the individual's age (if available).
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Outcome family, choose \code{gaussian} for continuous outcomes and \code{binomial} for binary outcomes (default \code{family="gaussian"})
#' @param SL.library Library of models/algorithms to include in the ensemble for the outcome (see the \code{\link[SuperLearner]{SuperLearner}} package for details).
#' @param g.SL.library Optional library of models/algorithms to model group assignment. Default is to use main terms logistic regression (SL.glm).
#' @param V Number of cross-validation folds for Super Learning to estimate outcome (Q) and treatment (g) models (default is \code{V=5}).
#' @param RFnodesize Optional argument to specify a range of minimum node sizes for the random Forest algorithm. If \code{SL.library} includes \code{SL.randomForest}, then the default is to search over node sizes of 15,20,...40. Specifying this option will override the default.
#' @param gamdf Optional argument to specify a range of degrees of freedom for natural smoothing splines in a generalized additive model. If \code{SL.library} includes \code{SL.gam}, then the default is to search over a range of df=2-10. Specifying this option will override the default.
#'
#' @details
#' The \code{tmleAb} function estimates adjusted means or differences in means in antibody measurements using targeted maximum likelihood estimation (TMLE). 
#' 
#' The function assumes a continuous outcome as the default (\code{family="gaussian"}). If you pass a binary outcome to the function with the \code{family="gaussian"} argument it will still estimate seroprevalence, but it will not necessarily bound predictions between 0 and 1. If you specify \code{family="binomial"} then the predictions will be bound between 0 and 1. Note that some estimation routines do not support binary outcomes (e.g., \code{SL.loess}), and you will see an error if you specify a binomial family with them in the library.  
#' 
#' Also note that if you specify \code{family="binomial"} for a binary outcome with a comparison group (\code{X=}), then other contrasts are available to you besides the difference in means (the default stored in \code{psi}). Specifically the relative risk (RR) and odds ratio (OR) will be saved in the \code{tmle_fit$estimates} list.
#'
#' @return \code{psi} Mean (if \code{X=NULL}) or difference
#' @return \code{se} Standard error of \code{psi}, estimated from the influence curve
#' @return \code{lb} Lower bound of the 95 percent confidence interval of \code{psi}
#' @return \code{ub} Upper bound of the 95 percent confidence interval of \code{psi}
#' @return \code{p} P-value for a test that \code{psi=0}
#' @return \code{tmle_fit} The original \code{tmle()} fit (see the \code{\link[tmle]{tmle}} package for details).
#'
#' @seealso
#' \code{\link[tmleAb]{agecurveAb}}, \code{\link[SuperLearner]{SuperLearner}},  \code{\link[tmle]{tmle}}
#'
#' @references Gruber S, van der Laan M. tmle: An R Package for Targeted Maximum Likelihood Estimation. J Stat Softw. 2012;51: 1–35.
#' @references van der Laan MJ, Polley EC, Hubbard AE. Super Learner. Stat Appl Genet Mol Biol. 2007;6: 1544–6115.
#'
#' @examples
#' \dontrun{
#' # load the Garki project serology data, subset to round 5
#' data("garki_sero")
#' garki_sero$village <- factor(garki_sero$village)
#' garki_sero$sex <- factor(garki_sero$sex)
#' garki_sero$tr01 <- ifelse(garki_sero$tr=="Intervention",1,0)
#' d <- subset(garki_sero,serosvy==5)
#'
#' # control and intervention village measurements
#' dc <- subset(d, tr=="Control")
#' di <- subset(d,tr=="Intervention")
#'
#' # estimate means in control and intervention villages
#' # set a seed for perfectly reproducible splits
#' # in the V-fold cross validation
#' set.seed(12345)
#' cmean <-tmleAb(Y=log10(dc$ifatpftitre+1),
#'                W=dc[,c("ageyrs","sex","village")],
#'                id=dc$id)
#' set.seed(12345)
#' imean <-tmleAb(Y=log10(di$ifatpftitre+1),
#'                W=di[,c("ageyrs","sex","village")],
#'                id=di$id)
#'
#' # estimate targeted difference in means
#' # bewteen control and intervention villages
#' # adjusted for age, sex and village
#' set.seed(12345)
#' psi_diff <-tmleAb(Y=log10(d$ifatpftitre+1),
#'                   X=d$tr01,
#'                   W=d[,c("ageyrs","sex","village")],
#'                   id=d$id)
#'}
#'
#' @export
#'
tmleAb <- function(Y,X=NULL,W=NULL,id=NULL,family="gaussian",SL.library=c("SL.mean","SL.glm","SL.gam","SL.loess"),g.SL.library=c("SL.glm"),V=5,RFnodesize=NULL,gamdf=NULL) {

  # ensure SuperLeaner and tmle packages are loaded
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("You need the SuperLearner package for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("tmle", quietly = TRUE)) {
    stop("You need the tmle package for this function to work. Please install it.",
         call. = FALSE)
  }
  require(SuperLearner)
  require(tmle)

  if (is.null(id)) {
    id <- 1:length(Y)
  }

  # restrict dataset to non-missing observations
  if (is.null(W)) {
    nullW <-TRUE
    fulld <- data.frame(id,Y)
  } else{
    nullW <-FALSE
    # convert W into a design matrix (SuperLearner currently does not acommodate factor variables)
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
	  cvRF <- ab_cvRF(Y=fitd$Y,X=fitW,id=fitd$id,family=family,SL.library=SL.library,cvControl=list(V = V),RFnodesize=RFnodesize)
	  SL.library <- cvRF$SL.library
	}

	# If SL.gam is included in the library,
	# select the optimal degrees of freedom for the smoothing splines using cross-validated risk
	# and then updated the ensemble library to include the optimal df
	if (length(grep("SL.gam",SL.library))>0) {
	  cvGAM <- ab_cvGAM(Y=fitd$Y,X=fitW,id=fitd$id,family=family,SL.library=SL.library,cvControl=list(V = V),df=gamdf)
	  SL.library <- cvGAM$SL.library
	}

	# estimate either the difference (A=fitd$X) or the mean (A=NULL), possibly adjusted for W
	if (nullX==FALSE & nullW==FALSE) {
		tmle_fit <- tmle::tmle(Y=fitd$Y,
			A=fitd$X,
			W=subset(fitW,select=c(-X)),
			id=fitd$id,
			Q.SL.library=SL.library,
			g.SL.library=g.SL.library,
			V=V,
			family=family,
			fluctuation = "logistic"
		)
		tmle::print.tmle(tmle_fit)
		psi  <- tmle_fit$estimates$ATE$psi
		se   <- sqrt(tmle_fit$estimates$ATE$var.psi)
		ci   <- tmle_fit$estimates$ATE$CI
		p    <- tmle_fit$estimates$ATE$pvalue
	}
  if (nullX==FALSE & nullW==TRUE) {
    # if W is null, then create a random normal variable
    # called emptyW -- it will just enable tmle to run (band-aid solution)
    emptyW <- data.frame(w1=rnorm(n=nrow(fitd)))
    tmle_fit <- tmle::tmle(Y=fitd$Y,
                           A=fitd$X,
                           W=emptyW,
                           id=fitd$id,
                           Q.SL.library=SL.library,
                           g.SL.library=g.SL.library,
                           V=V,
                           family=family,
                           fluctuation = "logistic"
    )
    tmle::print.tmle(tmle_fit)
    psi  <- tmle_fit$estimates$ATE$psi
    se   <- sqrt(tmle_fit$estimates$ATE$var.psi)
    ci   <- tmle_fit$estimates$ATE$CI
    p    <- tmle_fit$estimates$ATE$pvalue
  }
  if (nullX==TRUE & nullW==FALSE)  {
		tmle_fit <- tmle::tmle(Y=fitd$Y,
			A=NULL,
			W=fitW,
			id=fitd$id,
			Q.SL.library=SL.library,
			V=V,
			family=family,
			fluctuation = "logistic"
		)
		tmle::print.tmle(tmle_fit)
		psi  <- tmle_fit$estimates$EY1$psi
		se   <- sqrt(tmle_fit$estimates$EY1$var.psi)
		ci   <- tmle_fit$estimates$EY1$CI
		p    <- tmle_fit$estimates$EY1$pvalue
  }
  if (nullX==TRUE & nullW==TRUE)  {
    # if W is null, then create a random normal variable
    # called emptyW -- it will just enable tmle to run (band-aid solution)
    emptyW <- data.frame(w1=rnorm(n=nrow(fitd)))
    tmle_fit <- tmle::tmle(Y=fitd$Y,
                           A=NULL,
                           W=emptyW,
                           id=fitd$id,
                           Q.SL.library=SL.library,
                           V=V,
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






