

#----------------------------------------------------------------
#' Fit an age-dependent antibody curve with ensemble machine learning
#'
#' \code{agecurveAb} uses ensemble machine learning with \code{\link[SuperLearner]{SuperLearner}} to fit population mean antibody curves by age. If covariates \code{W} are included, the predicted curve is marginally adjusted for \code{W}.
#'
#' @param Y Antibody measurement. Must be a numeric vector.
#' @param Age Age of the individual at the time of measurement. Must be a numeric vector.
#' @param W An optional vector, matrix, or data.frame of covariates for each individual used to marginally adjust the curve
#' @param id An optional cluster or repeated measures id variable. For cross-validation splits, \code{id} forces observations in the same cluster or for the same individual to be in the same validation fold.
#' @param family Model family (gaussian for continuous outcomes, binomial for binary outcomes)
#' @param SL.library Library of algorithms to include in the ensemble (see the \code{\link[SuperLearner]{SuperLearner}} package for details).
#' @param RFnodesize Optional argument to specify a range of minimum node sizes for the random Forest algorithm. If \code{SL.library} includes \code{SL.randomForest}, then the default is to search over node sizes of 15,20,...40. Specifying this option will override the default.
#' @param gamdf Optional argument to specify a range of degrees of freedom for natural smoothing splines in a generalized additive model. If \code{SL.library} includes \code{SL.gam}, then the default is to search over a range of df=2-10. Specifying this option will override the default.
#'
#' @return \code{agecurveAb} returns a data.frame, which includes the dataset (feature matrix) used for estimation, along with fitted results (\code{pY}). Note that the estimation dataset excludes any observations with missing values in \code{Y}, \code{Age}, \code{W} (if not NULL), or \code{id} (if specified). Factors in \code{W} are converted to design-matrix-style indicator variables. Observations are sorted by \code{Age} for more convenient plotting. If covariates are included, then \code{pY} is the mean predicted antibody level at \code{Age=a}, averaged over the covariates \code{W}.
#'
#'
#'
#' @details The \code{agecurveAb} function is a wrapper for \code{\link[SuperLearner]{SuperLearner}} that provides a convenient interface for this specific estimation problem. If the \code{SL.library} argument includes just one model or algorithm, then there is no 'ensemble' but the function provides a standard interface for using single algorithms (e.g., \code{\link[stats]{SL.loess}})  Note that if \code{SL.randomForest} is included in the library, \code{agecurveAb} will select the minimum node size (between 15 and 40) with cross-validation to avoid over-fitting. If you wish to control the randomForest node size options using a range other than 15-40, you can do so by passing an argument \code{RFnodesize} through this function. Similarly, if \code{SL.gam} is included in the library, \code{agecurveAb} will select the optimal degrees of freedom for natural splines (between 2 and 10) with cross-validation to get the correct amount of smoothing.  If you wish to control the GAM df search, you can do so by passing an argument \code{gamdf} through this function.
#'
#' @seealso \code{\link{tmleAb}}
#' @seealso \code{\link[SuperLearner]{SuperLearner}}
#'
#' @references van der Laan MJ, Polley EC, Hubbard AE. Super Learner. Stat Appl Genet Mol Biol. 2007;6: 1544–6115. \link{http://www.ncbi.nlm.nih.gov/pubmed/17910531}
#'
#'
#' @export
#'
#' @examples
#' # load the Garki project serology data
#' data("garki_sero")
#' garki_sero$village <- factor(garki_sero$village)
#' garki_sero$sex <- factor(garki_sero$sex)
#'
#' # control village measurements in round 5
#' dc <- subset(garki_sero,serosvy==5 & tr=="Control")
#'
#' # intervention village measurements in round 5
#' di <- subset(garki_sero,serosvy==5 & tr=="Intervention")
#'
#' # fit an age-antibody curve in control and intervention villages
#' # adjusted for sex and village
#' # set a seed for perfectly reproducible splits in the V-fold cross validation
#' set.seed(12345)
#' ccurve <-agecurveAb(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,W=dc[,c("sex","village")],id=dc$id)
#' set.seed(12345)
#' icurve <-agecurveAb(Y=log10(di$ifatpftitre+1),Age=di$ageyrs,W=dc[,c("sex","village")],id=di$id)
#'
#' # plot the curves
#' plot(ccurve$Age,ccurve$pY,type="l",ylim=c(0,4),bty="l",las=1)
#' lines(icurve$Age,icurve$pY,col="blue")
#'
agecurveAb <-function(Y,Age,W=NULL,id=NULL,family=gaussian(),SL.library= c("SL.mean","SL.glm","SL.gam","SL.loess"), RFnodesize=NULL,gamdf=NULL) {

  # ensure SuperLeaner package is loaded
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("The SuperLearner package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.null(id)) id <- 1:length(Y)

  # restrict dataset to non-missing observations
  if (is.null(W)) {
    nullW <-TRUE
    fulld <- data.frame(id,Y,Age)
  } else{
    nullW <-FALSE
    # convert W into a design matrix (SuperLearner does not acommodate factor variables)
    Wdesign <- design_matrix(W)
    fulld <- data.frame(id,Y,Age,Wdesign)
  }
  fitd <- fulld[complete.cases(fulld),]

  # throw a warning if observations have missing data
  n.orig <- dim(fulld)[1]
  n.fit  <- dim(fitd)[1]
  if(n.orig>n.fit) warning(paste("\n\n",n.orig-n.fit,'observations were dropped due to missing values\n in the outcome, age, or adjustement covariates. \n The original dataset contained',n.orig,'observations,\n but agecurveAb is fitting the curve using',n.fit,'observations.'))

  #  matrix of features used for SuperLearner prediction (drop id and outcome, Y)
  X <- subset(fitd,select=-c(1:2) )

  # If SL.randomForest is included in the library,
  # select optimal node size (tree depth) using cross-validated risk
  # and then update the ensemble library to include the optimal node size
  if (length(grep("SL.randomForest",SL.library))>0) {
    cvRF <- ab_cvRF(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,RFnodesize=RFnodesize)
    SL.library <- cvRF$SL.library
  }

  # If SL.gam is included in the library,
  # select the optimal degrees of freedom for the smoothing splines using cross-validated risk
  # and then updated the ensemble library to include the optimal df
  if (length(grep("SL.gam",SL.library))>0) {
    cvGAM <- ab_cvGAM(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,df=gamdf)
    SL.library <- cvGAM$SL.library
  }

  # Fit SuperLearner
  SLfit <- SuperLearner::SuperLearner(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,family=family,method="method.NNLS")
  p_res <- cbind(SLfit$cvRisk,SLfit$coef)
  colnames(p_res) <- c("CV-Risk","Coef")
  cat("\nSummary of SuperLearner cross validated risk and \nweights for algorithms included in the library:\n\n")
  print(p_res)

  # obtain marginally averaged, predicted values of Y at A=a:  E_W[E(Y|A=a, X=x, W)]
  # with X=x implied by the subset of data used to fit the function
  # this parsing below is just to speed up computation in the case of
  # no adjustment covariates W=NULL
  # if W is not null, then need to do marginal averaging at each age
  # and merge that back to the analysis data frame, which is slower
  if(nullW==TRUE) {
    res <- fitd
    res$pY <- predict(SLfit)$pred
  } else{
    As <- unique(X$Age)
    pY <- rep(NA,length(As))
    for(i in 1:length(As)) {
      X$Age <- As[i]
      pYs <-predict(SLfit,newdata=X)$pred
      pY[i] <- mean(pYs)
    }
    res <- merge(fitd,data.frame(Age=As,pY=pY),by="Age",all.x=T,all.y=T)
  }

  # return the data frame with predicted values, sorted by Age for convenience
  res <- res[order(res$Age),]
  return(res)
}


