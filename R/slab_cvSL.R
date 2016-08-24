
#' SLAb wrapper function to estimate CV risk for the super learner ensemble predicting antibody measurements using age and potentially other covariates
#'
#' @param Y log antibody measurement
#' @param Age
#' @param W
#' @param id
#' @param family
#' @param V
#' @param SL.library
#'
#' @return
#' @export
#'
#' @examples
slab_cvSL <-function(Y,Age,W=NULL,id=1:length(Y),family=gaussian(),V=10,SL.library= c("SL.mean","SL.glm","SL.bayesglm","SL.loess","SL.gam","SL.randomForest")) {
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
	  Wdesign <- design.matrix(W)
	  fitd <- data.frame(id,Y,Age,Wdesign)
	}

	# restrict dataset to non-missing observations
	fitd <- fitd[complete.cases(fitd),]
	X <- subset(fitd,select=-c(1:2) )

  # If randomForest is included in the library,
  # select optimal node size (tree depth) using cross-validated risk
  # and then update the ensemble library to include the optimal node size
  if (grep("SL.randomForest",SL.library)>0) {
    cvRF <- SLAb.cvRF(Y=fitd$Y,X=X,id=fitd$id,SLlib=SL.library)
    SL.library <- cvRF$SLlib
  }

	# Fit CV.SuperLearner
	cvSL.fit <- CV.SuperLearner(Y=fitd$Y,X=X,id=fitd$id,SL.library=SL.library,V=V,family=family,control=list(saveFitLibrary=TRUE))
	return(cvSL.fit)
}



# --------------------------------------
# automatic transform of a covariate
# data.frame with factors into one
# with indicator variables (and an
# ommitted category)
# --------------------------------------
design.matrix <- function(W) {
  # W : data frame of covariates that might include factors
  if(class(W)!="matrix" & class(W)!="data.frame"){
    #cat("\n-----------------------------------------\nThe design matrix you supplied is not a matrix or a data.frame\nAssuming that it is a single variable\n-----------------------------------------\n")
    W <- data.frame(W)
  }
  ncolW <- ncol(W)
  flist <- numeric()
  for(i in 1:ncolW) {
    if(class(W[,i])!="factor"){
      next
    } else {
      flist <- c(flist,i)
      # strip out extra levels
      W[,i] <- factor(W[,i])
      # create a design matrix, remove the first level
      mm <- model.matrix(~-1+W[,i])
      mW <- mm[,-c(1)]
      # format the names of the indicator variables
      # and add them to the design matrix
      levs <- gsub(" ","",levels(W[,i]) )[-c(1)]
      if(length(levs)<2) mW <- matrix(mW,ncol=1)
      colnames(mW) <- paste(names(W)[i],levs,sep="")
      W <- data.frame(W,mW)
    }
  }
  # now drop the factors that have been replaced by indicators
  W <- W[,-c(flist)]
  return(W)
}


