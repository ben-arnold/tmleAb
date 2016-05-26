
#----------------------------------------------------------------
#' Convert a vector, matrix or data frame with factors into a design matrix
#'
#' An internal, convenience function that automatically transforms a vector, matrix, or data.frame with factors into a design matrix with indicator variables and an ommitted category.
#'
#' @param W A vector, matrix, or data.frame that includes numeric or factor variables.
#'
#' @return A design matrix version of \code{W} where factor variables have been converted into columns of indicator variables with the first level excluded.
#'
#' @examples
#'
design_matrix <- function(W) {
  # W : data frame of covariates that might include factors
  if(class(W)!="matrix" & class(W)!="data.frame"){
    cat("\n-----------------------------------------\nThe design matrix you supplied is not a matrix or a data.frame\nAssuming that it is a single variable\n-----------------------------------------\n")
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

