

#' Plots of cross-validated R-squared
#'
#'#' Plotting routine for cross-validated R-squared of super learner predictions of antibody data (along with its constituent algorithms)
#'
#' @param cvsl An object of class \code{CV.SuperLearner}
#' @param X The data.frame of features used in the CV.SuperLearner prediction
#' @param col (optional) Colors for plotting
#' @param ylab Name of the loss function (defaults uses generic "Risk") but it can be modified to be more specific, such as MSE, AUC, etc...
#' @param title  Title for the plot (optional)
#' @param ylim  set Y-axis limits, defaults to c(0,1)
#'
#' @return
#' @export
#'
#' @examples
ab_plot_cvR2 <- function(cvsl,X,col='blue',ylab="V-fold Cross Validated R-squared", title=NULL,ylim=c(0,1)) {

  # Compute cross-validated R-squared from the CV risk estimates
  lnames <- c("SuperLearner",cvsl$libraryNames)
  cvr2s <- matrix(NA,nrow=length(lnames),ncol=3)
  for(i in 1:length(lnames)) {
    r2obj <- r2weight(list(cvsl),list(X),whichAlgorithm=lnames[i])
    cvr2s[i,] <- c(r2obj$r2weight$cv.wR2,r2obj$r2weight$cv.wR2.ci)
  }
  cvr2 <- data.frame(Algorithm=lnames,r2=cvr2s[,1],lower=cvr2s[,2],upper=cvr2s[,3])


  # for some of the commonly used algorithms, clean up the names
  cvr2$Algorithm <- as.character(cvr2$Algorithm)
  cvr2$Algorithm[cvr2$Algorithm=="SuperLearner"] <- "Super Learner"
  cvr2$Algorithm[cvr2$Algorithm=="SL.mean_All"] <- "Simple Mean"
  cvr2$Algorithm[cvr2$Algorithm=="SL.glm_All"] <- "GLM"
  cvr2$Algorithm[cvr2$Algorithm=="SL.glm.interaction_All"] <- "GLM Interaction"
  cvr2$Algorithm[cvr2$Algorithm=="SL.bayesglm_All"] <- "Bayes GLM"
  cvr2$Algorithm[cvr2$Algorithm=="SL.loess_All"] <- "LOESS"
  cvr2$Algorithm[cvr2$Algorithm=="SL.polymars_All"] <- "MARS"
  cvr2$Algorithm[cvr2$Algorithm=="SL.glmnet_All"] <- "Lasso (glmnet)"
  cvr2$Algorithm[cvr2$Algorithm=="SL.Yman2016_All"] <- "Yman2016"
  cvr2$Algorithm[grep("SL.gam",cvr2$Algorithm)] <- "GAM"
  cvr2$Algorithm[grep("SL.randomForest",cvr2$Algorithm)] <- "Random Forest"
  cvr2$Algorithm <- as.factor(cvr2$Algorithm)

  # assign colors (first SL, then others in their original call order)
  nalg <- length(levels(cvr2$Algorithm))
  if(length(col)>1) {
    if(length(col)!=nalg) {
      warning("col must be either length 1 or the same as the number of learners +1 (for SL)\njust using the first color")
      col <- as.character(rep(col[1],nalg))
    }
  } else{
    col <- as.character(rep(col,nalg))
  }
  names(col) <- cvr2$Algorithm


  # sort algorithm by CV risk
  # always put Super Learner at the top
  cvr2$Algorithm <- reorder(cvr2$Algorithm, cvr2$r2)
  cvr2$Algorithm <- factor(cvr2$Algorithm, levels=c(levels(cvr2$Algorithm)[-grep("Super Learner",levels(cvr2$Algorithm))],"Super Learner") )

  # print R2 estimates
  print(cvr2,digits=4)

  # plot
  p <- ggplot(cvr2, aes_string(x = "Algorithm", y = "r2", ymin = "lower",  ymax = "upper", color = "Algorithm")) +
    geom_linerange() + geom_point(size=4) + coord_flip() +
    ggtitle(title) +
    ylab(ylab) + xlab("Method") +
    expand_limits(y=ylim) +
    scale_colour_manual(name = "Algorithm",values = col) +
    theme_bw() + theme(legend.position="none")

  return(p)

}
