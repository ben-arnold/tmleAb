
#--------------------------------
# generalized plotting routine
# for cross-validated super learner
# predictions of antibody data
#--------------------------------

#' Plotting routine for cross-validated risk of super learner predictions of antibody data (along with its constituent algorithms)
#'
#' @param x Object of class CV.SuperLearner
#' @param col (optional) Colors for plotting
#' @param title (optional) Title for the plot
#' @param xlim (optional) X-axis limits
#'
#' @return
#' @export
#'
#' @examples
slab_plot_cvSL <- function(x,col='blue',title=NULL,xlim=NULL) {

# load ggplot2
require(ggplot2)

# summarize the object
sumx <- summary(x)

# remove the discrete super learner
sumx$Table <- subset(sumx$Table,Algorithm!="Discrete SL")
sumx$Table$Algorithm <- factor(sumx$Table$Algorithm)

# for some of the commonly used algorithms, clean up the names
sumx$Table$Algorithm <- as.character(sumx$Table$Algorithm)
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.mean_All"] <- "Simple Mean"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.glm_All"] <- "GLM"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.bayesglm_All"] <- "Bayes GLM"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.gam_All"] <- "GAM"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.loess_All"] <- "LOESS"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.polymars_All"] <- "MARS"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.glmnet_All"] <- "Lasso (glmnet)"
sumx$Table$Algorithm[sumx$Table$Algorithm=="SL.Yman2016_All"] <- "Yman2016"
sumx$Table$Algorithm[grep("SL.randomForest",sumx$Table$Algorithm)] <- "Random Forest"
sumx$Table$Algorithm <- as.factor(sumx$Table$Algorithm)

# assign colors (first SL, then others in their original call order)
nalg <- length(levels(sumx$Table$Algorithm))
if(length(col)>1) {
  if(length(col)!=nalg) {
    warning("col must be either length 1 or the same as the number of learners +1 (for SL)\njust using the first color")
    col <- as.character(rep(col[1],nalg))
  }
} else{
  col <- as.character(rep(col,nalg))
}
names(col) <- sumx$Table$Algorithm


# sort algorithm by CV risk
# always put Super Learner at the top
sumx$Table$Algorithm <- reorder(sumx$Table$Algorithm, -sumx$Table$Ave)
sumx$Table$Algorithm <- factor(sumx$Table$Algorithm, levels=c(levels(sumx$Table$Algorithm)[-grep("Super Learner",levels(sumx$Table$Algorithm))],"Super Learner") )

# plot
Mean <- sumx$Table$Ave
se <- sumx$Table$se
Lower <- Mean - qnorm(0.975) * se
Upper <- Mean + qnorm(0.975) * se
assign("d", data.frame(Y = Mean, X = sumx$Table$Algorithm,
                       Lower = Lower, Upper = Upper))

  p <- ggplot(d, aes_string(x = "X", y = "Y", ymin = "Lower",  ymax = "Upper", color = "X")) +
    geom_linerange() + geom_point(size=4) + coord_flip() +
    ggtitle(title) +
    ylab("V-fold CV Risk Estimate") + xlab("Method") +
    scale_colour_manual(name = "X",values = col)  +
    theme_bw() + theme(legend.position="none")

  return(p)

}
