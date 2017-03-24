## ----algorithm libraries,include=FALSE-----------------------------------
# load other algorithm packages to clean up the output for agecurveAb(), below
require(tmle)
require(SuperLearner)
require(randomForest)
require(arm)
require(gam)
require(polspline)
require(glmnet)
require(scales)

## ----install,eval=FALSE--------------------------------------------------
#  library(devtools)
#  devtools::install_github("ben-arnold/tmleAb")

## ----read data-----------------------------------------------------------
# load the package into R
library(tmleAb)
# load the Garki serology data (included in the package)
data(garki_sero)
# and subset it to measurement rounds 3-5, ages <20, non-missing IFA Pf titres
d <- garki_sero
d <- subset(d,serosvy>=3 & serosvy<=5 & ageyrs<=20 & !is.na(ifatpftitre))
# create a log-transformed version of the IFA titres (add 1 for 0 values)
d$logpftitre <- log10(d$ifatpftitre+1)
# further subset the data by intervention group for convenience
dc <- d[d$tr=="Control" ,]
di <- d[d$tr=="Intervention",]

## ----fit EYax int,cache=TRUE---------------------------------------------
# list of models and algorithms to include:
mylib <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

# Intervention villages (X=1), fitted curve
set.seed(625234)
fit_i <- agecurveAb(Y=di$logpftitre,Age=di$ageyrs,id=di$id,SL.library=mylib)

## ----fit EYax int2,cache=TRUE--------------------------------------------
mylib2 <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess","SL.randomForest")
set.seed(625234)
fit_i2 <- agecurveAb(Y=di$logpftitre,Age=di$ageyrs,id=di$id,SL.library=mylib2)

## ----fit EYax control,cache=TRUE-----------------------------------------
# Control villages (X=0), fitted curve
fit_c <- agecurveAb(Y=dc$logpftitre,Age=dc$ageyrs,id=dc$id,SL.library=mylib)

## ----plot control and intervention---------------------------------------
# grab a color blind friendly palette
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- cbPalette[c(7,6)]
# plot the curves
plot(fit_c$Age,fit_c$pY,type="l",col=cols[1],
    main="Age-dependent antibody curves\nin control and intervention villages",
    xlim=c(0,20),xlab="Age, years",
    ylim=c(0,4),ylab="Log10 P. Falciparum IFA antibody titre",
    bty="n",las=1)
lines(fit_i$Age,fit_i$pY,col=cols[2],lty=2)
legend("topleft",legend=c("Control","Intervention"),col=cols,lty=c(1,2),bty="n",horiz=T)

## ----plot intervention---------------------------------------------------
# plot the curves
cols <- cbPalette[c(1,7)]
plot(fit_i2$Age,fit_i2$pY,type="l",col=cols[1],
    main="Age-dependent antibody curves\nin intervention villages with different ensembles",
    xlim=c(0,20),xlab="Age, years",
    ylim=c(0,4),ylab="Log10 P. Falciparum IFA antibody titre",
    bty="n",las=1)
lines(fit_i$Age,fit_i$pY,col=cols[2],lty=1)
legend("topleft",legend=c("Incl. random forest","Excl. random forest"),col=cols,lty=c(1,1),bty="n",horiz=T)

# add the underlying raw data points as well
library(scales)
points(fit_i$Age,fit_i$Y,pch=16,cex=0.2,col=alpha("black",alpha=0.5))

## ----fit EYax int adj,cache=TRUE-----------------------------------------
# fit the age-antibody curve in intervention villages, adjusting for sex, season, and village
fit_iadj <- agecurveAb(Y=di$logpftitre,Age=di$ageyrs,id=di$id,
                       W=subset(di,select=c("sex","wetseason","vname")),
                       SL.library=mylib)

## ----plot intervention adj-----------------------------------------------
# compare the results with the unadjusted curve
cols <- cbPalette[c(1,7)]
plot(fit_i$Age,fit_i$pY,type="l",col=cols[1],
    main="Unadjusted and marginally adjusted\nage-dependent antibody curves in intervention villages",
    xlim=c(0,20),xlab="Age, years",
    ylim=c(0,4),ylab="Log10 P. Falciparum IFA antibody titre",
    bty="n",las=1)
lines(fit_iadj$Age,fit_iadj$pY,col=cols[2],lty=1)
legend("topleft",legend=c("Unadjusted","Adjusted"),col=cols,lty=c(1,1),bty="n",horiz=T)


## ----fit EYx,cache=TRUE--------------------------------------------------
# list of models and algorithms to include:
mylib <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

# Intervention villages (X=1)
set.seed(625234)
Wi <- data.frame(Age=di$ageyrs)
mu_i <- tmleAb(Y=di$logpftitre,X=NULL,W=Wi,id=di$id,SL.library=mylib)

# Control villages (X=0)
Wc <- data.frame(Age=dc$ageyrs)
mu_c <- tmleAb(Y=dc$logpftitre,X=NULL,W=Wc,id=dc$id,SL.library=mylib)

## ----test for diff in EYx,cache=TRUE-------------------------------------
# make a binary indicator variable for intervention vs. control
d$tr01 <- ifelse(d$tr=="Intervention",1,0)
W <- data.frame(Age=d$ageyrs)
# estimate the difference between groups
psi_diff <- tmleAb(Y=d$logpftitre,X=d$tr01,W=W,id=d$id,SL.library=mylib)

