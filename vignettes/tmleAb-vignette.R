## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)

## ----read data-----------------------------------------------------------
# load the package into R
library(tmleAb)
library(SuperLearner)
# load the Garki serology data and subset it to measurement round 5
d <- garki_sero
d <- subset(d,serosvy==5)
# further subset the data by intervention group for convenience
dc <- d[d$tr=="Control" ,]
di <- d[d$tr=="Intervention",]

## ----fit EYax control,cache=TRUE-----------------------------------------
set.seed(625234)
# list of models and algorithms to include:
mylib <- c("SL.mean","SL.glm","SL.Yman2016","SL.gam","SL.loess")

# Control villages (X=0), fitted curve
fit_c <- agecurveAb(Y=log10(dc$ifatpftitre+1),Age=dc$ageyrs,id=dc$id,SL.library=mylib)

## ----fit EYax int,cache=TRUE,include=FALSE-------------------------------
# Intervention villages (X=1), fitted curve
fit_i <- agecurveAb(Y=log10(di$ifatpftitre+1),Age=di$ageyrs,id=di$id,SL.library=mylib)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

