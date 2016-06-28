#' Open Science Framework data on enteric disease in the United States
#'
#' The USA data include measurements from 86 anonymous blood donors in the USA. There is one record per child.
#' The codebooks for each file include variable descriptions and details.
#'
#' @docType data
#'
#' @usage data(usa_enterics)
#'
#' @format dataframe
#'
#' @keywords datasets
#'
#'        #----------TODO: I put the references I found in the README on OSF here. Edit if you'd like---------------
#'
#' @references Moore CR, Johnson LS, Kwak I-Y, Livny M, Broman KW,
#' Moss, Delynn M., Jeffrey W. Priest, Kathy Hamlin, Gordana Derado, Joel Herbein, William A. Petri Jr, and Patrick J. Lammie. 2014.
#' “Longitudinal Evaluation of Enteric Protozoa in Haitian Children by Stool Exam and Multiplex Serologic Assay.”
#' The American Journal of Tropical Medicine and Hygiene 90 (4): 653–60.
#'
#' @source \href{https://osf.io/7frw2/}{Open Science Framework Storage}
#'
#'
#'           #----------TODO: Put examples here ---------------
#' @examples
#' data(grav)
#' times <- attr(grav, "time")
#' phe <- grav$pheno
#' \donttest{
#' iplotCurves(phe, times, phe[,c(61,121)], phe[,c(121,181)],
#'             chartOpts=list(curves_xlab="Time (hours)", curves_ylab="Root tip angle (degrees)",
#'                            scat1_xlab="Angle at 2 hrs", scat1_ylab="Angle at 4 hrs",
#'                            scat2_xlab="Angle at 4 hrs", scat2_ylab="Angle at 6 hrs"))}
"usa_enterics"
