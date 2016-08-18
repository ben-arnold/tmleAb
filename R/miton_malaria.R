#' DESCRIPTION OF DATA HERE (what it is)
#'
#' DESCRIPTION OF DATA HERE (where its from, what kinds of things at the high-level are contained in it)
#'
#' @docType data
#'
#' @usage data(miton_malaria)
#'
#' @format dataframe
#' each row includes variable descriptions and details for a single individual
#'
#' @keywords datasets
#'
#'        #----------TODO: I put the references I found in the README on OSF here. Edit if you'd like---------------
#'
#' @references Arnold, Benjamin F., Jeffrey W. Priest, Katy L. Hamlin, Delynn M. Moss, John M. Colford Jr, and Patrick J. Lammie. 2014.
#' “Serological Measures of Malaria Transmission in Haiti: Comparison of Longitudinal and Cross-Sectional Methods.” PloS One 9 (4): e93684.
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
"miton_malaria"
