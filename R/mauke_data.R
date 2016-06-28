#' DESCRIPTION OF DATA HERE (what it is)
#'
#' DESCRIPTION OF DATA HERE (where its from, what kinds of things at the high-level are contained in it)
#'
#' @docType data
#'
#' @usage data(maukeWb123)
#'
#' @format dataframe
#' each row includes variable descriptions and details for a single individual
#'
#' @keywords datasets
#'
#'        #----------TODO: I put the references I found in the README on OSF here. Edit if you'd like---------------
#'
#' @references Steel, Cathy, Joseph Kubofcik, Eric A. Ottesen, and Thomas B. Nutman. 2012.
#' “Antibody to the Filarial Antigen Wb123 Reflects Reduced Transmission and Decreased Exposure in Children Born Following Single Mass Drug Administration (MDA).” PLoS Neglected Tropical Diseases 6 (12): e1940.
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
"maukeWb123"
