#' DESCRIPTION OF DATA HERE (what it is)
#'
#' DESCRIPTION OF DATA HERE (where its from, what kinds of things at the high-level are contained in it)
#'
#' @docType data
#'
#' @usage data(garki_sero)
#'
#' @format dataframe
#'
#' @keywords datasets
#'
#'        #----------TODO: I put the references I found in the README on OSF here. Edit if you'd like---------------
#'
#' @references Molineaux L, Gramiccia G, 1980. The Garki Project. Geneva: World Health Organisation.
#' http://whqlibdoc.who.int/publications/9241560614.pdf
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
"garki_sero"
