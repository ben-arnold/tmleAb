#' dataset: Mauke island lymphatic filariasis Wb123 antibody measurements
#'
#' Serological measurements from two surveys on Mauke, Cook Islands The dataset includes serum samples from two cross-sectional measurements of the permanent resident population; the first in 1975 (N=362, approximately 58% of the population) and the second in 1992 (N=553, approximately 88% percent of the population). The 1992 measurement occurred five years after a single, island-wide mass drug administration of diethylcarbamazine for all individuals >=5 years old. There were 115 individuals who were measured in both surveys (their \code{id75} value is repeated for both observations). Antibody levels to the Wb123 antigen were measured on the Luciferase Immunoprecipitation System (LIPS) and are in light units.
#'
#' @docType data
#'
#' @usage data(mauke_wb123)
#'
#' @format A data frame with 922 observations and 7 variables.
#' \describe{
#'  \item{\code{id75}}{individual ID in 1975}
#'  \item{\code{id92}}{individual ID in 1992}
#'  \item{\code{year}}{Year of the measurement (1975 or 1992)}
#'  \item{\code{age}}{age in years. Note that age records years completed.}
#'  \item{\code{CAg}}{Level of circulating filarial antigen (>32 is considered positive for circulating antigen)}
#'  \item{\code{wb123}}{response to the Wuchereria bancrofti Wb123 antigen}
#'  \item{\code{wb123pos}}{Seropositive, based on Wb123 using a cutoff of 10968, identified in Kubofcik et al. (2012)}
#'  }
#'
#'
#' @keywords datasets
#'
#'
#' @references C. Steel, J. Kubofcik, E. A. Ottesen, T. B. Nutman, Antibody to the filarial antigen Wb123 reflects reduced transmission and decreased exposure in children born following single mass drug administration (MDA). PLoS Negl. Trop. Dis. 6, e1940 (2012).
#' @references J. Kubofcik, D. L. Fink, T. B. Nutman, Identification of Wb123 as an early and specific marker of Wuchereria bancrofti infection. PLoS Negl. Trop. Dis. 6, e1930 (2012).
#'
#'
"mauke_wb123"
