#' dataset: Garki Project (Nigeria) serology data
#'
#' The \code{garki_sero} dataset includes a semi-processed version of publicly available serology data from the Garki Project, Nigeria. This file includes merged data from the \code{village}, \code{people}, and \code{serology} tables in the original Garki files. The original files are documented through the Garki Project's website (\url{http://garkiproject.nd.edu/}).  These data were downloaded originally from the Swiss Tropical Institute (\url{http://www.swisstph.ch/fr/ressources/epidemiological-databases.html}). The dictionary file from STI included important additional information about data labels. The Stata script that merged the original tables is located here: \url{https://osf.io/7frw2/}.
#'
#' @docType data
#'
#' @usage data(garki_sero)
#'
#' @format A data frame with 20,441 observations and 52 variables. Variable labels below have been copied from \url{http://garkiproject.nd.edu}. The original table's source is added as a prefix, else it was derived.
#' \describe{
#'  \item{\code{id}}{people: individual ID}
#'  \item{\code{village}}{village: village ID}
#'  \item{\code{compount}}{people: compound ID}
#'  \item{\code{serosvy}}{serology: Serological survey number}
#'  \item{\code{survey}}{village: Garki project survey number}
#'  \item{\code{area}}{village: Garki project study area}
#'  \item{\code{grid}}{village: Garki project grid number}
#'  \item{\code{fup}}{village: The follow-up cluster village is in}
#'  \item{\code{vname}}{village: Village name}
#'  \item{\code{vu}}{village: Code representing one village unit for a group}
#'  \item{\code{grouped1}}{village: number of compounds in compact part? CD}
#'  \item{\code{grouped2}}{village: unknown}
#'  \item{\code{scattered1}}{village: number of compounds in scattered part? CD}
#'  \item{\code{scattered2}}{village: unknown}
#'  \item{\code{population}}{village: Population of village (when?) CD}
#'  \item{\code{move}}{people: unknown}
#'  \item{\code{pers}}{people: unknown}
#'  \item{\code{sfreg}}{people: The number of the survey the person first registered with.}
#'  \item{\code{sex}}{people: Identifies the sex of the person. Coding unknown.}
#'  \item{\code{dob}}{people: unknown}
#'  \item{\code{dateb}}{people: The date of birth of the person.}
#'  \item{\code{batch_igg}}{serology: Batch number IgG.}
#'  \item{\code{batch_igm}}{serology: Batch number IgM.}
#'  \item{\code{plate_igg}}{serology: Plate number IgG.}
#'  \item{\code{plate_igm}}{serology: Plate number IgM.}
#'  \item{\code{diam_igg}}{serology: Specimen diameter IgG.}
#'  \item{\code{diam_igm}}{serology: Specimen diameter IgM.}
#'  \item{\code{eval_igg}}{serology: Specimen evaluation IgG.}
#'  \item{\code{eval_igm}}{serology: Specimen evaluation IgM.}
#'  \item{\code{perc_igg}}{serology: Percentage IgG.}
#'  \item{\code{perc_igm}}{serology: Percentage IgM.}
#'  \item{\code{precipitin}}{serology: Result of the precipitin reaction. Coded as follows: 0=Test not done, 1-7=Number of bands, 8=No bands}
#'  \item{\code{blank1}}{serology: unknown}
#'  \item{\code{ifat_pfa}}{serology: IFA test for Plasmodium falciparum, 0=no result, 1= Titre 20, 2=60, 3=180, 4=540, 5=1620,6=4860, 7=14580,8=titres >14580, 9= titre < 20 (for computations take 0); In general Titre=20*3**(code-1) }
#'  \item{\code{blank2}}{serology: unknown}
#'  \item{\code{ifat_pm}}{serology: IFA test for Plasmodium malariae, 0=no result, 1= Titre 20, 2=60, 3=180, 4=540, 5=1620,6=4860, 7=14580,8=titres >14580, 9= titre < 20 (for computations take 0); In general Titre=20*3**(code-1)  }
#'  \item{\code{haem_code}}{serology: 1 = AA; 2 = AS, 3 = SS, 4 =AC, 5 = CC, 6 = SC, 9 = not done}
#'  \item{\code{haem_phbs}}{serology: unknown}
#'  \item{\code{recd_iha}}{serology: Record identifier for IHA (code 116)}
#'  \item{\code{specno}}{serology: specimen number}
#'  \item{\code{titrcode}}{serology: IHA titre code, 4 to 16 = titre 2**4 to 2**16 (in general titre=2**code)}
#'  \item{\code{date}}{serology: unknown}
#'  \item{\code{tr}}{(derived): factor: control or intervention village}
#'  \item{\code{phase}}{(derived): factor: baseline, intervention, post-intervention phase}
#'  \item{\code{wetseason}}{(derived): factor: dry or wet season}
#'  \item{\code{bdate}}{(derived): individuals date of birth (standard date format)}
#'  \item{\code{mdate}}{(derived): measurement date}
#'  \item{\code{agedays}}{(derived): age in days at time of measurement}
#'  \item{\code{ageyrs}}{(derived): age in years at time of measurement}
#'  \item{\code{ifatpftitre}}{(derived): IFA test for P. falciparum antibody titres}
#'  \item{\code{ifatpmtitre}}{(derived): IFA test for P. malariae antibody titres}
#'  \item{\code{ihatitre}}{(derived): IHA test antibody titres}
#' }
#'
#' @keywords datasets
#'
#' @references Molineaux L, Gramiccia G, 1980. The Garki Project. Geneva: World Health Organisation.
#' http://whqlibdoc.who.int/publications/9241560614.pdf
#' @references http://garkiproject.nd.edu/
#' @references
#'
#' @source \url{http://www.swisstph.ch/fr/ressources/epidemiological-databases.html}
#'
#'
"garki_sero"
