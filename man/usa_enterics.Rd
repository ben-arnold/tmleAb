% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/usa_enterics.R
\docType{data}
\name{usa_enterics}
\alias{usa_enterics}
\title{dataset: United States enteric pathogen multiplex antibody data}
\format{A data frame with 86 observations and 10 variables.
\describe{
 \item{\code{id}}{individual ID}
 \item{\code{age}}{age in years. Note that age records years completed, except for children <1 year old.}
 \item{\code{cp17}}{response (MFI-bg) to Cryptosporidium parvum recombinant 17-kDa antigen}
 \item{\code{cp23}}{response (MFI-bg) to Cryptosporidium parvum recombinant 27-kDa antigen}
 \item{\code{vsp5}}{response (MFI-bg) to the VSP-5 fragment of Giardia intestinalis variant-specific surface protein 42e }
 \item{\code{leca}}{response (MFI-bg) to Entamoeba histolytica lectin adhesion molecule (LecA)}
 \item{\code{etec}}{response (MFI-bg) to Enterotoxigenic Escherichia coli (ETEC) heat labile toxin β subunit}
 \item{\code{salb}}{response (MFI-bg) to lipopolysaccharide (LPS) from Salmonella enterica serotype Typhimurium (Group B)}
 \item{\code{norogii}}{response (MFI-bg) to purified recombinant norovirus GII.4 New Orleans virus-like particles from a baculovirus expression system }
 \item{\code{norogi}}{response (MFI-bg) to purified recombinant norovirus GI.4 virus-like particles from a baculovirus expression system }
 }}
\usage{
data(usa_enterics)
}
\description{
IgG antibody responses to protozoan, bacterial, and viral enteric pathogens measured with a multiplex bead assay in sample sets collected from 86 anonymous blood donors in the USA. There is one record per child. All antibody levels were measured on the Luminex platform and are in units of median flourescence intensity minus background (MFI-bg). A comparable dataset (and run with the same bead set) for children in Haiti is also included in the package: \code{\link{haiti_enterics}}.
}
\references{
A. C. McDonald et al., Cryptosporidium parvum-specific antibody responses among children residing in Milwaukee during the 1993 waterborne outbreak. J. Infect. Dis. 183, 1373–1379 (2001).

D. M. Moss et al., Longitudinal Evaluation of Enteric Protozoa in Haitian Children by Stool Exam and Multiplex Serologic Assay.
The American Journal of Tropical Medicine and Hygiene 90 (4): 653–60 (2014).

D. M. Moss, J. M. Montgomery, S. V. Newland, J. W. Priest, P. J. Lammie, Detection of cryptosporidium antibodies in sera and oral fluids using multiplex bead assay. J. Parasitol. 90, 397–404 (2004).

J. W. Priest et al., Multiplex assay detection of immunoglobulin G antibodies that recognize Giardia intestinalis and Cryptosporidium parvum antigens. Clin. Vaccine Immunol. 17, 1695–1707 (2010).

E. Houpt et al., Prevention of intestinal amebiasis by vaccination with the Entamoeba histolytica Gal/GalNac lectin. Vaccine. 22, 611–617 (2004).

J. Flores et al., Enterotoxigenic Escherichia coli Heat-Labile Toxin Seroconversion in US Travelers to Mexico. J. Travel Med. 15, 156–161 (2008).

M. A. Strid, T. Dalby, K. Mølbak, K. A. Krogfelt, Kinetics of the human antibody response against Salmonella enterica Serovars Enteritidis and Typhimurium determined by lipopolysaccharide enzyme-linked immunosorbent assay. Clin. Vaccine Immunol. 14, 741–747 (2007).

E. Vega et al., Novel surveillance network for norovirus gastroenteritis outbreaks, United States. Emerg. Infect. Dis. 17, 1389–1395 (2011).

X. Jiang, M. Wang, D. Y. Graham, M. K. Estes, Expression, self-assembly, and antigenicity of the Norwalk virus capsid protein. J. Virol. 66, 6527–6532 (1992).
}
\seealso{
\code{\link{haiti_enterics}}
}
\keyword{datasets}
