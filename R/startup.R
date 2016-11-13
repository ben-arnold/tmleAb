

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the tmleAb package\nTargeted maximum likelihood estimation for antibody measurements.\n(Version 0.2.0)\n\nPeriodically check for the latest development version using \ndevtools::install_github('ben-arnold/tmleAb')  \n\nThis software is based on work funded by the National Institute of Allergy and Infectius Diseases grant K01-AI119180\n\n")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/tmleAb",
    devtools.install.args = "",
    devtools.name = "Ben Arnold",
    devtools.desc.author = '"Ben Arnold <benarnold@berkeley.edu> [aut, cre]"',
    devtools.desc.license = "GNU General Public License v3.0",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  invisible()
}
