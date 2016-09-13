

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the tmleAb package\nTargeted maximum likelihood estimation for antibody measurements.\n(Version 0.1.1)\n\nCheck github.com/ben-arnold/tmleAb for the most current version\n(we have plans for some useful additions)\n\n")
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
