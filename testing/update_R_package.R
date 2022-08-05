### Code for re-generating R package files whenever there has been additional
### functions added or changes to man pages.

rm(list = ls(all.names = TRUE))

library(devtools)
library(roxygen2)

setwd("/home/gdouglas/scripts/FuncDiv/")

# NAMESPACE should be deleted if want new lines to be added.
# Otherwise it wont be changed, since manual lines were added.

usethis::use_rcpp()
usethis::use_testthat()
devtools::document()

# Note that I added these lines manually after the NAMESPACE was generated:
# useDynLib(FuncDiv)
# exportPattern("ˆ[[:alpha:]]+")
# importFrom(Rcpp, evalCpp)


devtools::load_all(path = "/home/gdouglas/scripts/FuncDiv/")
devtools::test()
#devtools::build()


