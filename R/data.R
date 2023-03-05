#' Function to create Rcpp pointer for Jensen-Shannon divergence function
#'
#' @description
#' ParallelDist can take custom distance metrics in as Rcpp pointers. However,
#' this was causing issues when defined as a normal R object in the FuncDiv
#' package, so instead, the pointer will be created as needed internally with
#' this packaged function.
#' 
#' This also serves as an example of how one can specify other custom pointer
#' functions to compute beta diversity (see below).
#'
#' @format ## `create_jensen_shannon_divergence_FuncPtr()`
#' Function that a returns a RcppXPtrUtils::cppXPtr pointer.
#' Call without parentheses to see the code: `create_jensen_shannon_divergence_FuncPtr`
"create_jensen_shannon_divergence_FuncPtr"
