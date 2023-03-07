faiths_pd <- function(tips_in_sample, tree) {

  # Simplified version of equivalent functions from picante.
  # Requires tree to be rooted and will only compute distances
  # based on edges in tree after pruning to tips in samples
  # (i.e., distance to overall root is not included).

  # NOTE: Since this function is expected to be run with many different input vectors and the same tree,
  # this function does not perform checks on the input tree itself!

  if (length(tips_in_sample) == 1) { return(0) }

  if (length(which(! tips_in_sample %in% tree$tip.label) > 0)) {
     stop("Stopping - some features in sample are not found as tips in the tree.")
  }

  if (length(which(! tree$tip.label %in% tips_in_sample) > 0)) {
    tree <- ape::drop.tip(phy = tree,
                          tip = tree$tip.label[which(! tree$tip.label %in% tips_in_sample)],
                          trim.internal = TRUE)
  }

  return(sum(tree$edge.length))

}

richness <- function(x) {
  length(which(x > 0))
}

shannon_index <- function(x, log_base=2) {
  x <- x[which(x > 0)]
  x <- x / sum(x)
  return(-1 * sum(x * log(x, log_base)))
}

berger_parker_dominance <- function(x) {
  max(x, na.rm = TRUE) / sum(x, na.rm = TRUE)
}

ENS_pie <- function(x) {
  x <- x / sum(x, na.rm = TRUE)
  return(1 / sum(x ** 2, na.rm = TRUE))
}

fishers_alpha <- function(x, min_unique = 10) {

  # Simplified version of function implemented in vegan.
  # Returns NA if fit cannot be performed
  # (which is assumed to always be the case when there are
  # fewer than 'min_unique' unique observations, which is 10 by default).

  x <- x[which(x > 0)]

  if (length(x) < min_unique) { return(NA) }

  unique_obs <- length(x)
  total_obs <- sum(x)

  func2scan <- function(x, unique_obs, total_obs) { x * log(1 + total_obs / x) - unique_obs }

  uniroot_out <- try(stats::uniroot(func2scan,
                                    c(1, 50),
                                    extendInt = "upX",
                                    unique_obs = unique_obs, 
                                    total_obs = total_obs)$root,
                     silent = TRUE)

  if (inherits(uniroot_out, "try-error")) {
    return(NA)
  } else {
    return(uniroot_out)
  }
}

heips_evenness <- function(x) {
  richness_out <- richness(x)
  if (richness_out <= 1) { return(0) }
  (exp(shannon_index(x = x, log_base = exp(1))) - 1) / (richness_out - 1)
}

margalefs_richness <- function(x) {
  (richness(x) - 1) / log(sum(x, na.rm = TRUE))
}

mcintoshs_dominance <- function(x) {
  total_obs <- sum(x, na.rm = TRUE)
  (total_obs - sqrt(sum(x ** 2, na.rm = TRUE))) / (total_obs - sqrt(total_obs))
}

mcintoshs_evenness <- function(x) {
  x <- x[which(x > 0)]
  total_obs <- sum(x, na.rm = TRUE)
  unique_obs <- length(x)
  sqrt(sum(x ** 2, na.rm = TRUE)) / sqrt((total_obs - unique_obs + 1)**2 + unique_obs - 1)
}

menhinicks_richness <- function(x) {
  x <- x[which(x > 0)]
  length(x) / sqrt(sum(x, na.rm = TRUE))
}

pielous_evenness <- function(x) {
  x <- x[which(x > 0)]
  if (length(x) == 1) { return(0) }
  shannon_index(x = x, log_base = exp(1)) / log(length(x))
}

gini_simpson_index <- function(x) {
  x <- x / sum(x, na.rm = TRUE)
  return(1 - sum(x ** 2, na.rm = TRUE))
}

simpsons_evenness <- function(x) {
  x <- x / sum(x, na.rm = TRUE)
  return(inverse_simpson_index(x) / richness(x))
}

inverse_simpson_index <- function(x) {
  x <- x / sum(x, na.rm = TRUE)
  return(1 / sum(x ** 2, na.rm = TRUE))
}


#' List object containing the functions to compute the default alpha diversity metrics
#' 
#' These functions are used by the `alpha_div_contrib` function to compute contributional diversity,
#' but can be used for any arbitrary input vector as well to compute standard alpha diversity.
#' 
#' The metrics were primarily taken from definitions provided by `scikit-bio` Python package, as well as the `vegan` and `picante` R packages.
#' The functions are provided as elements of this list, so that it is more convenient to call them programatically.
#' All available alpha diversity metrics can be seen by typing `names(FuncDiv_alpha_metrics)`.
#' The code to compute each metric can be inspected for each function, for instance, for richness, by typing: `FuncDiv_alpha_metrics$richness`.
#'
#' These functions all have a single input: a numeric vector containing taxa abundances within a given sample.
#' The exception is for `faiths_pd`, which expects a character vector of taxa labels that are present, as well as a tree (phylo object),
#' which must contain all these specified taxa labels as tip labels.
#'
#' Note that not all these metrics are appropriate for relative abundance data. In particular, these metrics expect count data (e.g., read counts) 
#' corresponding to the number of occurrences of each category (e.g., each microbe): `menhinicks_richness`, `mcintoshs_evenness`, `mcintoshs_dominance`,
#' `margalefs_richness`, and `fishers_alpha`.
#' 
#' @return numeric vector with alpha diversity value.
#'
#' @examples
#' # Most metrics just require an input vector of abundances.
#' test_abun <- c(0, NA, 1, 2, 10, 4)
#' FuncDiv_alpha_metrics[["richness"]](test_abun)
#'
#' # Note that the input for computing Faith's PD is different.
#' # Get a randomly generated tree:
#' test_tree <- ape::rtree(n = 50)
#' test_present_tips <- c('t1', 't2', 't3')
#' FuncDiv_alpha_metrics[["faiths_pd"]](test_present_tips, test_tree)
#'
#' @export
FuncDiv_alpha_metrics <- list()
FuncDiv_alpha_metrics[["richness"]] <- richness
FuncDiv_alpha_metrics[["shannon_index"]] <- shannon_index
FuncDiv_alpha_metrics[["berger_parker_dominance"]] <- berger_parker_dominance
FuncDiv_alpha_metrics[["ENS_pie"]] <- ENS_pie
FuncDiv_alpha_metrics[["faiths_pd"]] <- faiths_pd
FuncDiv_alpha_metrics[["fishers_alpha"]] <- fishers_alpha
FuncDiv_alpha_metrics[["heips_evenness"]] <- heips_evenness
FuncDiv_alpha_metrics[["margalefs_richness"]] <- margalefs_richness
FuncDiv_alpha_metrics[["mcintoshs_dominance"]] <- mcintoshs_dominance
FuncDiv_alpha_metrics[["mcintoshs_evenness"]] <- mcintoshs_evenness
FuncDiv_alpha_metrics[["menhinicks_richness"]] <- menhinicks_richness
FuncDiv_alpha_metrics[["pielous_evenness"]] <- pielous_evenness
FuncDiv_alpha_metrics[["gini_simpson_index"]] <- gini_simpson_index
FuncDiv_alpha_metrics[["simpsons_evenness"]] <- simpsons_evenness
FuncDiv_alpha_metrics[["inverse_simpson_index"]] <- inverse_simpson_index

#' Convenience function for running default alpha diversity metrics on a single vector input
#' 
#' This is a simple wrapper for `FuncDiv_alpha_metrics`, and you can see more details with `?FuncDiv_alpha_metrics`.
#' 
#' These functions all have a single input: a numeric vector containing taxa abundances within a given sample.
#' The exception is for `faiths_pd`, which expects a character vector of taxa labels that are present, as well as a tree (phylo object),
#' which must contain all these specified taxa labels as tip labels.
#' 
#' @param x input vector. Either class numeric (representing abundance of categories \[e.g., microbes\]) or character (indicating which taxa are present, which is required for `faiths_pd`). 
#' @param metric alpha diversity metric to compute. Must be one of `names(FuncDiv_alpha_metrics)`.
#' @param ... included so that functions with single arguments will not throw errors if `tree` is included (and ignored). This should be a phylo object to use in case of `faiths_pd`. 
#'
#' @return numeric vector with alpha diversity value.
#'
#' @examples
#' # Most metrics just require an input vector of abundances.
#' test_abun <- c(0, NA, 1, 2, 10, 4)
#' compute_alpha_div(x = test_abun, metric = "richness")
#'
#' # Note that the input for computing Faith's PD is different.
#' # Get a randomly generated tree:
#' test_tree <- ape::rtree(n = 50)
#' test_present_tips <- c('t1', 't2', 't3')
#' compute_alpha_div(x = test_present_tips, metric = "faiths_pd", tree = test_tree)
#' 
#' @export
compute_alpha_div <- function(x, metric, ...) {
  FuncDiv_alpha_metrics[[metric]](x, ...)
}
