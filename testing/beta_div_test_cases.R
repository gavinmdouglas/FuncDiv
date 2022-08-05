# Test that overall workflow works for calculating a variety of metrics.
# Also make sure that the computations work with either the two-table or
# contributional table options.

# Make sure that expected errors are returned under expected conditions as well.

library(FuncDiv)

# For testing:
rm(list = ls(all.names = TRUE))
setwd("/home/gdouglas/scripts/FuncDiv/tests/testthat/")

parDist_methods_test_set <- c("bhjattacharyya", "bray", "canberra", "chord", 
                              "divergence", "dtw", "euclidean", "fJaccard", "geodesic", 
                              "hellinger", "kullback",  "manhattan", 
                              "maximum", "minkowski", "podani", "soergel", "wave", 
                              "whittaker", "binary", "braun-blanquet", "dice", "fager", 
                              "faith", "hamman", "kulczynski1", "kulczynski2", "michael", 
                              "mountford", "mozley", "ochiai", "phi", "russel", "simple matching", 
                              "simpson", "stiles", "tanimoto", "yule", "yule2", "cosine", "hamming")

# Read in test input files.

func_tab <- read.table("../../example_files/func_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

abun_tab <- read.table("../../example_files/taxa_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

all_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                            func_tab = func_tab["K01956", , drop = FALSE],
                            abun_tab = abun_tab[ , c("ERR1190790", "ERR2013618")],
                   ncores = 1, 
                   return_objects = TRUE)

metric_out <- c()
for (m in parDist_methods_test_set) {
  metric_out <- c(metric_out, all_out[[m]]$K01956["ERR1190790", "ERR2013618"])
}


paste(as.character(metric_out), collapse = ", ")




all_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                            func_tab = func_tab["K00228", , drop = FALSE],
                            abun_tab = abun_tab[ , c("ERR1190875", "SRR3506419")],
                            ncores = 1, 
                            return_objects = TRUE)



metric_out <- c()
for (m in parDist_methods_test_set) {
  metric_out <- c(metric_out, all_out[[m]]$K00228["ERR1190875", "SRR3506419"])
}


paste(as.character(metric_out), collapse = ", ")
