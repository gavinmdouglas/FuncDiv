# Test that overall workflow works for calculating a variety of metrics.
# Also make sure that the computations work with either the two-table or
# contributional table options.

# Make sure that expected errors are returned under expected conditions as well.

library(FuncDiv)

# For testing:
rm(list = ls(all.names = TRUE))
setwd("/home/gdouglas/scripts/FuncDiv/tests/testthat/")



# Read in test input files.
contrib_tab <- read.table("../../example_files/contrib_input.tsv.gz",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

func_tab <- read.table("../../example_files/func_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

abun_tab <- read.table("../../example_files/taxa_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

# Small sets for quicker tests.
func_tab_subset <- func_tab[1:2, ,]
abun_tab_subset <- abun_tab[, 1:10]

test_tree <- ape::read.tree("../../example_files/taxa.tree")
test_tree_num.added <- test_tree
test_tree_num.added$tip.label <- paste("111", test_tree_num.added$tip.label, sep = "")

K16052_contributors <- colnames(func_tab_subset)[which(func_tab_subset["K16052", ] > 0)]

K16052_contributors_in_SRR4481738 <- K16052_contributors[which(abun_tab_subset[K16052_contributors, "SRR4481738"] > 0)]

K16052_contributors_in_SRR4481738_abun <- abun_tab_subset[K16052_contributors_in_SRR4481738, "SRR4481738"]

abun_alpha_div(metric = "richness", x = K16052_contributors_in_SRR4481738_abun)
abun_alpha_div(metric = "faiths_pd", x = K16052_contributors_in_SRR4481738, tree = test_tree)

picante::pd(samp = t(abun_tab_subset[K16052_contributors_in_SRR4481738, "SRR4481738", drop = FALSE]), tree = test_tree)


metric_out <- c()

for (m in all_alpha_metrics) {

  if (m == "faiths_pd") {
    metric_out <- c(metric_out, abun_alpha_div(metric = m, x = K16052_contributors_in_SRR4481738, tree = test_tree))
  } else {
    metric_out <- c(metric_out, abun_alpha_div(metric = m, x = K16052_contributors_in_SRR4481738_abun))
  }
}