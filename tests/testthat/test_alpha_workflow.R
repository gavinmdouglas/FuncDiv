# Test that overall workflow works for calculating a variety of metrics.
# Also make sure that the computations work with either the two-table or
# contributional table options.

# Make sure that expected errors are returned under expected conditions as well.

library(FuncDiv)

# Read in test input files.
contrib_tab <- read.table("../../example_files/contrib_input.tsv.gz",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

func_tab <- read.table("../../example_files/func_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

abun_tab <- read.table("../../example_files/taxa_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

# Small sets for quicker tests.
func_tab_subset <- func_tab[1:2, ,]
abun_tab_subset <- abun_tab[, 1:10]

contrib_tab_subset <- contrib_tab
contrib_tab_subset <- contrib_tab_subset[which(contrib_tab_subset$samp %in% colnames(abun_tab_subset)), ]
contrib_tab_subset <- contrib_tab_subset[which(contrib_tab_subset$func %in% rownames(func_tab_subset)), ]


func_tab_subset_num.added <- func_tab_subset
rownames(func_tab_subset_num.added) <- paste("111", rownames(func_tab_subset_num.added), sep = "")
colnames(func_tab_subset_num.added) <- paste("111", colnames(func_tab_subset_num.added), sep = "")

abun_tab_subset_num.added <- abun_tab_subset
rownames(abun_tab_subset_num.added) <- paste("111", rownames(abun_tab_subset_num.added), sep = "")
colnames(abun_tab_subset_num.added) <- paste("111", colnames(abun_tab_subset_num.added), sep = "")

contrib_tab_subset_num.added <- contrib_tab_subset
contrib_tab_subset_num.added$samp <- paste("111", contrib_tab_subset_num.added$samp, sep = "")
contrib_tab_subset_num.added$func <- paste("111", contrib_tab_subset_num.added$func, sep = "")
contrib_tab_subset_num.added$tax <- paste("111", contrib_tab_subset_num.added$tax, sep = "")

# Tree needed for calculating Faith's PD.
test_tree <- ape::read.tree("../../example_files/taxa.tree")
test_tree_num.added <- test_tree
test_tree_num.added$tip.label <- paste("111", test_tree_num.added$tip.label, sep = "")

# Re-make richness as a custom metric
custom_metric_functions <- list()
custom_metric_functions[["test"]] <- function(x) { length(which(x > 0)) }

all_alpha_metrics <- c("richness", "shannon_index", "berger_parker_dominance", "ENS_pie",
                       "faiths_pd", "fishers_alpha", "heips_evenness", "margalefs_richness",
                       "mcintoshs_dominance", "mcintoshs_evenness", "menhinicks_richness",
                       "pielous_evenness", "gini_simpson_index", "simpsons_evenness",
                       "inverse_simpson_index")

non_relabun_metrics <- c("menhinicks_richness", "mcintoshs_evenness", "mcintoshs_dominance",
                         "margalefs_richness", "fishers_alpha")

test_that("alpha_div_contrib returns same values with either multi table or contributional table input", {
  
  multitab_out <- alpha_div_contrib(metrics = c(all_alpha_metrics),
                                          func_tab = func_tab,
                                          abun_tab = abun_tab,
                                          in_tree = test_tree,
                                          ncores = 1)

  contrib_out <- alpha_div_contrib(metrics = c(all_alpha_metrics),
                                         contrib_tab = contrib_tab,
                                         in_tree = test_tree,
                                         ncores = 1,
                                         samp_colname = "samp",
                                         func_colname = "func",
                                         taxon_colname = "tax",
                                         abun_colname = "tax_abun")
  
  expect_equal(multitab_out, contrib_out)
  
})


test_that("alpha_div_contrib (multi-tab input) returns same values when abun is relative or absolute", {
  
  relabun_compatible_metrics_only <- all_alpha_metrics[which(! all_alpha_metrics %in% non_relabun_metrics)]
    
  
  nonrel_out <- alpha_div_contrib(metrics = relabun_compatible_metrics_only,
                                  func_tab = func_tab_subset,
                                  abun_tab = abun_tab_subset,
                                  in_tree = test_tree,
                                  ncores = 1)
  
  abun_tab_subset_rel <- data.frame(sweep(abun_tab_subset, 2, colSums(abun_tab_subset), '/'), check.names = FALSE)
  
  rel_out <- alpha_div_contrib(metrics = relabun_compatible_metrics_only,
                                  func_tab = func_tab_subset,
                                  abun_tab = abun_tab_subset_rel,
                                  in_tree = test_tree,
                                  ncores = 1)
  
  expect_equal(nonrel_out, rel_out)
  
})


test_that("alpha_div_contrib (contrib input) returns same values when abun is relative or absolute", {
  
  relabun_compatible_metrics_only <- all_alpha_metrics[which(! all_alpha_metrics %in% non_relabun_metrics)]
  
  nonrel_out <- alpha_div_contrib(metrics = relabun_compatible_metrics_only,
                                  func_tab = NULL,
                                  abun_tab = NULL,
                                  contrib_tab = contrib_tab_subset,
                                  in_tree = test_tree,
                                  ncores = 1,
                                  samp_colname = "samp",
                                  func_colname = "func",
                                  taxon_colname = "tax",
                                  abun_colname = "tax_abun")
  
  abun_tab_subset_rel <- data.frame(sweep(abun_tab_subset, 2, colSums(abun_tab_subset), '/'), check.names = FALSE)
  
  contrib_rel <- multitab_to_contrib(func_tab = func_tab_subset,
                                     abun_tab_subset_rel,
                                     samp_colname = "samp",
                                     func_colname = "func",
                                     taxon_colname = "tax",
                                     abun_colname = "tax_abun")

  
  rel_out <- alpha_div_contrib(metrics = relabun_compatible_metrics_only,
                               contrib_tab = contrib_rel,
                               in_tree = test_tree,
                               ncores = 1,
                               samp_colname = "samp",
                               func_colname = "func",
                               taxon_colname = "tax",
                               abun_colname = "tax_abun")
  
  expect_equal(nonrel_out, rel_out)
  
})




test_that("alpha_div_contrib (multi-tab input) returns same values when numbers at start of function, sample, and taxa ids.", {
  
  orig_out <- alpha_div_contrib(metrics = all_alpha_metrics,
                                  func_tab = func_tab_subset,
                                  abun_tab = abun_tab_subset,
                                  in_tree = test_tree,
                                  ncores = 1)
  

  number_out <- alpha_div_contrib(metrics = all_alpha_metrics,
                               func_tab = func_tab_subset_num.added,
                               abun_tab = abun_tab_subset_num.added,
                               in_tree = test_tree_num.added,
                               ncores = 1)
  
  for (m in names(number_out)) {
    colnames(number_out[[m]]) <- gsub("^111", "", colnames(number_out[[m]]))
    rownames(number_out[[m]]) <- gsub("^111", "", rownames(number_out[[m]]))
  }
  
  expect_equal(orig_out, number_out)
  
})


test_that("alpha_div_contrib (contrib input) returns same values when numbers at start of function, sample, and taxa ids.", {
  
  orig_out <- alpha_div_contrib(metrics = all_alpha_metrics,
                                contrib_tab = contrib_tab_subset,
                                in_tree = test_tree,
                                ncores = 1,
                                samp_colname = "samp",
                                func_colname = "func",
                                taxon_colname = "tax",
                                abun_colname = "tax_abun")
  
  
  number_out <- alpha_div_contrib(metrics = all_alpha_metrics,
                                  contrib_tab = contrib_tab_subset_num.added,
                                  in_tree = test_tree_num.added,
                                  ncores = 1,
                                  samp_colname = "samp",
                                  func_colname = "func",
                                  taxon_colname = "tax",
                                  abun_colname = "tax_abun")
  
  for (m in names(number_out)) {
    colnames(number_out[[m]]) <- gsub("^111", "", colnames(number_out[[m]]))
    rownames(number_out[[m]]) <- gsub("^111", "", rownames(number_out[[m]]))
  }
  
  expect_equal(orig_out, number_out)
  
})



test_that("alpha_div_contrib returns expected error if none of the input tables provided", {

  testthat::expect_error(object = alpha_div_contrib(metrics = all_alpha_metrics,
                                                    in_tree = test_tree,
                                                    ncores = 1),
                         regexp = "Stopping - either \"func_tab\" and \"abun_tab\" \\*or\\* \"contrib_tab\" must be specified")
})



test_that("alpha_div_contrib returns expected error if only func table input.", {
  
  testthat::expect_error(object = alpha_div_contrib(metrics = all_alpha_metrics,
                                                    func_tab = func_tab_subset_num.added,
                                                    in_tree = test_tree,
                                                    ncores = 1),
                         regexp = "Stopping - either both \"func_tab\" and \"abun_tab\" or neither must be specified")
})

test_that("alpha_div_contrib returns expected error if only abun table input.", {
  
  testthat::expect_error(object = alpha_div_contrib(metrics = all_alpha_metrics,
                                                    abun_tab = abun_tab_subset_num.added,
                                                    in_tree = test_tree,
                                                    ncores = 1),
                         regexp = "Stopping - either both \"func_tab\" and \"abun_tab\" or neither must be specified")
})



test_that("alpha_div_contrib returns expected error if all three input tables specified.", {
  
  testthat::expect_error(object = alpha_div_contrib(metrics = all_alpha_metrics,
                                                    abun_tab = abun_tab_subset_num.added,
                                                    func_tab = func_tab_subset_num.added,
                                                    contrib_tab = contrib_tab,
                                                    in_tree = test_tree,
                                                    ncores = 1),
                         regexp = "Stopping - either \"func_tab\" and \"abun_tab\" \\*or\\* \"contrib_tab\" must be specified, not all!")
})


test_that("alpha_div_contrib works with custom alpha diversity metric.", {
  
  custom_out <- alpha_div_contrib(metrics = "test",
                                           abun_tab = abun_tab_subset,
                                           func_tab = func_tab_subset,
                                           ncores = 1,
                                           custom_metric_functions = custom_metric_functions)
  
  orig_out <- alpha_div_contrib(metrics = "richness",
                                  abun_tab = abun_tab_subset,
                                  func_tab = func_tab_subset,
                                  ncores = 1)
  
  custom_out$test[is.na(custom_out$test)] <- 0
  
  expect_equal(custom_out$test, orig_out$richness)
  
})
