# Test that overall workflow works for calculating representative beta diversity metrics.
# Also make sure that the computations work with either the two-table or
# contributional table options.

# Make sure that expected errors are returned under expected conditions as well.

test_that("beta_div_contrib returns same values no matter the input table type.", {
  
  multitab_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                                   func_tab = func_tab,
                                   abun_tab = abun_tab,
                                   ncores = 1, 
                                   return_objects = TRUE)
  
  
  contrib_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                                  contrib_tab = contrib_tab,
                                  ncores = 1,
                                  return_objects = TRUE,
                                  samp_colname = "samp",
                                  func_colname = "func",
                                  taxon_colname = "tax",
                                  abun_colname = "tax_abun")
  
  expect_equal(multitab_out, contrib_out)
  
})


test_that("beta_div_contrib returns same values no matter whether saved to RDS or returned as object.", {
  
  obj_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                              func_tab = func_tab,
                              abun_tab = abun_tab,
                              ncores = 1, 
                              return_objects = TRUE)
  
  out_dir <- paste(tempdir(), "output", sep = "/")
  unlink(out_dir, recursive = TRUE)
  beta_div_contrib(metrics = parDist_methods_test_set,
                   func_tab = func_tab,
                   abun_tab = abun_tab,
                   ncores = 1, 
                   write_outfiles = TRUE,
                   outdir = out_dir)
  
  all_matches <- as.integer()
  
  for (m in parDist_methods_test_set) {

    m_counter <- 0
    
    for (func_id in rownames(func_tab)) {
     
      outfile <- paste(out_dir, '/', m, '/', func_id, '.tsv', sep = "")
      
      intab <- read.table(file = outfile, header = TRUE, sep = "\t", check.names = FALSE)
      
      if (ncol(intab) == 1 && nrow(intab) == 0) {
        intab[1, 1] <- NA 
      }
      
      rownames(intab) <- colnames(intab)
      intab[, 1] <- as.numeric(intab[, 1])
  
      if (sum(apply(obj_out[[m]][[func_id]], 2, is.nan)) > 0) {
        tmp <- obj_out[[m]][[func_id]]
        tmp[apply(tmp, 2, is.nan)] <- NA
        
        if (nrow(intab) > 1) {
          intab <- as.data.frame(sapply(intab, as.numeric))
          rownames(intab) <- colnames(intab)
        }
        
        if (all.equal(intab, tmp)) {
          m_counter <- m_counter + 1
        }
        
      } else {

        if (all.equal(intab, obj_out[[m]][[func_id]])) {
          m_counter <- m_counter + 1
        }
  
      }
      
    }
    
    all_matches <- c(all_matches, m_counter)
    
  }
  
  # Delete temporary files.
  unlink(out_dir, recursive = TRUE)
  
  expect_equal(all_matches, rep(x = 50, length(parDist_methods_test_set)))
  
})


test_that("beta_div_contrib returns same values no matter whether input abundances are relabun or not.", {
  
  orig_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                               func_tab = func_tab,
                               abun_tab = abun_tab,
                               ncores = 1, 
                               return_objects = TRUE)
  
  abun_tab_rel <- data.frame(sweep(abun_tab, 2, colSums(abun_tab), '/'), check.names = FALSE)
  
  abun_tab_rel <- abun_tab_rel[, -which(is.nan(colSums(abun_tab_rel)))]
  
  rel_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                              func_tab = func_tab,
                              abun_tab = abun_tab_rel,
                              ncores = 1,
                              return_objects = TRUE)
  
  expect_equal(rel_out, orig_out)
  
})


test_that("beta_div_contrib check dist values for random sample / KO combo #1 output over all tested metrics.", {
  
  all_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                              func_tab = func_tab["K00228", , drop = FALSE],
                              abun_tab = abun_tab[ , c("ERR1190875", "SRR3506419")],
                               ncores = 1, 
                               return_objects = TRUE)
  
  observed_out <- c()
  for (m in parDist_methods_test_set) {
    observed_out <- c(observed_out, all_out[[m]]$K00228["ERR1190875", "SRR3506419"])
  }
  
  expected_out <- c(1.4142135623731, 1, 2, 1.4142135623731, 2, 2, 1.4142135623731, 1,
                    1.5707963267949, 1.4142135623731, NaN, 2, 1, 1.4142135623731, 2, 1,
                    2, 1, 1, 1, 1, 0.5, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, -Inf, 1, 0, 0, 1, 1)
  
  expect_equal(observed_out, expected_out)
  
})



test_that("beta_div_contrib check dist values for random sample / KO combo #2 output over all tested metrics.", {
  
  all_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                              func_tab = func_tab["K01956", , drop = FALSE],
                              abun_tab = abun_tab[ , c("ERR1190790", "ERR2013618")],
                              ncores = 1, 
                              return_objects = TRUE)
  
  observed_out <- c()
  for (m in parDist_methods_test_set) {
    observed_out <- c(observed_out, all_out[[m]]$K01956["ERR1190790", "ERR2013618"])
  }
  
  expected_out <- c(1.39742611897941, 0.991934259881616, 13.7908235504253, 1.41295542219403,
                    13.6254018879072, 0.997656868608457, 0.79162467299016, 0.99595080004354,
                    1.56901783841108, 1.39742611897941, NaN, 1.98386851976323, 0.466577047332018,
                    0.79162467299016, 1.91208791208791, 0.99595080004354, 13.8831953882195,
                    0.991934259881616, 0.928571428571429, 0.875, 0.866666666666667, -0.189245034576083,
                    0.928571428571429, 0.142857142857143, 0.923076923076923, 0.866071428571429,
                    0.0117647058823529, 0.979381443298969, 0.75, 0.866369379043788, 0.133974596215561,
                    0.928571428571429, 0.928571428571429, 0.857142857142857, -0.986732143575568,
                    0.962962962962963, 0, 0, 0.998221512553748, 1)
  
  expect_equal(expected_out, observed_out)
  
})



test_that("beta_div_contrib works correctly with func_ids input, repeating earlier test with slight change", {
  
  all_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                              func_tab = func_tab,
                              abun_tab = abun_tab[ , c("ERR1190790", "ERR2013618")],
                              func_ids = "K01956",
                              ncores = 1, 
                              return_objects = TRUE)
  
  observed_out <- c()
  for (m in parDist_methods_test_set) {
    observed_out <- c(observed_out, all_out[[m]]$K01956["ERR1190790", "ERR2013618"])
  }
  
  expected_out <- c(1.39742611897941, 0.991934259881616, 13.7908235504253, 1.41295542219403,
                    13.6254018879072, 0.997656868608457, 0.79162467299016, 0.99595080004354,
                    1.56901783841108, 1.39742611897941, NaN, 1.98386851976323, 0.466577047332018,
                    0.79162467299016, 1.91208791208791, 0.99595080004354, 13.8831953882195,
                    0.991934259881616, 0.928571428571429, 0.875, 0.866666666666667, -0.189245034576083,
                    0.928571428571429, 0.142857142857143, 0.923076923076923, 0.866071428571429,
                    0.0117647058823529, 0.979381443298969, 0.75, 0.866369379043788, 0.133974596215561,
                    0.928571428571429, 0.928571428571429, 0.857142857142857, -0.986732143575568,
                    0.962962962962963, 0, 0, 0.998221512553748, 1)
  
  expect_equal(expected_out, observed_out)
  
})


test_that("beta_div_contrib works correctly with func_ids input, repeating earlier test with contrib input as well", {
  
  contrib_tab_subset <- contrib_tab[which(contrib_tab$samp %in% c("ERR1190790", "ERR2013618")), ]
  
  all_out <- beta_div_contrib(metrics = parDist_methods_test_set,
                              contrib_tab = contrib_tab_subset,
                              func_ids = "K01956",
                              ncores = 1, 
                              return_objects = TRUE,
                              samp_colname = "samp",
                              func_colname = "func",
                              taxon_colname = "tax",
                              abun_colname = "tax_abun")
  
  observed_out <- c()
  for (m in parDist_methods_test_set) {
    observed_out <- c(observed_out, all_out[[m]]$K01956["ERR1190790", "ERR2013618"])
  }
  
  expected_out <- c(1.39742611897941, 0.991934259881616, 13.7908235504253, 1.41295542219403,
                    13.6254018879072, 0.997656868608457, 0.79162467299016, 0.99595080004354,
                    1.56901783841108, 1.39742611897941, NaN, 1.98386851976323, 0.466577047332018,
                    0.79162467299016, 1.91208791208791, 0.99595080004354, 13.8831953882195,
                    0.991934259881616, 0.928571428571429, 0.875, 0.866666666666667, -0.189245034576083,
                    0.928571428571429, 0.142857142857143, 0.923076923076923, 0.866071428571429,
                    0.0117647058823529, 0.979381443298969, 0.75, 0.866369379043788, 0.133974596215561,
                    0.928571428571429, 0.928571428571429, 0.857142857142857, -0.986732143575568,
                    0.962962962962963, 0, 0, 0.998221512553748, 1)
  
  expect_equal(observed_out, expected_out)
  
})


test_that("beta_div_contrib produces correct weighted UniFrac values", {
  
  contrib_tab_subset <- contrib_tab[which(contrib_tab$samp %in% c("ERR1190790", "ERR1190797", "ERR1190798", "ERR1190806",
                                                                  "ERR1190816", "ERR1305892", "ERR1190946", "SRR5963251",
                                                                  "SRR3466404", "SRR5963379", "SRR1825367", "SRR6257471")), ]
  
  all_out <- beta_div_contrib(metrics = "weighted_unifrac",
                              contrib_tab = contrib_tab_subset,
                              func_ids = "K11070",
                              in_tree = test_tree,
                              ncores = 1, 
                              return_objects = TRUE,
                              samp_colname = "samp",
                              func_colname = "func",
                              taxon_colname = "tax",
                              abun_colname = "tax_abun")

  observed_out <- c(all_out$weighted_unifrac$K11070["ERR1190790", "ERR1190797"],
                    all_out$weighted_unifrac$K11070["ERR1190798", "ERR1190806"],
                    all_out$weighted_unifrac$K11070["ERR1190816", "ERR1305892"],
                    all_out$weighted_unifrac$K11070["ERR1190946", "SRR5963251"],
                    all_out$weighted_unifrac$K11070["SRR3466404", "SRR5963379"],
                    all_out$weighted_unifrac$K11070["SRR1825367", "SRR6257471"])

  expected_out <- c(2.860657, 3.134838, 2.945104, 2.239138, 2.384578, 2.922436)
  
  expect_equal(observed_out, expected_out, tolerance = 1e-06)
  
})



test_that("beta_div_contrib produces correct unweighted UniFrac values", {
  
  contrib_tab_subset <- contrib_tab[which(contrib_tab$samp %in% c("ERR1190790", "ERR1190797", "ERR1190798", "ERR1190806",
                                                                  "ERR1190816", "ERR1305892", "ERR1190946", "SRR5963251",
                                                                  "SRR3466404", "SRR5963379", "SRR1825367", "SRR6257471")), ]
  
  all_out <- beta_div_contrib(metrics = "unweighted_unifrac",
                              contrib_tab = contrib_tab_subset,
                              func_ids = "K11070",
                              in_tree = test_tree,
                              ncores = 1, 
                              return_objects = TRUE,
                              samp_colname = "samp",
                              func_colname = "func",
                              taxon_colname = "tax",
                              abun_colname = "tax_abun")
  
  observed_out <- c(all_out$unweighted_unifrac$K11070["ERR1190790", "ERR1190797"],
                    all_out$unweighted_unifrac$K11070["ERR1190798", "ERR1190806"],
                    all_out$unweighted_unifrac$K11070["ERR1190816", "ERR1305892"],
                    all_out$unweighted_unifrac$K11070["ERR1190946", "SRR5963251"],
                    all_out$unweighted_unifrac$K11070["SRR3466404", "SRR5963379"],
                    all_out$unweighted_unifrac$K11070["SRR1825367", "SRR6257471"])
  
  expected_out <- c(0.7570614, 0.7911921, 0.6419828, 0.9786591, 0.7802325, 0.9739450)
  
  expect_equal(observed_out, expected_out, tolerance = 1e-07)
  
})


test_that("beta_div_contrib produces correct Jensen-Shannon divergence values", {
  
  all_out <- beta_div_contrib(metrics = "jensen_shannon_div",
                              func_tab = func_tab,
                              abun_tab = abun_tab[, c("ERR1190790", "ERR1190797", "ERR1190798", "ERR1190806",
                                                      "ERR1190816", "ERR1305892", "ERR1190946", "SRR5963251",
                                                      "SRR3466404", "SRR5963379", "SRR1825367", "SRR6257471")] + 1,
                              func_ids = "K11070",
                              ncores = 1, 
                              return_objects = TRUE,
                              samp_colname = "samp",
                              func_colname = "func",
                              taxon_colname = "tax",
                              abun_colname = "tax_abun")
  
  observed_out <- c(all_out$jensen_shannon_div$K11070["ERR1190790", "ERR1190797"],
                    all_out$jensen_shannon_div$K11070["ERR1190798", "ERR1190806"],
                    all_out$jensen_shannon_div$K11070["ERR1190816", "ERR1305892"],
                    all_out$jensen_shannon_div$K11070["ERR1190946", "SRR5963251"],
                    all_out$jensen_shannon_div$K11070["SRR3466404", "SRR5963379"],
                    all_out$jensen_shannon_div$K11070["SRR1825367", "SRR6257471"])
  
  expected_out <- c(0.4104305, 0.3295498, 0.3248088, 0.3115956, 0.3159707, 0.2849792)
  
  expect_equal(observed_out, expected_out, tolerance = 1e-7)
  
})


test_that("beta_div_contrib (multi-tab input) returns same values when numbers at start of function, sample, and taxa ids.", {
  
  orig_out <- beta_div_contrib(metrics = c("bray", "binary"),
                               func_tab = func_tab,
                               abun_tab = abun_tab,
                               ncores = 1, 
                               return_objects = TRUE)
  
  
  number_out <- beta_div_contrib(metrics = c("bray", "binary"),
                                 func_tab = func_tab_num.added,
                                 abun_tab = abun_tab_num.added,
                                 ncores = 1, 
                                 return_objects = TRUE)
  
  for (m in names(number_out)) {
    
    names(number_out[[m]]) <- gsub("^111", "", names(number_out[[m]]))
    
    for(func_id in names(number_out[[m]])) {
      colnames(number_out[[m]][[func_id]]) <- gsub("^111", "", colnames(number_out[[m]][[func_id]]))
      rownames(number_out[[m]][[func_id]]) <- gsub("^111", "", rownames(number_out[[m]][[func_id]]))
    }
  }
  
  expect_equal(number_out, orig_out)
  
})


test_that("beta_div_contrib (contrib input) returns same values when numbers at start of function, sample, and taxa ids.", {
  
  orig_out <- beta_div_contrib(metrics = c("bray", "binary"),
                               contrib_tab = contrib_tab,
                               ncores = 1, 
                               return_objects = TRUE,
                               samp_colname = "samp",
                               func_colname = "func",
                               taxon_colname = "tax",
                               abun_colname = "tax_abun")
  
  
  number_out <- beta_div_contrib(metrics = c("bray", "binary"),
                                 contrib_tab = contrib_tab_num.added,
                                 ncores = 1, 
                                 return_objects = TRUE,
                                 samp_colname = "samp",
                                 func_colname = "func",
                                 taxon_colname = "tax",
                                 abun_colname = "tax_abun")
  
  for (m in names(number_out)) {
    
    names(number_out[[m]]) <- gsub("^111", "", names(number_out[[m]]))
    
    for(func_id in names(number_out[[m]])) {
      colnames(number_out[[m]][[func_id]]) <- gsub("^111", "", colnames(number_out[[m]][[func_id]]))
      rownames(number_out[[m]][[func_id]]) <- gsub("^111", "", rownames(number_out[[m]][[func_id]]))
    }
  }
  
  expect_equal(number_out, orig_out)
  
})


test_that("beta_div_contrib returns expected error if none of the input tables provided", {
  
  testthat::expect_error(object = beta_div_contrib(metrics = "binary",
                                                    ncores = 1),
                         regexp = "Stopping - either \"func_tab\" and \"abun_tab\" \\*or\\* \"contrib_tab\" must be specified")
})



test_that("beta_div_contrib returns expected error if only func table input.", {
  
  testthat::expect_error(object = beta_div_contrib(metrics = "binary",
                                                    func_tab = func_tab,
                                                    ncores = 1),
                         regexp = "Stopping - either both \"func_tab\" and \"abun_tab\" or neither must be specified")
})

test_that("beta_div_contrib returns expected error if only abun table input.", {
  
  testthat::expect_error(object = beta_div_contrib(metrics = "binary",
                                                   abun_tab = abun_tab,
                                                   ncores = 1),
                         regexp = "Stopping - either both \"func_tab\" and \"abun_tab\" or neither must be specified")
})



test_that("beta_div_contrib returns expected error if all three input tables specified.", {
  
  testthat::expect_error(object = beta_div_contrib(metrics = "binary",
                                                   abun_tab = abun_tab,
                                                   func_tab = func_tab,
                                                   contrib_tab = contrib_tab,
                                                   ncores = 1),
                         regexp = "Stopping - either \"func_tab\" and \"abun_tab\" \\*or\\* \"contrib_tab\" must be specified, not all!")
})


