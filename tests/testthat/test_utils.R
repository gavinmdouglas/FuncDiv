# Check functions for converting between multi-table and contributional formatted intables.

test_that("converting from contributional to multi-table format works as expected.", {

  multi_tab <- contrib_to_multitab(contrib_tab = contrib_tab_ordered,
                                   samp_colname = "samp",
                                   func_colname = "func",
                                   taxon_colname = "tax",
                                   abun_colname = "tax_abun",
                                   copy.num_colname = "tax_func_copy_number")
  
  # Note that the subset_func_and_abun_tables function is also tested as it is included here.
  subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab,
                                                  abun_table = abun_tab)
  
  expected <- list(taxon_abun = subsetted_tables$abun[sort(rownames(subsetted_tables$abun)),
                                                      sort(colnames(subsetted_tables$abun))],
                   function_copy_num = subsetted_tables$func[sort(rownames(subsetted_tables$func)),
                                                             sort(colnames(subsetted_tables$func))])
  
  expect_equal(multi_tab, expected)
  
})


test_that("converting from multi-table to contributional format works as expected.", {
  
  contrib_formatted <- multitab_to_contrib(func_tab = func_tab,
                                           abun_tab = abun_tab,
                                           ncores = 1,
                                           samp_colname = "samp",
                                           func_colname = "func",
                                           taxon_colname = "tax",
                                           abun_colname = "tax_abun",
                                           copy.num_colname = "tax_func_copy_number")
  
  contrib_formatted <- contrib_formatted[order(contrib_formatted$samp, contrib_formatted$func, contrib_formatted$tax), ]
  
  expect_equal(contrib_formatted, contrib_tab_ordered)
  
})




test_that("converting from contributional to multi-table format works as expected, even with taxa/samples/functions with numbers at start.", {
  
  multi_tab <- contrib_to_multitab(contrib_tab = contrib_tab_ordered_num.added,
                                   samp_colname = "samp",
                                   func_colname = "func",
                                   taxon_colname = "tax",
                                   abun_colname = "tax_abun",
                                   copy.num_colname = "tax_func_copy_number")
  
  
  subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab,
                                                  abun_table = abun_tab)
  
  expected <- list(taxon_abun = subsetted_tables$abun[sort(rownames(subsetted_tables$abun)),
                                                      sort(colnames(subsetted_tables$abun))],
                   function_copy_num = subsetted_tables$func[sort(rownames(subsetted_tables$func)),
                                                             sort(colnames(subsetted_tables$func))])
  
  for (m in names(multi_tab)) {
    colnames(multi_tab[[m]]) <- gsub("^111", "", colnames(multi_tab[[m]]))
    rownames(multi_tab[[m]]) <- gsub("^111", "", rownames(multi_tab[[m]]))
  }
  
  expect_equal(multi_tab, expected)
  
})


test_that("converting from multi-table to contributional format works as expected, even with taxa/samples/functions with numbers at start", {

  contrib_formatted <- multitab_to_contrib(func_tab = func_tab_num.added,
                                           abun_tab = abun_tab_num.added,
                                           ncores = 1,
                                           samp_colname = "samp",
                                           func_colname = "func",
                                           taxon_colname = "tax",
                                           abun_colname = "tax_abun",
                                           copy.num_colname = "tax_func_copy_number")
  
  contrib_formatted$samp <- gsub("^111", "", contrib_formatted$samp)
  contrib_formatted$func <- gsub("^111", "", contrib_formatted$func)
  contrib_formatted$tax <- gsub("^111", "", contrib_formatted$tax)
  
  contrib_formatted <- contrib_formatted[order(contrib_formatted$samp, contrib_formatted$func, contrib_formatted$tax), ]
  
  expect_equal(contrib_formatted, contrib_tab_ordered)
  
})


test_that("table subsetting function works as expected when specific func_id specified", {
  
  subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab,
                                                  abun_table = abun_tab,
                                                  func_ids = c("K04749"))
  
  exp_out <- list(func = func_tab["K04749", rownames(abun_tab)], abun = abun_tab[, which(colSums(abun_tab) > 0)] )

  expect_equal(subsetted_tables, exp_out)
  
})


test_that("cross-product of taxa abundance and function tables works as expected", {
  
  func_test <- func_tab[1:10, 1:10]
  abun_test <- abun_tab[colnames(func_tab)[1:10], 1:2]
  
  obs_out <- func_abun_crossproduct(func_test, abun_test)
  
  exp_out <- data.frame(matrix(0.0, nrow = 10, ncol = 2))
  colnames(exp_out) <- colnames(abun_test)
  rownames(exp_out) <- rownames(func_test)
  
  for (func_id in rownames(func_test)) {
    
    func_contributors <- colnames(func_test)[which(as.numeric(func_test[func_id, ]) > 0)]
    
    for (samp in colnames(abun_test)) {
      
      exp_out[func_id, samp] <- sum(abun_test[func_contributors, samp])
      
    }

  }
  
  expect_equal(obs_out, exp_out)
  
})
