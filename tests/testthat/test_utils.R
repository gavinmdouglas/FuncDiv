# Check functions for converting between multi-table and contributional formatted intables.

library(FuncDiv)

# Read in test input files.
contrib_tab <- read.table("../../example_files/contrib_input.tsv.gz",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
contrib_tab <- contrib_tab[order(contrib_tab$samp, contrib_tab$func, contrib_tab$tax), ]


func_tab <- read.table("../../example_files/func_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
abun_tab <- read.table("../../example_files/taxa_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1)


func_tab_num.added <- func_tab
rownames(func_tab_num.added) <- paste("111", rownames(func_tab_num.added), sep = "")
colnames(func_tab_num.added) <- paste("111", colnames(func_tab_num.added), sep = "")

abun_tab_num.added <- abun_tab
rownames(abun_tab_num.added) <- paste("111", rownames(abun_tab_num.added), sep = "")
colnames(abun_tab_num.added) <- paste("111", colnames(abun_tab_num.added), sep = "")

contrib_tab_num.added <- contrib_tab
contrib_tab_num.added$samp <- paste("111", contrib_tab_num.added$samp, sep = "")
contrib_tab_num.added$func <- paste("111", contrib_tab_num.added$func, sep = "")
contrib_tab_num.added$tax <- paste("111", contrib_tab_num.added$tax, sep = "")


test_that("Converting from contributional to multi-table format works as expected.", {

  multi_tab <- contrib_to_multitab(contrib_tab = contrib_tab,
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
  
  expect_equal(multi_tab, expected)
  
})


test_that("Converting from multi-table to contributional format works as expected.", {
  
  contrib_formatted <- multitab_to_contrib(func_tab = func_tab,
                                           abun_tab = abun_tab,
                                           ncores = 1,
                                           samp_colname = "samp",
                                           func_colname = "func",
                                           taxon_colname = "tax",
                                           abun_colname = "tax_abun",
                                           copy.num_colname = "tax_func_copy_number")
  
  contrib_formatted <- contrib_formatted[order(contrib_formatted$samp, contrib_formatted$func, contrib_formatted$tax), ]
  
  expect_equal(contrib_formatted, contrib_tab)
  
})




test_that("Converting from contributional to multi-table format works as expected, even with taxa/samples/functions with numbers at start.", {
  
  multi_tab <- contrib_to_multitab(contrib_tab = contrib_tab_num.added,
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


test_that("Converting from multi-table to contributional format works as expected, even with taxa/samples/functions with numbers at start.", {
  
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
  
  expect_equal(contrib_formatted, contrib_tab)
  
})



