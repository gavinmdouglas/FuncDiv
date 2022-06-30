rm(list = ls(all.names = TRUE))

library("collapse")

samp_colname <- "sample"
func_colname <- "function."
abun_colname <- "taxon_abun"
taxon_colname <- "taxon"

test_file <- read.table("/Users/gavin/Google Drive/My Drive/postdoc/FuncDiv/isabel_project/data/testfiles/test_contrib_file/hmp_picrust2_v2.4.0-b_path_abun_contrib.tsv.gz",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

test_file <- test_file[, c(samp_colname, func_colname, taxon_colname, abun_colname)]

test_file <- test_file[c(1, 2, 3, 100, 660857, 190002, 190001, 660858, 2000, 190003), ]

test_file <- test_file[collapse::radixorder(test_file[, samp_colname], test_file[, func_colname]), ]

group_starts <- attr(radixorder(test_file[, samp_colname], test_file[, func_colname], starts = TRUE, sort = FALSE), "starts")

sample_set <- test_file[group_starts, samp_colname]
func_set <- test_file[group_starts, func_colname]

abun_values <- list()
num_groups <- length(group_starts)
for (i in 1:(num_groups - 1)) {
  abun_values[[i]] <- test_file[group_starts[i]:(group_starts[i + 1] - 1), abun_colname]
}
abun_values[[i]] <- test_file[group_starts[num_groups]:nrow(test_file), abun_colname]


taxa_labels <- list()
num_groups <- length(group_starts)
for (i in 1:(num_groups - 1)) {
  taxa_labels[[i]] <- test_file[group_starts[i]:(group_starts[i + 1] - 1), taxon_colname]
}
taxa_labels[[i]] <- test_file[group_starts[num_groups]:nrow(test_file), taxon_colname]

