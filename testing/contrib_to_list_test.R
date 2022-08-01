rm(list = ls(all.names = TRUE))

library("collapse")
library("data.table")

samp_colname <- "sample"
func_colname <- "function."
abun_colname <- "taxon_abun"
taxon_colname <- "taxon"

contrib_tab <- read.table("/home/gdouglas/tmp/hmp_picrust2_v2.4.0-b_path_abun_contrib.tsv.gz",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

contrib_tab <- contrib_tab[, c(samp_colname, func_colname, taxon_colname, abun_colname)]

contrib_tab_split <- collapse::rsplit(x = contrib_tab, by = contrib_tab$function.)

tmp_in <- ontrib_tab_split[[1]]

tmp <- data.table::dcast.data.table(data = data.table(contrib_tab_split[[1]]),
                                    formula =  contrib_tab_split[[1]]$taxon ~ contrib_tab_split[[1]]$sample,
                                    value.var = 'taxon_abun')

tmp[is.na(tmp)] <- 0

rownames(tmp) <- tmp$contrib_tab_split

tmp <- tmp[, -1]

tmp <- sweep(x = tmp, MARGIN = 2, STATS = colSums(tmp), '/') * 100

contrib_tab <- contrib_tab[collapse::radixorder(contrib_tab[, samp_colname], contrib_tab[, func_colname]), ]

group_starts <- attr(radixorder(contrib_tab[, samp_colname], contrib_tab[, func_colname], starts = TRUE, sort = FALSE), "starts")

sample_set <- contrib_tab[group_starts, samp_colname]
func_set <- contrib_tab[group_starts, func_colname]

abun_values <- list()
num_groups <- length(group_starts)
for (i in 1:(num_groups - 1)) {
  abun_values[[i]] <- contrib_tab[group_starts[i]:(group_starts[i + 1] - 1), abun_colname]
}
abun_values[[num_groups]] <- contrib_tab[group_starts[num_groups]:nrow(contrib_tab), abun_colname]


taxa_labels <- list()
num_groups <- length(group_starts)
for (i in 1:(num_groups - 1)) {
  taxa_labels[[i]] <- contrib_tab[group_starts[i]:(group_starts[i + 1] - 1), taxon_colname]
}
taxa_labels[[i]] <- contrib_tab[group_starts[num_groups]:nrow(contrib_tab), taxon_colname]

