rm(list = ls(all.names = TRUE))

library("collapse")
library("data.table")

samp_colname <- "sample"
func_colname <- "function."
abun_colname <- "taxon_abun"
taxon_colname <- "taxon"
copy.num_colname <- "genome_function_count"

contrib_tab <- read.table("/home/gdouglas/tmp/hmp_picrust2_v2.4.0-b_path_abun_contrib.tsv.gz",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")


source("/home/gdouglas/scripts/FuncDiv/R/utils.R")

multi_tab <- contrib_to_multitab(contrib_tab = contrib_tab,
                                 samp_colname = "sample",
                                 func_colname = "function.",
                                 abun_colname = "taxon_abun",
                                 taxon_colname = "taxon",
                                 copy.num_colname = "genome_function_count")

contrib_tab_remade <- multitab_to_contrib(func_tab = multi_tab$function_copy_num,
                                          abun_tab = multi_tab$taxon_abun,
                                          ncores = 20,
                                          samp_colname = "sample",
                                          func_colname = "function.",
                                          abun_colname = "taxon_abun",
                                          taxon_colname = "taxon",
                                          copy.num_colname = "genome_function_count")

contrib_tab_orig <- contrib_tab[, colnames(contrib_tab_remade)]

contrib_tab_orig <- contrib_tab_orig[order(contrib_tab_orig$sample, contrib_tab_orig$function., contrib_tab_orig$taxon), ]
contrib_tab_remade <- contrib_tab_remade[order(contrib_tab_remade$sample, contrib_tab_remade$function., contrib_tab_remade$taxon), ]

rownames(contrib_tab_orig) <- NULL
rownames(contrib_tab_remade) <- NULL

all.equal(contrib_tab_remade, contrib_tab_orig)






func_tab_orig <- read.table("/home/gdouglas/scripts/FuncDiv/example_files/func_input.tsv.gz",
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

taxa_tab_orig <- read.table("/home/gdouglas/scripts/FuncDiv/example_files/taxa_input.tsv.gz",
                            header = TRUE, sep = "\t", row.names = 1)

contrib_tab_out <- multitab_to_contrib(func_tab = func_tab_orig,
                                          abun_tab = taxa_tab_orig,
                                          ncores = 20,
                                          samp_colname = "sample",
                                          func_colname = "function.",
                                          abun_colname = "taxon_abun",
                                          taxon_colname = "taxon",
                                          copy.num_colname = "genome_function_count")


contrib_tab_out[which(contrib_tab_out$samp == "ERR1190790" & 
                        contrib_tab_out$func == "K01591"), "taxon"]

rownames(taxa_tab_orig)[which(taxa_tab_orig$ERR1190790 > 0)]





write.table(x = contrib_tab_out, file = "/home/gdouglas/scripts/FuncDiv/example_files/contrib_input.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab_orig,
                                                abun_table = taxa_tab_orig)
func_tab_orig_subset <- subsetted_tables$func
taxa_tab_orig_subset <- subsetted_tables$abun

rm(subsetted_tables)

mutlitab_remade <- contrib_to_multitab(contrib_tab = contrib_tab_out, 
                                       samp_colname = "sample",
                                       func_colname = "function.",
                                       abun_colname = "taxon_abun",
                                       taxon_colname = "taxon",
                                       copy.num_colname = "genome_function_count")


taxon_abun_remade <- mutlitab_remade$taxon_abun[rownames(taxa_tab_orig_subset), colnames(taxa_tab_orig_subset)]
func_tab_remade <- mutlitab_remade$function_copy_num[rownames(func_tab_orig_subset), colnames(func_tab_orig_subset)]

all.equal(taxon_abun_remade, taxa_tab_orig_subset)
all.equal(func_tab_remade, func_tab_orig_subset)


