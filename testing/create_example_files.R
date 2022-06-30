rm(list = ls(all.names = TRUE))

library(ape)


# Create example files that should fun fairly quickly.
almeida_ko <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_annot/kegg_summary.tsv.gz",
                         sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)

almeida_abun <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_abun/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- almeida_abun[, -which(colSums(almeida_abun) == 0)]
almeida_abun <- almeida_abun[-which(rowSums(almeida_abun) == 0), ]

intersecting_taxa <- colnames(almeida_ko)[which(colnames(almeida_ko) %in% rownames(almeida_abun))]

almeida_ko <- almeida_ko[, intersecting_taxa]

almeida_abun <- almeida_abun[intersecting_taxa, ]

almeida_tree <- read.tree("/data1/gdouglas/projects/contrib_div/data/Almeida_2019/phylogenies/raxml_hgr-umgs_phylogeny.nwk")

almeida_tree <- phangorn::midpoint(almeida_tree)

top_100_samples <- names(sort(colSums(almeida_abun), decreasing = TRUE))[1:100]

top_100_samples_subset <- almeida_abun[, top_100_samples]
top_100_samples_subset <- top_100_samples_subset[-which(rowSums(top_100_samples_subset) == 0), ]
top_100_samples_subset <- top_100_samples_subset[sample(1:nrow(top_100_samples_subset), size = 200), ]

almeida_ko_subset <- almeida_ko[, rownames(top_100_samples_subset)]
almeida_ko_subset <- almeida_ko_subset[-which(rowSums(almeida_ko_subset) == 0), ]
almeida_ko_subset <- almeida_ko_subset[sample(1:nrow(almeida_ko_subset), size = 50), ]

almeida_tree_subset <- drop.tip(phy = almeida_tree, tip = almeida_tree$tip.label[which(! almeida_tree$tip.label %in% rownames(top_100_samples_subset))])


write.table(x = almeida_ko_subset, file = "/home/gdouglas/scripts/FuncDiv/example_files/func_input.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = top_100_samples_subset, file = "/home/gdouglas/scripts/FuncDiv/example_files/taxa_input.tsv",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.tree(phy = almeida_tree_subset, file = "/home/gdouglas/scripts/FuncDiv/example_files/taxa.tree")


