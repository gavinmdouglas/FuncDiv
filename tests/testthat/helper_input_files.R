# Alpha metrics test objects.
# Toy example input vectors.
in_vec <- c(3, 3, 1, 1, 4, 5)
in_vec_complex <- c(NA, 3, 0, 1, 4, 0)

# Test tree
test_tree <- ape::read.tree("../example_files/taxa.tree")

# Alpha workflow test objects.
# Read in test input files.
contrib_tab <- read.table("../example_files/contrib_input.tsv.gz",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

func_tab <- read.table("../example_files/func_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

abun_tab <- read.table("../example_files/taxa_input.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

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


# Additional beta workflow test objects.
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

parDist_methods_test_set <- c("bhjattacharyya", "bray", "canberra", "chord", 
                              "divergence", "dtw", "euclidean", "fJaccard", "geodesic", 
                              "hellinger", "kullback",  "manhattan", 
                              "maximum", "minkowski", "podani", "soergel", "wave", 
                              "whittaker", "binary", "braun-blanquet", "dice", "fager", 
                              "faith", "hamman", "kulczynski1", "kulczynski2", "michael", 
                              "mountford", "mozley", "ochiai", "phi", "russel", "simple matching", 
                              "simpson", "stiles", "tanimoto", "yule", "yule2", "cosine", "hamming")
# Note: excluded mahalanobis as it caused this error: "inv(): matrix is singular"


# Then for utils tests (note that this table is ordered differently than above).
contrib_tab_ordered <- contrib_tab[order(contrib_tab$samp, contrib_tab$func, contrib_tab$tax), ]
contrib_tab_ordered_num.added <- contrib_tab_ordered
contrib_tab_ordered_num.added$samp <- paste("111", contrib_tab_ordered_num.added$samp, sep = "")
contrib_tab_ordered_num.added$func <- paste("111", contrib_tab_ordered_num.added$func, sep = "")
contrib_tab_ordered_num.added$tax <- paste("111", contrib_tab_ordered_num.added$tax, sep = "")
