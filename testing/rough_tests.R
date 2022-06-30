rm(list = ls(all.names = TRUE))

setwd("/home/gdouglas/scripts/FuncDiv/testing/")

test_tree <- ape::read.tree("../example_files/taxa.tree")
func_input <- read.table("../example_files/func_input.tsv.gz", row.names = 1, header = TRUE, sep = "\t")
taxa_input <- read.table("../example_files/taxa_input.tsv.gz", row.names = 1, header = TRUE, sep = "\t")


# Alpha workflow.

FuncDiv_abun_alpha_metrics[["TEST"]] <- function(x) { length(x) ** 2 }

alpha_out <- alpha_div_contrib(metrics = c("richness", "shannon_index", "TEST"),
                               func_tab = func_input,
                               abun_tab = taxa_input,
                               replace_NA = TRUE,
                               custom_metric_functions = FuncDiv_abun_alpha_metrics)

CRC_taxa <- readRDS("/home/gdouglas/tmp/2021-10-14.YachidaS_2019_metaphlan2.rds")


# Long-format input

combined_tables <- readRDS(file = "/data1/gdouglas/projects/contrib_div/data/curatedMetagenomicData_ML/input/pathway_tables.rds")


