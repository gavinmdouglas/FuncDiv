# Should use slightly different options to make sure nothing got hard-coded

# func_tab <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_annot/kegg_summary.tsv.gz",
#                        sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
# abun_tab <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_abun/bwa_depth_min25coverage.tsv.gz",
#                        header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")
# 
# metrics = c("binary", "kullback")
# func_ids = c("K03298", "K14698", "K03895")
# return_objects = FALSE
# write_outfiles = TRUE
# outdir = "/home/gdouglas/tmp/test4"
# parDist_func = NULL
# ncores <- 20

# 
# beta_div_contrib <- function(func_tab,
#                              abun_tab,
#                              metrics = c("binary", "kullback"),
#                              func_ids = NULL,
#                              return_objects = FALSE,
#                              write_outfiles = FALSE,
#                              outdir = NULL,
#                              ncores = 1,
#                              parDist_func = NULL) 