
#' @export
alpha_div_contrib <- function(metrics, func_tab, abun_tab) {
  
  abun_tab <- abun_tab[which(rowSums(abun_tab) > 0), ]
  abun_tab <- abun_tab[, which(colSums(abun_tab) > 0)]
  
  func_tab <- func_tab[, rownames(abun_tab), drop = FALSE]
  func_tab <- func_tab[which(rowSums(func_tab) > 0), ]
  func_tab <- func_tab[, which(colSums(func_tab) > 0)]
  
  abun_tab <- abun_tab[colnames(func_tab), ]
  
  func_set <- as.character()
  sample_set <- as.character()
  for (i in 1:nrow(func_tab)) {
    func_set <- c(func_set, rep(rownames(func_tab)[i], ncol(abun_tab)))
    sample_set <- c(sample_set, colnames(abun_tab))
  }
  
  div_metric_out <- list()
  
  prepped_abun <- prep_all_sample_func_vec(abun_tab = as.matrix(abun_tab),
                                           func_tab = as.matrix(func_tab))
  
  # Get list of taxa ids that are present in each sample / func combination if faiths_pd specified as well.
  if ("faiths_pd" %in% metrics) {
    prepped_taxa_present <- prep_all_sample_func_taxa_vec(abun_tab = as.matrix(abun_tab),
                                                          func_tab = as.matrix(func_tab))
    
    div_metric_out[[m]] <- data.frame(matrix(NA, nrow = length(prepped_abun), ncol = 3))
    colnames(div_metric_out[[m]]) <- c("sample", "func", "value")
    div_metric_out[[m]]$sample <- sample_set
    div_metric_out[[m]]$func <- func_set
    
    div_metric_out[[m]]$value <- sapply(prepped_taxa_present, calc_diversity[[m]], tree = in_tree) 
    
    div_metric_out[[m]] <- reshape2::dcast(data = div_metric_out[[m]],
                                           formula = func ~ sample)
    rownames(div_metric_out[[m]]) <- div_metric_out[[m]]$func
    div_metric_out[[m]] <- div_metric_out[[m]][, -1]
    
    metrics <- metrics[-which(metrics == "faiths_pd")]
  }
  
  for (m in metrics) {
    div_metric_out[[m]] <- data.frame(matrix(NA, nrow = length(prepped_abun), ncol = 3))
    colnames(div_metric_out[[m]]) <- c("sample", "func", "value")
    div_metric_out[[m]]$sample <- sample_set
    div_metric_out[[m]]$func <- func_set
    
    div_metric_out[[m]]$value <- sapply(prepped_abun, calc_diversity[[m]])
    
    div_metric_out[[m]] <- reshape2::dcast(data = div_metric_out[[m]],
                                           formula = func ~ sample)
    rownames(div_metric_out[[m]]) <- div_metric_out[[m]]$func
    div_metric_out[[m]] <- div_metric_out[[m]][, -1]
  }
  
  return(div_metric_out)
}


#' @export
beta_div_contrib <- function(func_tab,
                             abun_tab,
                             metrics = c("binary", "kullback"),
                             return_objects = FALSE,
                             write_outfiles = FALSE,
                             outdir = NULL,
                             parDist_func = NULL) {
  
  if (! return_objects && ! write_outfiles) {
    stop("Stopping - either the return_objects or write_outfiles arguments must be set to TRUE, otherwise this function wont actually return or write out anything.") 
  }

  parDist_methods <- c("bhjattacharyya", "bray", "canberra", "chord", 
                       "divergence", "dtw", "euclidean", "fJaccard", "geodesic", 
                       "hellinger", "kullback", "mahalanobis", "manhattan", 
                       "maximum", "minkowski", "podani", "soergel", "wave", 
                       "whittaker", "binary", "braun-blanquet", "dice", "fager", 
                       "faith", "hamman", "kulczynski1", "kulczynski2", "michael", 
                       "mountford", "mozley", "ochiai", "phi", "russel", "simple matching", 
                       "simpson", "stiles", "tanimoto", "yule", "yule2", "cosine", 
                       "hamming", "custom")
  
  if (length(which(! metrics %in% parDist_methods)) > 0) {
    stop("Stopping - at least one specified distance/divergence metric is not available through parDist.") 
  }
    
  if ("custom" %in% metrics && is.null(parDist_func)) {
    stop("Stopping - the parDist_func option must be specified to use a custom metric.") 
  }
  
  almeida_ko <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_annot/kegg_summary.tsv.gz",
                           sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", header = TRUE, check.names = FALSE, row.names = 1)
  
  almeida_abun <- read.table(file = "/data1/gdouglas/projects/contrib_div/data/Almeida_2019/MAG_abun/bwa_depth_min25coverage.tsv.gz",
                             header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")
  
  abun_tab <- almeida_abun
  func_tab <- almeida_ko

  intersecting_features <- colnames(func_tab)[which(colnames(func_tab) %in% rownames(abun_tab))]
  
  if (length(intersecting_features) == 0) {
    stop("Stopping - no features intersect between the two input tables (i.e., between the column names of the function table and row names of the abundance table).")
  }
  
  func_tab <- func_tab[, intersecting_features, drop = FALSE]
  abun_tab <- abun_tab[intersecting_features, , drop = FALSE]
  
  abun_tab <- abun_tab[which(rowSums(abun_tab) > 0), , drop = FALSE]
  
  func_tab <- func_tab[, rownames(abun_tab), drop = FALSE]
  func_tab <- func_tab[which(rowSums(func_tab) > 0), , drop = FALSE]
  func_tab <- func_tab[, which(colSums(func_tab) > 0), drop = FALSE]
  
  if (ncol(func_tab) == 0) {
   stop("Stopping - no features remaining after filtering out those that are all 0's.") 
  }
  
  if (nrow(func_tab) == 0) {
    stop("Stopping - no functions remaining after filtering out those that are all 0's.") 
  }
  
  abun_tab <- abun_tab[colnames(func_tab), , drop = FALSE]
  abun_tab <- abun_tab[, which(colSums(abun_tab) > 0), drop = FALSE]
  
  if (ncol(abun_tab) == 0) {
    stop("Stopping - no samples remaining after filtering out those that are all 0's.") 
  }
  
  # Note that "PATH" will only be defined if this script itself was sourced.
  # Source this file each time the function is called otherwise can get odd behaviour (sometimes).
  # But you also get odd behaviour if this is run separately, so perhaps should be commented out.
  PATH <- "/home/gdouglas/scripts/contrib_div_rough/R_package/"
  Rcpp_filepath <- paste(PATH, "prep_input_list.cpp", sep = "/")
  Rcpp::sourceCpp(Rcpp_filepath)
  
  prepped_contrib <- prep_func_contributor_dimnames(abun_tab = as.matrix(abun_tab),
                                                    func_tab = as.matrix(func_tab))

  
  
  tmp <- abun_tab[which(func_tab[3, ] > 0), ]
  

  
  for (func_i in 1:length(prepped_contrib)) {
  
    func_contrib_abun <- collapse::ss(x = abun_tab,
                                      i = prepped_contrib[[func_i]][[1]] + 1,
                                      j = prepped_contrib[[func_i]][[2]] + 1)
    
    func_contrib_abun <- t(sweep(x = func_contrib_abun,
                                 MARGIN = 2,
                                 STATS = colSums(func_contrib_abun),
                                 FUN = '/'))
    
    func_contrib_abun_bray <- parDist(func_contrib_abun, method = "bray")
    
    func_contrib_abun_bray <- as.matrix(func_contrib_abun_bray)
    
    diag(func_contrib_abun_bray) <- NA
    func_contrib_abun_bray[lower.tri(func_contrib_abun_bray)] <- NA

    data.table::fwrite(x = func_contrib_abun_bray,
                       file = "/home/gdouglas/tmp/test",
                       sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
  }
  
  
  func_contrib_subsets <- parallel::mclapply(1:length(prepped_contrib),
                                       function(func_i) {
                                         as.matrix(collapse::ss(x = abun_tab,
                                            i = prepped_contrib[[func_i]][[1]] + 1,
                                            j = prepped_contrib[[func_i]][[2]] + 1))
                                       },
                                       mc.cores = 10)
  
  
  bray_dists <- parallel::mclapply(1:length(func_contrib_subsets),
                         function(func_i) {
                           parallelDist::parDist(func_contrib_subsets[[func_i]],
                                   method = "bray")
                         },
                         mc.cores = 10)
  
}


#' @export
run_div_metric_strat_long <- function(div_metric, strat_tab, in_tree = NULL) {
  
  unique_funcs <- unique(strat_tab$func)
  unique_samples <- unique(strat_tab$sample)
  
  out <- data.frame(matrix(NA, nrow = length(unique_funcs), ncol = length(unique_samples)))
  rownames(out) <- unique_funcs
  colnames(out) <- unique_samples
  
  for (func in unique_funcs) {
    
    strat_tab_func <- strat_tab[which(strat_tab$func == func), ]
    
    for (s in unique_samples) {
      
      strat_tab_func_sample <- strat_tab_func[which(strat_tab_func$sample == s), ]
      
      if (nrow(strat_tab_func_sample) == 0) { next }
      
      if (div_metric == "faiths_pd") {
        
        out[func, s] <- calc_diversity[[div_metric]](tips_in_sample = strat_tab_func_sample$taxon, tree = in_tree)
        
      } else {
  
        out[func, s] <- calc_diversity[[div_metric]](strat_tab_func_sample$relabun)
        
      }

    }

  }
  
  return(out)

}
