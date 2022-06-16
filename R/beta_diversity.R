
#' @export
beta_div_contrib <- function(func_tab,
                             abun_tab,
                             metrics = c("binary", "kullback"),
                             func_ids = NULL,
                             return_objects = FALSE,
                             write_outfiles = FALSE,
                             outdir = NULL,
                             ncores = 1,
                             parDist_func = NULL) {
  
  if (sum(c(return_objects, write_outfiles) != 1)) {
    stop("Stopping - one (but not both) of the return_objects or write_outfiles arguments must be set to TRUE.") 
  }

  if (write_outfiles) {
    if (is.null(outdir)) {
      stop("Stopping - outdir needs to be specified when write_outfiles option specified.")
    } else if (dir.exists(outdir)) {
      stop("Stopping - outdir already exists and will not be overwritten. Please either delete this folder or specify a different folder name to be created.")
    } else if (file.exists(outdir)) {
      stop("Stopping - there is a file with the name of the output directory that you have specified. Please use a different output folder name to be created.")
    } else {
      dir.create(outdir)
    }
  } else if (! is.null(outdir)) {
    stop("Stopping - outdir argument cannot be set unless write_outfiles = TRUE.") 
  }

  if (is.null(func_ids)) {
    message("The func_ids argument was not set, so beta diversity will be computed based on the taxonomic contributors for *every* function (i.e., row) of the input function table.")
    message("Note that if there are many functions that this can result in long running times and either extremely large objects returned or many files written.")
  } else {
    
    if (length(which(! func_ids %in% rownames(func_tab))) > 0) {
      stop("Stopping - the following function ids specified in func_ids are not present as rows in the function table:\n",
           paste(func_ids[which(! func_ids %in% rownames(func_tab))], collapse = " "))
    }
    
    func_tab <- func_tab[func_ids, ]
  
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
   stop("Stopping - no features remaining after restricting to taxa in the taxa abundance table, and filtering out those that are all 0's.") 
  }
  
  if (is.null(func_ids) && nrow(func_tab) == 0) {
    stop("Stopping - no functions remaining after restricting to taxa in the taxa abundance table, and filtering out those that are all 0's.") 
  } else if(! is.null(func_ids) && nrow(func_tab) < length(func_ids)) {
    stop("Stopping - the above function ids specified in func_ids are not present as rows in the function table after restricting to taxa in the taxa abundance table, and filtering out rows that are all 0's:\n",
         paste(func_ids[which(! func_ids %in% rownames(func_tab))], collapse = " "))
  }
  
  abun_tab <- abun_tab[colnames(func_tab), , drop = FALSE]
  abun_tab <- abun_tab[, which(colSums(abun_tab) > 0), drop = FALSE]
  
  if (ncol(abun_tab) == 0) {
    stop("Stopping - no samples remaining after filtering out those that are all 0's.") 
  }
  
  prepped_contrib <- prep_func_contributor_dimnames(abun_tab = as.matrix(abun_tab),
                                                    func_tab = as.matrix(func_tab))
  
  func_contrib_subsets <- parallel::mclapply(1:length(prepped_contrib),
                                       function(func_i) {
                                         func_contrib_abun <- as.matrix(collapse::ss(x = abun_tab,
                                                                                     i = prepped_contrib[[func_i]][[1]] + 1,
                                                                                     j = prepped_contrib[[func_i]][[2]] + 1))
                                         
                                         t(sweep(x = func_contrib_abun,
                                                 MARGIN = 2,
                                                 STATS = colSums(func_contrib_abun),
                                                 FUN = '/'))
                                         
                                       },
                                       mc.cores = ncores)
  
  if (return_objects) {

    all_dists <- list()
    
    for (metric in metrics) {
    
      all_dists[[metric]] <- parallel::mclapply(1:length(func_contrib_subsets),
    
                                     function(func_i) {
            
                                         func_dist <- as.matrix(parallelDist::parDist(func_contrib_subsets[[func_i]],
                                                                                      method = metric,
                                                                                      func = parDist_func))
                                         
                                         diag(func_dist) <- NA
                                         
                                         func_dist[lower.tri(func_dist)] <- NA
                                         
                                         return(func_dist)
                                       
                                       }, mc.cores = ncores)

      names(all_dists[[metric]]) <- rownames(func_tab)
      
    }
    
    return(all_dists)
  
  } else if (write_outfiles) {
 
    for (metric in metrics) {

      dir.create(paste(outdir, metric, sep = "/"))
      
      null_return <- parallel::mclapply(1:length(func_contrib_subsets),
                                        
                                                function(func_i) {
                                                  
                                                    func_id <- rownames(func_tab)[func_i]
                                                    
                                                    outfile <- paste(outdir, "/", metric, "/", func_id, ".tsv", sep = "")
                                                    
                                                    func_dist <- as.matrix(parallelDist::parDist(func_contrib_subsets[[func_i]],
                                                                                                 method = metric,
                                                                                                 func = parDist_func))
                                                    diag(func_dist) <- NA
                                                    
                                                    func_dist[lower.tri(func_dist)] <- NA
                                                    
                                                    data.table::fwrite(x = data.table::data.table(func_dist),
                                                                       file = outfile,
                                                                       sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
                                                    
                                                    return(NULL)
                                                },
                                                mc.cores = 10)
    }
    
    return("Results written to specified output directory.")
  
  }

}
