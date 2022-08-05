
#' @export
alpha_div_contrib <- function(metrics,
                              func_tab = NULL,
                              abun_tab = NULL,
                              contrib_tab = NULL,
                              in_tree = NULL,
                              ncores = 1,
                              replace_NA = FALSE,
                              custom_metric_functions = NULL,
                              samp_colname = "sample",
                              func_colname = "function.",
                              abun_colname = "taxon_abun",
                              taxon_colname = "taxon") {
  
  if (class(metrics) != "character") {
    stop("Stopping - \"metrics\" input argument must be a character vector.")
  }
  
  if (is.null(func_tab) && is.null(abun_tab) && is.null(contrib_tab)) {
    stop("Stopping - either \"func_tab\" and \"abun_tab\" *or* \"contrib_tab\" must be specified.")
  }
  
  if (! is.null(func_tab) && ! is.null(abun_tab) && ! is.null(contrib_tab)) {
    stop("Stopping - either \"func_tab\" and \"abun_tab\" *or* \"contrib_tab\" must be specified, not all!")
  }
  
  if ((is.null(func_tab) && ! is.null(abun_tab)) || (! is.null(func_tab) && is.null(abun_tab))) {
    stop("Stopping - either both \"func_tab\" and \"abun_tab\" or neither must be specified")
  }
  
  if ((is.null(func_tab) && ! is.null(abun_tab)) || (! is.null(func_tab) && is.null(abun_tab))) {
    stop("Stopping - either both \"func_tab\" and \"abun_tab\" or neither must be specified")
  }

  if (! is.null(func_tab) && ! is.null(abun_tab)) {
    
    workflow_type <- "multi_tab"

    if (class(func_tab) != "data.frame") {
      stop("Stopping - \"func_tab\" input argument must be a data.frame.")
    }

    if (class(abun_tab) != "data.frame") {
      stop("Stopping - \"abun_tab\" input argument must be a data.frame.")
    }

  } else if (! is.null(contrib_tab)) {

    workflow_type <- "contrib_tab"

    if (class(contrib_tab) != "data.frame") {
      stop("Stopping - \"contrib_tab\" input argument must be a data.frame.")
    }
    
    if (length(which(c(samp_colname, func_colname, taxon_colname, abun_colname) %in% colnames(contrib_tab))) < 4) {
      stop("Stopping - at least one of the specified \"samp_colname\", \"func_colname\", \"taxon_colname\", or \"abun_colname\" options is not present as a column name in \"contrib_tab\".")
    }

  }

  if (class(replace_NA) != "logical" || length(replace_NA) != 1) {
    stop("Stopping - \"replace_NA\" input argument must be a logical vector of length one.")
  }

  if (is.null(custom_metric_functions)) {
    
    metric_functions <- FuncDiv_abun_alpha_metrics
    
    if (length(which(! metrics %in% names(FuncDiv_abun_alpha_metrics))) > 0) {
      stop("Stopping - the following specified metrics are not names in the \"FuncDiv_abun_alpha_metrics\" object:\n   ",
           paste(metrics[which(! metrics %in% names(FuncDiv_abun_alpha_metrics))], collapse = " "), "\n\n",
           "  You can see all available alpha diversity metrics by typing \"names(FuncDiv_abun_alpha_metrics)\".")
    }

  } else {
    
    metric_functions <- custom_metric_functions
    
    if (class(custom_metric_functions) != "list") {
      stop("Stopping - \"custom_metric_functions\" argument must be a list (if specified).")
    } else if (length(which(! metrics %in% names(custom_metric_functions))) > 0) {
      stop("Stopping - the following specified metrics are not names in the list of custom metric functions that was input:\n   ",
           paste(metrics[which(! metrics %in% names(custom_metric_functions))], collapse = " "))
    }
  }
  
  if ("faiths_pd" %in% metrics) {
    if (is.null(in_tree)) {
      stop("Stopping - a phylo (i.e., a tree) object must be passed to the \"in_tree\" argument to calculate faiths_pd.")
    } else if (class(in_tree) != "phylo") {
      stop("Stopping - \"in_tree\" is not a phylo object.")
    } else if (is.null(in_tree$edge.length)) {
      stop("Stopping - no branch lengths in tree. faiths_pd will not be possible to compute.")
    } else if (! ape::is.rooted(in_tree)) {
      stop("Stopping - tree is unrooted - please add root and re-run.")
    }
  } else if(! is.null(in_tree)) {
    message("Note that it was not necessary to specify the \"in_tree\" argument, because faiths_pd was not amongst the selected metrics.")
  }

  if (workflow_type == "multi_tab") {
  
    subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab,
                                                    abun_table = abun_tab)
    func_tab <- subsetted_tables$func
    abun_tab <- subsetted_tables$abun
    rm(subsetted_tables)
    
    func_set <- as.character()
    sample_set <- as.character()
    for (i in 1:nrow(func_tab)) {
      func_set <- c(func_set, rep(rownames(func_tab)[i], ncol(abun_tab)))
      sample_set <- c(sample_set, colnames(abun_tab))
    }
    
  } else if (workflow_type == "contrib_tab") {
   
    contrib_tab <- contrib_tab[, c(samp_colname, func_colname, taxon_colname, abun_colname)]
    
    contrib_tab <- contrib_tab[collapse::radixorder(contrib_tab[, samp_colname], contrib_tab[, func_colname]), ]
    
    group_starts <- attr(collapse::radixorder(contrib_tab[, samp_colname], contrib_tab[, func_colname], starts = TRUE, sort = FALSE), "starts")
    num_groups <- length(group_starts)
    
    sample_set <- contrib_tab[group_starts, samp_colname]
    func_set <- contrib_tab[group_starts, func_colname]
    
  }
  
  div_metric_out <- list()
  
  # Get list of taxa ids that are present in each sample / func combination if faiths_pd specified as well.
  if ("faiths_pd" %in% metrics) {
    
    if (workflow_type == "multi_tab") {
  
      prepped_taxa_present <- prep_all_sample_func_taxa_vec(abun_tab = as.matrix(abun_tab),
                                                            func_tab = as.matrix(func_tab))
      
      prepped_taxa_present_nonzero_i <- which(sapply(prepped_taxa_present, length) > 0)
      
      prepped_taxa_present <- prepped_taxa_present[prepped_taxa_present_nonzero_i]
      sample_set_taxa_present <- sample_set[prepped_taxa_present_nonzero_i]
      func_set_taxa_present <- func_set[prepped_taxa_present_nonzero_i]
    
    } else if (workflow_type == "contrib_tab") {
      
      prepped_taxa_present <- parallel::mclapply(1:(num_groups - 1),
                                                function(i) {
                                                  contrib_tab[group_starts[i]:(group_starts[i + 1] - 1), taxon_colname]
                                                },
                                                mc.cores = ncores)
      
      prepped_taxa_present[[num_groups]] <- contrib_tab[group_starts[num_groups]:nrow(contrib_tab), taxon_colname]
      
      sample_set_taxa_present <- sample_set
      func_set_taxa_present <- func_set
  
    }
    
    div_metric_out[["faiths_pd"]] <- data.frame(matrix(NA, nrow = length(prepped_taxa_present), ncol = 3))
    colnames(div_metric_out[["faiths_pd"]]) <- c("sample", "func", "value")
    div_metric_out[["faiths_pd"]]$sample <- sample_set_taxa_present
    div_metric_out[["faiths_pd"]]$func <- func_set_taxa_present
    
    div_metric_out[["faiths_pd"]]$value <- unlist(parallel::mclapply(prepped_taxa_present,
                                                                     metric_functions[["faiths_pd"]],
                                                                     tree = in_tree,
                                                                     mc.cores = ncores))
                    
    div_metric_out[["faiths_pd"]] <- reshape2::dcast(data = div_metric_out[["faiths_pd"]],
                                                     formula = func ~ sample)
    
    rownames(div_metric_out[["faiths_pd"]]) <- div_metric_out[["faiths_pd"]]$func
    
    div_metric_out[["faiths_pd"]] <- div_metric_out[["faiths_pd"]][, -1]
    
    metrics <- metrics[-which(metrics == "faiths_pd")]

  }
  
  if (length(which(metrics != "faiths_pd")) > 0) {
  
    if (workflow_type == "multi_tab") {
      
      prepped_abun <- prep_all_sample_func_vec(abun_tab = as.matrix(abun_tab),
                                               func_tab = as.matrix(func_tab))
      
      prepped_abun_nonzero_i <- which(sapply(prepped_abun, length) > 0)
      
      prepped_abun <- prepped_abun[prepped_abun_nonzero_i]
      sample_set_nonzero_abun <- sample_set[prepped_abun_nonzero_i]
      func_set_nonzero_abun <- func_set[prepped_abun_nonzero_i]
      
    } else if (workflow_type == "contrib_tab") {
      
      prepped_abun <- parallel::mclapply(1:(num_groups - 1),
                                         function(i) {
                                           contrib_tab[group_starts[i]:(group_starts[i + 1] - 1), abun_colname]
                                         },
                                         mc.cores = ncores)
                
      prepped_abun[[num_groups]] <- contrib_tab[group_starts[num_groups]:nrow(contrib_tab), abun_colname]
      
      sample_set_nonzero_abun <- sample_set
      func_set_nonzero_abun <- func_set
    
    }
    
    for (m in metrics) {
  
      div_metric_out[[m]] <- data.frame(matrix(NA, nrow = length(prepped_abun), ncol = 3))
      colnames(div_metric_out[[m]]) <- c("sample", "func", "value")
      div_metric_out[[m]]$sample <- sample_set_nonzero_abun
      div_metric_out[[m]]$func <- func_set_nonzero_abun
      
      div_metric_out[[m]]$value <- unlist(parallel::mclapply(prepped_abun,
                                                             metric_functions[[m]],
                                                             mc.cores = ncores))
      
      div_metric_out[[m]] <- reshape2::dcast(data = div_metric_out[[m]],
                                             formula = func ~ sample)
      rownames(div_metric_out[[m]]) <- div_metric_out[[m]]$func
      div_metric_out[[m]] <- div_metric_out[[m]][, -1]
      
    }
    
    # If richness was one of the metrics specified, then replace NA's with 0's for that table
    if ("richness" %in% names(div_metric_out)) {
      div_metric_out[["richness"]][is.na(div_metric_out[["richness"]])] <- 0
    }
  
    if (replace_NA) {
      for (m in names(div_metric_out)[which(names(div_metric_out) != "richness")]) {
        div_metric_out[[m]][is.na(div_metric_out[[m]])] <- 0
      }
    }
    
  }
  
  return(div_metric_out)
}
