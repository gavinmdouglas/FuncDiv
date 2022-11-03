# Convenience object for programmatically computing beta diversity metrics.
compute_betadiv <- list()

parDist_methods <- c("bhjattacharyya", "bray", "canberra", "chord", 
                     "divergence", "dtw", "euclidean", "fJaccard", "geodesic", 
                     "hellinger", "kullback", "mahalanobis", "manhattan", 
                     "maximum", "minkowski", "podani", "soergel", "wave", 
                     "whittaker", "binary", "braun-blanquet", "dice", "fager", 
                     "faith", "hamman", "kulczynski1", "kulczynski2", "michael", 
                     "mountford", "mozley", "ochiai", "phi", "russel", "simple matching", 
                     "simpson", "stiles", "tanimoto", "yule", "yule2", "cosine", "hamming")


for (m in parDist_methods) {
  func_dist_cmd_char <- paste("function(in_tab, ...) {\n",
                              "func_dist <- as.matrix(parallelDist::parDist(in_tab, method = \"", eval(m), "\"))\n",
                              "func_dist[lower.tri(func_dist, diag = TRUE)] <- NA\n",
                              "return(func_dist)\n}", sep = "")
  
  compute_betadiv[[m]] <- eval(parse(text = func_dist_cmd_char))
    
}

compute_betadiv[["weighted_unifrac"]] <- function(in_tab, in_phylo) {
  fast_weighted_UniFrac(in_tab, in_phylo)
}

compute_betadiv[["unweighted_unifrac"]] <- function(in_tab, in_phylo) {
  fast_unweighted_UniFrac(in_tab, in_phylo)
}

compute_betadiv[["jensen_shannon_div"]] <- function(in_tab, ...) {
  func_dist <- as.matrix(parallelDist::parDist(in_tab, method = "custom", func = jensen_shannon_divergence_FuncPtr))
  func_dist[lower.tri(func_dist, diag = TRUE)] <- NA
  return(func_dist)
}

#' Main function for computing contributional **beta** diversity
#' 
#' Based on joint taxa-function input data (i.e., contributional data), the beta diversity (i.e., inter-sample distance or divergence)
#' will be computed for the subset of taxa encoding each individual function separately. A large List object containing all these tables 
#' can be returned, or alternatively these tables will be written to the disk as plain-text files.
#' 
#' Input data can be either a separate function copy number and taxonomic abundance table, or a joint contributional table.
#' Metrics must be one of "weighted_unifrac", "unweighted_unifrac", "jensen_shannon_div", or a default metric available through the `parallelDist::parDist` function. See `?parallelDist::parDist` for a description of all default metrics.
#' 
#' The taxonomic abundances will be converted to relative abundances prior to computing inter-sample distances.
#' 
#' @param metrics beta diversity metrics to compute. Must be default metric computed by `parallelDist::parDist` or one of \"weighted_unifrac\", \"unweighted_unifrac\", or \"jensen_shannon_div\".
#' @param func_tab data.frame object containing function copy numbers, with rows as functions and columns as taxa. Required if `abun_tab` is specified, and is mutually exclusive with `contrib_tab`.
#' @param abun_tab data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples. Required if `func_tab` is specified, and is mutually exclusive with `contrib_tab`.
#' @param contrib_tab data.frame object containing combined taxa abundances and function copy numbers across taxa. Must contain columns corresponding to the sample ids, function ids, taxa ids, and taxa 
#' abundances within samples. These column names are specified by the `samp_colname`, `func_colname`, `taxon_colname`, and `abun_colname`, respectively.Mutually exclusive with `abun_tab` and `func_tab`. 
#' @param in_tree Phylo object to use if `weighted_unifrac` or `unweighted_unifrac` are specified.
#' @param func_ids character vector specifying subset of function ids to include for analysis. Will analyze all functions present if this is not specified.
#' @param return_objects Boolean vector of length one, specifying whether function should return a list of all output distance tables (nested by metric name, and then by function id). Incompatible with `write_outfiles`.
#' @param write_outfiles Boolean vector of length one, specifying whether function write all distance tables to plain-text files in the specified `outdir` location. Incompatible with `return_objects`.
#' @param outdir character vector of length one, indicating where to save output files if `write_outfiles = TRUE`.
#' @param ncores integer indicating number of cores to use for parallelizable steps.
#' @param samp_colname sample id column name of `contrib_tab` input data.frame.
#' @param func_colname function id column name of `contrib_tab` input data.frame.
#' @param taxon_colname taxon id column name of `contrib_tab` input data.frame.
#' @param abun_colname taxonomic abundance (within each sample) column name of `contrib_tab` input data.frame.
#'
#' @export
beta_div_contrib <- function(metrics = NULL,
                             func_tab = NULL,
                             abun_tab = NULL,
                             contrib_tab = NULL,
                             in_tree = NULL,
                             func_ids = NULL,
                             return_objects = FALSE,
                             write_outfiles = FALSE,
                             outdir = NULL,
                             ncores = 1,
                             samp_colname = "sample",
                             func_colname = "function.",
                             taxon_colname = "taxon",
                             abun_colname = "taxon_abun") {
  
  if (length(metrics) == 0) {
    stop("Stopping - at least one metric must be specified.")  
  }

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
    
    num_func <- nrow(func_tab)
    
  } else if (! is.null(contrib_tab)) {
    
    workflow_type <- "contrib_tab"
    
    if (class(contrib_tab) != "data.frame") {
      stop("Stopping - \"contrib_tab\" input argument must be a data.frame.")
    }
    
    if (length(which(c(samp_colname, func_colname, taxon_colname, abun_colname) %in% colnames(contrib_tab))) < 4) {
      stop("Stopping - at least one of the specified \"samp_colname\", \"func_colname\", \"taxon_colname\", or \"abun_colname\" options is not present as a column name in \"contrib_tab\".")
    }
    
    num_func <- length(unique(contrib_tab[, func_colname]))
    
  }
  
  if (class(return_objects) != "logical" || length(return_objects) != 1) {
    stop("Stopping - \"return_objects\" input argument must be a logical vector of length one.")
  }
  
  if (class(write_outfiles) != "logical" || length(write_outfiles) != 1) {
    stop("Stopping - \"write_outfiles\" input argument must be a logical vector of length one.")
  }
  
  if (class(return_objects) != "logical" || length(return_objects) != 1) {
    stop("Stopping - \"return_objects\" input argument must be a logical vector of length one.")
  }
  
  if ((class(ncores) != "numeric" && class(ncores) != "integer") && length(ncores) != 1) {
    stop("Stopping - \"ncores\" input argument must be an integer vector of length one.")
  } else {
   ncores <- as.integer(ncores) 
  }
  
  if (sum(c(return_objects, write_outfiles)) != 1) {
    stop("Stopping - one (but not both) of the return_objects or write_outfiles arguments must be set to TRUE.") 
  }

  if (write_outfiles) {
    if (is.null(outdir)) {
      stop("Stopping - \"outdir\" argument needs to be specified when \"write_outfiles\" option specified.")
    } else if (class(outdir) != "character" || length(outdir) != 1) {
      stop("Stopping - \"outdir\" argument must be a character vector of length one (if specified).")
    } else if (dir.exists(outdir)) {
      stop("Stopping - specified output directory already exists and will not be overwritten. Please either delete this folder or specify a different folder name to be created.")
    } else if (file.exists(outdir)) {
      stop("Stopping - there is a file with the name of the output directory that you have specified. Please use a different output folder name to be created.")
    } else {
      dir.create(outdir)
    }
  } else if (! is.null(outdir)) {
    stop("Stopping - \"outdir\" argument cannot be set unless \"write_outfiles = TRUE\".") 
  }

  if (is.null(func_ids) && num_func > 100) {
    
    message("The \"func_ids\" argument was not set, so beta diversity will be computed based on the taxonomic contributors for *every* function (i.e., row) of the input function table.")
    message("Note that if there are many functions that this can result in long running times and either extremely large objects returned or many files written.")
  
  } else if (! is.null(func_ids)) {
    
    if (class(func_ids) != "character") {
      stop("Stopping - \"func_ids\" input argument must either be NULL or a character vector.")
    }
    
    if (workflow_type == "multi_tab") {
      if (length(which(! func_ids %in% rownames(func_tab))) > 0) {
        stop("Stopping - the following function ids specified in \"func_ids\" are not present as rows in the function table:\n.  ",
             paste(func_ids[which(! func_ids %in% rownames(func_tab))], collapse = " "))
      }
    
      func_tab <- func_tab[func_ids, ]

    } else if (workflow_type == "contrib_tab") {
      if (length(which(! func_ids %in% contrib_tab[, func_colname])) > 0) {
        stop("Stopping - the following function ids specified in \"func_ids\" are not present as contributional table function id column:\n.  ",
             paste(func_ids[which(! func_ids %in% contrib_tab[, func_colname])], collapse = " "))
      }
      
      contrib_tab <- contrib_tab[which(contrib_tab[, func_colname] %in% func_ids), ]
      
    }
  
  }
  
  other_methods <- c("weighted_unifrac", "unweighted_unifrac", "jensen_shannon_div")
  
  all_methods <- c(parDist_methods, other_methods)
  
  if (length(which(! metrics %in% all_methods) > 0)) {
    stop("Stopping - the following specified distance/divergence metrics are not amongst the possible options for this function (including the metrics available through parallelDist::parDist):\n   ",
         paste(metrics[which(! metrics %in% all_methods)], collapse = " "))
  }
  
  if (any(c("weighted_unifrac", "unweighted_unifrac") %in% metrics)) {

    # Run sanity checks on input tree if UniFrac approach specified.
    if (is.null(in_tree)) {
      stop("Stopping - \"in_tree\" parameter must be set when \"weighted_unifrac\" or \"unweighted_unifrac\" are specified.") 
    }
    
    if (is.null(in_tree$edge.length)) {
      stop("Stopping - input tree has no branch lengths, so UniFrac distances cannot be computed.")
    }
    
    if (! ape::is.rooted(in_tree)) {
      stop("Stopping - rooted input tree required for UniFrac calculation.")
    }
    
  }
  
  if (workflow_type == "multi_tab") {
    subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab,
                                                    abun_table = abun_tab,
                                                    func_ids = func_ids)
    func_tab <- subsetted_tables$func
    abun_tab <- subsetted_tables$abun
    rm(subsetted_tables)
    
    # Sort function table in ascending order by function, and abundance table by sample id.
    func_tab <- func_tab[sort(rownames(func_tab), decreasing = FALSE), ]
    abun_tab <- abun_tab[, sort(colnames(abun_tab), decreasing = FALSE)]
    
    prepped_contrib <- prep_func_contributor_dimnames(abun_tab = as.matrix(abun_tab),
                                                      func_tab = as.matrix(func_tab))
    
    func_contrib_subsets <- parallel::mclapply(1:length(prepped_contrib),
                                         function(func_i) {
                                           func_contrib_abun <- as.matrix(collapse::ss(x = abun_tab,
                                                                                       i = prepped_contrib[[func_i]][[1]] + 1,
                                                                                       j = prepped_contrib[[func_i]][[2]] + 1))
                                           
                                           func_contrib_abun <- t(sweep(x = func_contrib_abun,
                                                                        MARGIN = 2,
                                                                        STATS = colSums(func_contrib_abun),
                                                                        FUN = '/'))
                                           
                                           # Sort column names to be alphabetical order, which is needed to make sure the
                                           # dtw metric produces the same result with either the multi-table or contrib table inputs.
                                           func_contrib_abun[, sort(colnames(func_contrib_abun), decreasing = FALSE), drop = FALSE]
                                           
                                         },
                                         mc.cores = ncores)
    
    func_ordered <- rownames(func_tab)

  } else if (workflow_type == "contrib_tab") {
    
    contrib_tab <- contrib_tab[, c(samp_colname, func_colname, taxon_colname, abun_colname)]
    
    contrib_tab_split <- collapse::rsplit(x = contrib_tab, by = contrib_tab[, func_colname])
    
    func_contrib_subsets <- parallel::mclapply(contrib_tab_split,
                                        function(x) {
                                          tab <- data.table::dcast.data.table(data = data.table::data.table(x),
                                                                      formula = x[, taxon_colname] ~ x[, samp_colname],
                                                                      value.var = abun_colname)
                                          
                                          tab[is.na(tab)] <- 0
                                          
                                          rownames(tab) <- tab$x
                                          
                                          tab <- tab[, -1, drop = FALSE]
                                          
                                          tab <- t(sweep(x = tab, MARGIN = 2, STATS = colSums(tab), '/'))
                                          
                                          return(tab)
                                        },
                                        mc.cores = ncores)

    func_ordered <- names(contrib_tab_split)
    
  }
  
  if (return_objects) {

    all_dists <- list()
    
    for (metric in metrics) {

      all_dists[[metric]] <- parallel::mclapply(1:length(func_contrib_subsets),
    
                                     function(func_i) {
            
                                         func_dist <- compute_betadiv[[metric]](func_contrib_subsets[[func_i]], in_tree)
                                         
                                         return(data.frame(func_dist, check.names = FALSE))
                                       
                                       }, mc.cores = ncores)

      names(all_dists[[metric]]) <- func_ordered

    }
    
    return(all_dists)
  
  } else if (write_outfiles) {
 
    for (metric in metrics) {

      dir.create(paste(outdir, metric, sep = "/"))
      
      null_return <- parallel::mclapply(1:length(func_contrib_subsets),
                                          
                                                  function(func_i) {
                                                    
                                                      func_id <- rownames(func_tab)[func_i]
                                                      
                                                      outfile <- paste(outdir, "/", metric, "/", func_id, ".tsv", sep = "")
                                                      
                                                      func_dist <- compute_betadiv[[metric]](func_contrib_subsets[[func_i]], in_tree)

                                                      data.table::fwrite(x = data.table::data.table(func_dist),
                                                                         file = outfile,
                                                                         sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
                                                      
                                                      return(NULL)
                                                  },
                                                  mc.cores = 10)
    }
    
    return(paste("Results written to:", outdir, sep = " "))
  
  }

}
