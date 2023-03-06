#' Main function for computing contributional **alpha** diversity
#' 
#' Based on joint taxa-function input data (i.e., contributional data), a dataframe
#' will be returned for each specified metric, which will contain the metric values for all function and sample combinations.
#' 
#' Input data can be either a separate function copy number and taxonomic abundance table, or a joint contributional table.
#' By default, specified metrics must be one of `names(FuncDiv_alpha_metrics)`. However, custom alpha diversity metric functions
#' can be specified with the `custom_metric_functions` parameter.
#' 
#' Note that the taxonomic abundances can be relative abundance, read counts, or transformed in another way. However, note that some default metrics
#' are only compatible with count data (see `?FuncDiv_alpha_metrics`).
#' 
#' @param metrics alpha diversity metrics to compute. Must either be names of functions in `FuncDiv_alpha_metrics`, or alternatively in `custom_metric_functions`, if specified. 
#' @param func_tab data.frame object containing function copy numbers, with rows as functions and columns as taxa. Required if `abun_tab` is specified, and is mutually exclusive with `contrib_tab`.
#' @param abun_tab data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples. Required if `func_tab` is specified, and is mutually exclusive with `contrib_tab`.
#' @param contrib_tab data.frame object containing combined taxa abundances and function copy numbers across taxa. Must contain columns corresponding to the sample ids, function ids, taxa ids, and taxa 
#' abundances within samples. These column names are specified by the `samp_colname`, `func_colname`, `taxon_colname`, and `abun_colname`, respectively. Mutually exclusive with `abun_tab` and `func_tab`. 
#' @param in_tree phylo object to use if `faiths_pd` is specified.
#' @param ncores integer indicating number of cores to use for parallelizable steps.
#' @param replace_NA Boolean vector of length one, indicating whether all NA's in the output of all metrics should be converted to 0's. Note that this done automatically done for `richness` either way.
#' @param custom_metric_functions List object containing custom alpha diversity metric functions. This overrides `FuncDiv_alpha_metrics` when specified. The list element names must correspond to at
#' least the names indicated by the `metrics` parameter.
#' @param samp_colname sample id column name of `contrib_tab` input data.frame.
#' @param func_colname function id column name of `contrib_tab` input data.frame.
#' @param taxon_colname taxon id column name of `contrib_tab` input data.frame.
#' @param abun_colname taxonomic abundance (within each sample) column name of `contrib_tab` input data.frame.
#'
#' @return a list, containing one dataframe for each specified alpha diversity metric.
#' In each dataframe, rows are functions and samples are columns.
#'
#' @examples
#' # First, simulate some (non-realistic) data.
#' set.seed(123)
#' test_tree <- ape::rtree(100)
#' test_abun <- data.frame(matrix(rnorm(500), nrow = 100, ncol = 5))
#' rownames(test_abun) <- test_tree$tip.label
#' colnames(test_abun) <- c("sample1", "sample2", "sample3", "sample4", "sample5")
#' test_abun[test_abun < 0] <- 0
#' test_func <- data.frame(matrix(sample(c(0L, 1L), 200, replace = TRUE),
#'                                nrow = 2, ncol = 100))
#' colnames(test_func) <- test_tree$tip.label
#' rownames(test_func) <- c("func1", "func2")
#'
#' # Compute alpha diversity, based on (observed) richness, Faith's phylogenetic
#' # diversity, and the Gini-Simpson Index.
#' contrib_alpha <- alpha_div_contrib(metrics = c("richness",  "faiths_pd", "gini_simpson_index"),
#'                                    func_tab = test_func,
#'                                    abun_tab = test_abun,
#'                                    in_tree = test_tree,
#'                                    ncores = 1)
#'
#' # Print out computed Gini-Simpson Index values.
#' contrib_alpha$gini_simpson_index
#'
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
                              taxon_colname = "taxon",
                              abun_colname = "taxon_abun") {
  
  if (! inherits(metrics, "character")) {
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

    if (! inherits(func_tab, "data.frame")) {
      stop("Stopping - \"func_tab\" input argument must be a data.frame.")
    }

    if (! inherits(abun_tab, "data.frame")) {
      stop("Stopping - \"abun_tab\" input argument must be a data.frame.")
    }

  } else if (! is.null(contrib_tab)) {

    workflow_type <- "contrib_tab"

    if (! inherits(contrib_tab, "data.frame")) {
      stop("Stopping - \"contrib_tab\" input argument must be a data.frame.")
    }
    
    if (length(which(c(samp_colname, func_colname, taxon_colname, abun_colname) %in% colnames(contrib_tab))) < 4) {
      stop("Stopping - at least one of the specified \"samp_colname\", \"func_colname\", \"taxon_colname\", or \"abun_colname\" options is not present as a column name in \"contrib_tab\".")
    }

  }

  if (! inherits(replace_NA, "logical") || length(replace_NA) != 1) {
    stop("Stopping - \"replace_NA\" input argument must be a logical vector of length one.")
  }

  if (is.null(custom_metric_functions)) {
    
    metric_functions <- FuncDiv_alpha_metrics
    
    if (length(which(! metrics %in% names(FuncDiv_alpha_metrics))) > 0) {
      stop("Stopping - the following specified metrics are not names in the \"FuncDiv_alpha_metrics\" object:\n   ",
           paste(metrics[which(! metrics %in% names(FuncDiv_alpha_metrics))], collapse = " "), "\n\n",
           "  You can see all available alpha diversity metrics by typing \"names(FuncDiv_alpha_metrics)\".")
    }

  } else {
    
    metric_functions <- custom_metric_functions
    
    if (! inherits(custom_metric_functions, "list")) {
      stop("Stopping - \"custom_metric_functions\" argument must be a list (if specified).")
    } else if (length(which(! metrics %in% names(custom_metric_functions))) > 0) {
      stop("Stopping - the following specified metrics are not names in the list of custom metric functions that was input:\n   ",
           paste(metrics[which(! metrics %in% names(custom_metric_functions))], collapse = " "))
    }
  }
  
  if ("faiths_pd" %in% metrics) {
    if (is.null(in_tree)) {
      stop("Stopping - a phylo (i.e., a tree) object must be passed to the \"in_tree\" argument to calculate faiths_pd.")
    } else if (! inherits(in_tree, "phylo")) {
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
   
    contrib_tab <- contrib_tab[which(contrib_tab[, abun_colname] > 0), , drop = FALSE]
    
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
                    
    div_metric_out[["faiths_pd"]] <- data.frame(data.table::dcast.data.table(data = data.table::data.table(div_metric_out[["faiths_pd"]]),
                                                                             formula = func ~ sample), check.names = FALSE)
    
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
      
      div_metric_out[[m]] <- data.frame(data.table::dcast.data.table(data = data.table::data.table(div_metric_out[[m]]),
                                                                     formula = func ~ sample), check.names = FALSE)
        
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
