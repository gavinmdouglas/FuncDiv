#' Utility function to convert from contributional to multi-table input objects
#' 
#' Converts from contributional-type table (i.e., a single, long table with joint taxa/function information) to separate taxa abundance and function copy number tables.
#' 
#' @param contrib_tab data.frame object containing combined taxa abundances and function copy numbers across taxa. Must contain columns corresponding to the sample ids, function ids, taxa ids, and taxa 
#' abundances within samples. These column names are specified by the `samp_colname`, `func_colname`, `taxon_colname`, `abun_colname`, and `copy.num_colname`, respectively. 
#' @param samp_colname sample id column name of `contrib_tab` input data.frame.
#' @param func_colname function id column name of `contrib_tab` input data.frame.
#' @param taxon_colname taxon id column name of `contrib_tab` input data.frame.
#' @param abun_colname taxonomic abundance (within each sample) column name of `contrib_tab` input data.frame.
#' @param copy.num_colname function copy number column name of `contrib_tab` input data.frame.
#'
#' @return list with taxon abundance (`taxon_abun`) and function copy number (`function_copy_num`) data.frames as separate elements.
#'
#' @export
contrib_to_multitab <- function(contrib_tab,
                                samp_colname = "sample",
                                func_colname = "function.",
                                abun_colname = "taxon_abun",
                                taxon_colname = "taxon",
                                copy.num_colname = "genome_function_count") {
  
  if (! inherits(contrib_tab, "data.frame")) {
    stop("Stopping - \"contrib_tab\" input argument must be a data.frame.")
  }
  
  if (length(which(c(samp_colname, func_colname, taxon_colname, abun_colname) %in% colnames(contrib_tab))) < 4) {
    stop("Stopping - at least one of the specified \"samp_colname\", \"func_colname\", \"taxon_colname\", \"abun_colname\", \"copy.num_colname\" options is not present as a column name in \"contrib_tab\".")
  }
  
  contrib_tab <- contrib_tab[, c(samp_colname, func_colname, taxon_colname, abun_colname, copy.num_colname)]
  
  
  # Taxon abundance table
  contrib_taxa_abun <- contrib_tab[, c(samp_colname, taxon_colname, abun_colname)]
  
  contrib_taxa_abun <- contrib_taxa_abun[which(! duplicated(contrib_taxa_abun)), ]
  
  duplicated_sample_taxa_combos_i <- which(duplicated(contrib_taxa_abun[, c(samp_colname, taxon_colname)]))
  
  if (length(duplicated_sample_taxa_combos_i) > 0) {
    stop("Stopping - at least one combination of samples and taxa has a different specificied abundance across different functions. This could indicate an issue with the contributional format in the first place, and means that it is impossible to convert to the two table format.")
  }
  
  contrib_taxa_abun_wide <- data.frame(data.table::dcast.data.table(data = data.table::data.table(contrib_taxa_abun),
                                                                    formula = stats::as.formula(paste(taxon_colname, "~", samp_colname, sep = " ")),
                                                                    value.var = abun_colname), check.names = FALSE)
  
  rownames(contrib_taxa_abun_wide) <- as.character(contrib_taxa_abun_wide[, taxon_colname])
  contrib_taxa_abun_wide <- contrib_taxa_abun_wide[, -1]
  contrib_taxa_abun_wide[is.na(contrib_taxa_abun_wide)] <- 0
  
  
  # Function copy number table.
  contrib_tab_func_copy_num <- contrib_tab[, c(func_colname, taxon_colname, copy.num_colname)]
  
  contrib_tab_func_copy_num <- contrib_tab_func_copy_num[which(! duplicated(contrib_tab_func_copy_num)), ]
  
  duplicated_func_taxa_combos_i <- which(duplicated(contrib_tab_func_copy_num[, c(func_colname, taxon_colname)]))
  
  if (length(duplicated_func_taxa_combos_i) > 0) {
    stop("Stopping - at least one combination of taxa / functions has a different specificied copy number across samples. This can happen for instance if pathway levels per taxon are computed based on how much they contribute to the community-wide pathway abundance. This means that it is impossible to convert to the two table format.")
  }
  
  contrib_tab_func_copy_num_wide <- data.frame(data.table::dcast.data.table(data = data.table::data.table(contrib_tab_func_copy_num),
                                                                            formula = stats::as.formula(paste(func_colname, "~", taxon_colname, sep = " ")),
                                                                            value.var = copy.num_colname), check.names = FALSE)
  
  rownames(contrib_tab_func_copy_num_wide) <- as.character(contrib_tab_func_copy_num_wide[, func_colname])
  contrib_tab_func_copy_num_wide <- contrib_tab_func_copy_num_wide[, -1]
  contrib_tab_func_copy_num_wide[is.na(contrib_tab_func_copy_num_wide)] <- 0
  
  return(list(taxon_abun = contrib_taxa_abun_wide,
              function_copy_num = contrib_tab_func_copy_num_wide))
  
}


#' Utility function to convert from multi-table objects to contributional table
#' 
#' Converts from separate taxa abundance and function copy number table input style to contributional-type table (i.e., a single, long table with joint taxa/function information). 
#' @param func_tab data.frame object containing function copy numbers, with rows as functions and columns as taxa.
#' @param abun_tab data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples.
#' @param ncores integer specifying number of cores to use for parallizable steps.
#' @param samp_colname sample id column name of `contrib_tab` output data.frame.
#' @param func_colname function id column name of `contrib_tab` output data.frame.
#' @param taxon_colname taxon id column name of `contrib_tab` output data.frame.
#' @param abun_colname taxonomic abundance (within each sample) column name of `contrib_tab` output data.frame.
#' @param copy.num_colname function copy number (within each taxa) column name of `contrib_tab` output data.frame.
#'
#' @return data.frame in contributional format (i.e., single, long-format version of both input tables).
#'
#' @export
multitab_to_contrib <- function(func_tab,
                                abun_tab,
                                ncores = 1,
                                samp_colname = "sample",
                                func_colname = "function.",
                                abun_colname = "taxon_abun",
                                taxon_colname = "taxon",
                                copy.num_colname = "genome_function_count") {

  if (! inherits(func_tab, "data.frame")) {
    stop("Stopping - \"func_tab\" input argument must be a data.frame.")
  }
  
  if (! inherits(abun_tab, "data.frame")) {
    stop("Stopping - \"abun_tab\" input argument must be a data.frame.")
  }
  
  subsetted_tables <- subset_func_and_abun_tables(func_table = func_tab,
                                                  abun_table = abun_tab)
  func_tab <- subsetted_tables$func
  abun_tab <- subsetted_tables$abun
  rm(subsetted_tables)
  
  contrib_blocks <- parallel::mclapply(X = rownames(func_tab),
                                       FUN = function(func_id) {
                                         
                                         taxa_encoders <- colnames(func_tab)[which(func_tab[func_id, ]  > 0)]
                                         
                                         abun_tab_subset <- abun_tab[taxa_encoders, , drop = FALSE]
                                         abun_tab_subset <- abun_tab_subset[, which(colSums(abun_tab_subset) > 0), drop = FALSE]
                                         
                                         abun_tab_subset[, taxon_colname] <- rownames(abun_tab_subset)
                                         
                                         contrib_block <- data.frame(data.table::melt.data.table(data = data.table::data.table(abun_tab_subset),
                                                                                                 id.vars = taxon_colname,
                                                                                                 variable.name = samp_colname,
                                                                                                 value.name = abun_colname), check.names = FALSE)

                                         contrib_block <- contrib_block[which(contrib_block[, abun_colname] > 0), , drop = FALSE]
                                         
                                         contrib_block[, samp_colname] <- as.character(contrib_block[, samp_colname])
                                         
                                         contrib_block[, func_colname] <- func_id
                                         
                                         contrib_block[, copy.num_colname] <- as.numeric(func_tab[func_id, contrib_block[, taxon_colname]])
                                         
                                         return(contrib_block)
                                         
                                       },
                                       mc.cores = ncores)
  
  contrib_tab <- do.call(rbind, contrib_blocks)
  
  rownames(contrib_tab) <- NULL
  
  return(contrib_tab[, c(samp_colname, func_colname, taxon_colname, copy.num_colname, abun_colname), drop = FALSE])

}


#' Utility function to get community-wide function abundance table
#' 
#' Takes in table of function copy numbers across taxa and table of taxa abundances across samples.
#' I.e., it represents the multiplication of the function copy numbers by the abundances of the taxa within each sample.
#'
#' @param func_tab data.frame object containing function copy numbers, with rows as functions and columns as taxa.
#' @param abun_tab data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples.
#'
#' @return data.frame representing the *unnormalized* community-wide abundances of functions across samples.
#'
#' @export
func_abun_crossproduct <- function(func_tab, abun_tab) {
  
  # Check that all column names in function table.
  if(length(which(! rownames(abun_tab) %in% colnames(func_tab))) > 0) {
    stop("Stoppings - some rows in abundance table not found in function table.")
  }
  
  func_tab <- func_tab[, rownames(abun_tab)]
  
  return(data.frame(t(crossprod(as.matrix(abun_tab), t(func_tab))), check.names = FALSE))
  
}


#' Utility function to subset function copy number and taxonomic abundance tables
#' 
#' The input tables will be returned except subset to the same taxa ids. Any functions and / or samples that are totally absent after this step will be dropped.
#' 
#' @param func_table data.frame object containing function copy numbers, with rows as functions and columns as taxa.
#' @param abun_table data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples.
#' @param func_ids optional character vector of function ids to retain (all other rows of `func_tab` will be removed).
#'
#' @return list containing subsetted function and abundance data.frames as separate elements.
#'
#' @export
subset_func_and_abun_tables <- function(func_table, abun_table, func_ids = NULL) {

  intersecting_features <- colnames(func_table)[which(colnames(func_table) %in% rownames(abun_table))]
  
  if (length(intersecting_features) == 0) {
    stop("Stopping - no features intersect between the two input tables (i.e., between the column names of the function table and row names of the abundance table).")
  }

  func_table <- func_table[, intersecting_features, drop = FALSE]
  abun_table <- abun_table[intersecting_features, , drop = FALSE]
  
  abun_table <- abun_table[which(rowSums(abun_table) > 0), , drop = FALSE]
  
  func_table <- func_table[, rownames(abun_table), drop = FALSE]
  func_table <- func_table[which(rowSums(func_table) > 0), , drop = FALSE]
  func_table <- func_table[, which(colSums(func_table) > 0), drop = FALSE]
  
  if (ncol(func_table) == 0) {
    stop("Stopping - no features remaining after restricting to taxa in the taxa abundance table, and filtering out those that are all 0's.") 
  }
  
  if (is.null(func_ids) && nrow(func_table) == 0) {
    stop("Stopping - no functions remaining after restricting to taxa in the taxa abundance table, and filtering out those that are all 0's.") 
  } else if(! is.null(func_ids) && length(which(! func_ids %in% rownames(func_table))) > 0) {
    stop("Stopping - the above function ids specified in func_ids are not present as rows in the function table after restricting to taxa in the taxa abundance table, and filtering out rows that are all 0's:\n",
         paste(func_ids[which(! func_ids %in% rownames(func_table))], collapse = " "))
  } else if(! is.null(func_ids)) {
    func_table <- func_table[func_ids, , drop = FALSE]
  }
  
  abun_table <- abun_table[colnames(func_table), , drop = FALSE]
  abun_table <- abun_table[, which(colSums(abun_table) > 0), drop = FALSE]
  
  if (ncol(abun_table) == 0) {
    stop("Stopping - no samples remaining after filtering out those that are all 0's.") 
  }
  
  return(list(func = func_table, abun = abun_table))
  
}
