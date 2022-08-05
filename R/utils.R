#' @export
contrib_to_multitab <- function(contrib_tab,
                                samp_colname = "sample",
                                func_colname = "function.",
                                abun_colname = "taxon_abun",
                                taxon_colname = "taxon",
                                copy.num_colname = "genome_function_count") {
  
  if (class(contrib_tab) != "data.frame") {
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
    stopping("Stopping - at least one combination of samples and taxa has a different specificied abundance across different functions. This could indicate an issue with the contributional format in the first place, and means that it is impossible to convert to the two table format.")
  }
  
  contrib_taxa_abun_wide <- data.frame(data.table::dcast.data.table(data = data.table(contrib_taxa_abun),
                                                                    formula = as.formula(paste(taxon_colname, "~", samp_colname, sep = " ")),
                                                                    value.var = abun_colname), check.names = FALSE)
  
  rownames(contrib_taxa_abun_wide) <- as.character(contrib_taxa_abun_wide[, taxon_colname])
  contrib_taxa_abun_wide <- contrib_taxa_abun_wide[, -1]
  contrib_taxa_abun_wide[is.na(contrib_taxa_abun_wide)] <- 0
  
  
  # Function copy number table.
  contrib_tab_func_copy_num <- contrib_tab[, c(func_colname, taxon_colname, copy.num_colname)]
  
  contrib_tab_func_copy_num <- contrib_tab_func_copy_num[which(! duplicated(contrib_tab_func_copy_num)), ]
  
  duplicated_func_taxa_combos_i <- which(duplicated(contrib_tab_func_copy_num[, c(func_colname, taxon_colname)]))
  
  if (length(duplicated_func_taxa_combos_i) > 0) {
    stopping("Stopping - at least one combination of taxa / functions has a different specificied copy number across samples. This can happen for instance if pathways levels per taxon are computed based on how much they contribute to the community-wide pathway abundance. This means that it is impossible to convert to the two table format.")
  }
  
  contrib_tab_func_copy_num_wide <- data.frame(data.table::dcast.data.table(data = data.table(contrib_tab_func_copy_num),
                                                                            formula = as.formula(paste(func_colname, "~", taxon_colname, sep = " ")),
                                                                            value.var = copy.num_colname), check.names = FALSE)
  
  rownames(contrib_tab_func_copy_num_wide) <- as.character(contrib_tab_func_copy_num_wide[, func_colname])
  contrib_tab_func_copy_num_wide <- contrib_tab_func_copy_num_wide[, -1]
  contrib_tab_func_copy_num_wide[is.na(contrib_tab_func_copy_num_wide)] <- 0
  
  return(list(taxon_abun = contrib_taxa_abun_wide,
              function_copy_num = contrib_tab_func_copy_num_wide))
  
}


# Convert from multi-table format to contributional format.
#' @export
multitab_to_contrib <- function(func_tab,
                                abun_tab,
                                ncores = 1,
                                samp_colname = "sample",
                                func_colname = "function.",
                                abun_colname = "taxon_abun",
                                taxon_colname = "taxon",
                                copy.num_colname = "genome_function_count") {

  if (class(func_tab) != "data.frame") {
    stop("Stopping - \"func_tab\" input argument must be a data.frame.")
  }
  
  if (class(abun_tab) != "data.frame") {
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
                                         
                                         contrib_block <- reshape2::melt(data = abun_tab_subset,
                                                                         id = taxon_colname,
                                                                         variable.name = samp_colname,
                                                                         value.name = abun_colname)
                                         
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

#' @export
func_abun_crossproduct <- function(in_abun, in_func) {
  
  # Check that all rows are found in function table.
  if(length(which(! rownames(in_abun) %in% rownames(in_func))) > 0) {
    stop("Stoppings - some rows in abundance table not found in function table.")
  }
  
  in_func <- in_func[rownames(in_abun), ]
  
  return(data.frame(t(crossprod(as.matrix(in_abun), as.matrix(in_func))), check.names = FALSE))
  
}


# Get intersecting taxa in function and abundance table.
# Remove rows and columns that are all 0's in both tables.
# Do some basic sanity checks at key steps during this process.
# Can specify a subset of func ids that are expected to be present after filtering, otherwise will only throw error if no functions remain.
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
  } else if(! is.null(func_ids) && nrow(func_table) < length(func_ids)) {
    stop("Stopping - the above function ids specified in func_ids are not present as rows in the function table after restricting to taxa in the taxa abundance table, and filtering out rows that are all 0's:\n",
         paste(func_ids[which(! func_ids %in% rownames(func_table))], collapse = " "))
  }
  
  abun_table <- abun_table[colnames(func_table), , drop = FALSE]
  abun_table <- abun_table[, which(colSums(abun_table) > 0), drop = FALSE]
  
  if (ncol(abun_table) == 0) {
    stop("Stopping - no samples remaining after filtering out those that are all 0's.") 
  }
  
  return(list(func = func_table, abun = abun_table))
  
}
