#' @export
calc_func_abun <- function(in_abun, in_func, ncores=1) {
  
  out_df <- data.frame(matrix(NA, nrow=ncol(in_func), ncol=ncol(in_abun)))
  colnames(out_df) <- colnames(in_abun)
  rownames(out_df) <- colnames(in_func)
  
  # Check that all rows are found in function table.
  if(length(which(! rownames(in_abun) %in% rownames(in_func))) > 0) {
    stop("Stoppings - some rows in abundance table not found in function table.")
  }
  
  in_func <- in_func[rownames(in_abun), ]
  
  out_sample_func_abun <- mclapply(colnames(in_abun), function(x) { return(colSums(in_abun[, x] * in_func)) }, mc.cores=ncores)
  names(out_sample_func_abun) <- colnames(in_abun)
  
  for(sample in colnames(in_abun)) {
    out_df[, sample] <- out_sample_func_abun[[sample]]
  }
  
  return(out_df)
  
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
