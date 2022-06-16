
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
