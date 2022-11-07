#' Fast Weighted UniFrac
#' 
#' Compute pairwise weighted UniFrac distances between samples based on taxa abundances and input tree. Returns matrix with distances in upper triangle.
#' 
#' This is a modified version of the fastUniFrac function taken from v1.39.1 of phyloseq (https://github.com/joey711/phyloseq).
#' This code was distributed under the AGPL-3 license: https://www.gnu.org/licenses/agpl-3.0.en.html.
#'
#' Keys changes:
#' * Separate functions for weighted vs unweighted UniFrac
#' * Distances are always normalized to be within 0 - 1 (i.e., the "normalized" option was removed).
#' * They are compatible with non-phyloseq input objects.
#' * The parallel package (and mclapply) is used rather than foreach, just to be compatible with the rest of FuncDiv.
#' * Matrix of distances rather than a dist object will be returned.
#' 
#' The original function was adapted by the phyloseq authors based on
#' The ISME Journal (2010) 4, 17-27; doi:10.1038/ismej.2009.97;
#' http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html
#' 
#' @param tips_abun data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples.
#' @param tree phylo object, which must contain all row names of `tips_abun` as tip labels.
#' @param ncores integer indicating number of cores to use for parallelizable steps.
#' 
#' @export
fast_weighted_UniFrac <- function(tips_abun, tree, ncores = 1) {
  
  tips_abun <- t(tips_abun)
  
  # Subset tree to be just rownames of intable.
  if (length(which(! rownames(tips_abun) %in% tree$tip.label)) > 0) {
    stop("Stopping - some features in sample are not found as tips in the tree.")
  }
  
  if (length(which(! tree$tip.label %in% rownames(tips_abun))) > 0) {
    tree <- ape::drop.tip(phy = tree,
                          tip = tree$tip.label[which(! tree$tip.label %in% rownames(tips_abun))],
                          trim.internal = TRUE)
  }
  
  
  # Make sure tree is in postorder order for weighted UniFrac.
  tree = ape::reorder.phylo(tree, order = "postorder")
  
  # Re-order rows to match tip label order if needed.
  if (! identical(rownames(tips_abun), tree$tip.label)){
    tips_abun <- tips_abun[tree$tip.label, ]
  }
  
  tips_abun <- as.matrix(tips_abun)
  
  # Create N x 2 matrix of all pairwise combinations of samples.
  sample_combos <- combn(colnames(tips_abun), 2, simplify = FALSE)
  
  ########################################
  # Build the requisite matrices as defined 
  # in the Fast UniFrac article.
  ########################################
  ## This only needs to happen once in a call to UniFrac.
  ## Notice that A and B do not appear in this section.
  # Begin by building the edge descendants matrix (edge-by-sample)
  # `edge_array`
  
  # Create a list of descendants, starting from the first internal node (root)
  ntip <- length(tree$tip.label)
  num_samples <- ncol(tips_abun)
  
  # Create a matrix that maps each internal node to its 2 descendants
  # This matrix doesn't include the tips, so must use node#-ntip to index into it
  node.desc <- matrix(tree$edge[order(tree$edge[, 1]), 2], byrow = TRUE, ncol = 2)
  
  # Define the edge_array object
  # Right now this is a node_array object, each row is a node (including tips)
  # It will be subset and ordered to match tree$edge later
  edge_array <- matrix(0, nrow = ntip + tree$Nnode, ncol = num_samples, 
                       dimnames = list(NULL, sample_names = colnames(tips_abun)))
  
  # Load the tip counts in directly
  edge_array[1:ntip, ] <- as.matrix(tips_abun)
  
  # Get a list of internal nodes ordered by increasing depth
  ord.node <- order(ape::node.depth(tree))[(ntip + 1):(ntip + tree$Nnode)]
  
  # Loop over internal nodes, summing their descendants to get the nodes count
  for (i in ord.node) {
    edge_array[i, ] <- colSums(edge_array[node.desc[i - ntip, ], , drop = FALSE], na.rm = TRUE)
  }
  
  # Keep only those with a parental edge (i.e., drop root) and order to match tree$edge
  edge_array <- edge_array[tree$edge[, 2], , drop = FALSE]
  
  # Remove unneeded variable.
  rm(node.desc)
  
  # For denominator in the normalized distance (i.e., for distances to range from 0 - 1),
  # we need the age of each tip (and due to how this is calculated and returned,
  # postorder formatted tree is needed, which was done above).
  tipAges = ape::node.depth.edgelength(tree)
  
  # Keep only the tips, and add the tip labels in case `z` order differs from `tree`
  tipAges <- tipAges[1:length(tree$tip.label)]

  names(tipAges) <- tree$tip.label
  
  samplesums = colSums(tips_abun)
  
  distlist <- parallel::mclapply(sample_combos,
                       function(i) {
                                    A  <- i[1]
                                    B  <- i[2]
                                    AT <- samplesums[A]
                                    BT <- samplesums[B]

                                    # Weighted UniFrac
                                    wUF_branchweight <- abs(edge_array[, A] / AT - edge_array[, B] / BT)
                                    
                                    # calculate the w-UF numerator
                                    numerator <- sum((tree$edge.length * wUF_branchweight), na.rm = TRUE)
                                    
                                    # denominator (assumes tree-indices and abun table indices are same order)
                                    denominator <- sum((tipAges * (tips_abun[, A] / AT + tips_abun[, B] / BT)), na.rm = TRUE)
                                    
                                    return(numerator / denominator)
                                  },
                       mc.cores = ncores)

  # Initialize UniFracMat with NAs
  UniFracMat <- matrix(NA_real_, num_samples, num_samples)
  rownames(UniFracMat) <- colnames(UniFracMat) <- colnames(tips_abun)

  # Matrix-assign upper-triangle of UniFracMat. Then coerce to dist and return.
  matIndices <- do.call(rbind, sample_combos)

  # Take care of edge case where there are two samples -> 1 pair of indices -> rbind doesn't return a matrix
  if (! is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol = 2)

  UniFracMat[matIndices] <- unlist(distlist)

  return(UniFracMat)

}

#' Fast Unweighted UniFrac
#' 
#' Compute pairwise unweighted UniFrac distances between samples based on taxa abundances and input tree. Returns matrix with distances in upper triangle.
#' 
#' This is a modified version of the fastUniFrac function taken from v1.39.1 of phyloseq (https://github.com/joey711/phyloseq).
#' This code was distributed under the AGPL-3 license: https://www.gnu.org/licenses/agpl-3.0.en.html.
#' function taken from v1.39.1 of phyloseq (https://github.com/joey711/phyloseq).
#' This code was distributed under the AGPL-3 license: https://www.gnu.org/licenses/agpl-3.0.en.html.
#'
#' Keys changes:
#' * Separate functions for weighted vs unweighted UniFrac
#' * Distances are always normalized to be within 0 - 1 (i.e., the "normalized" option was removed).
#' * They are compatible with non-phyloseq input objects.
#' * The parallel package (and mclapply) is used rather than foreach, just to be compatible with the rest of FuncDiv.
#' * Matrix of distances rather than a dist object will be returned.
#' 
#' The original function was adapted by the phyloseq authors based on
#' The ISME Journal (2010) 4, 17-27; doi:10.1038/ismej.2009.97;
#' http://www.nature.com/ismej/journal/v4/n1/full/ismej200997a.html
#' 
#' 
#' @param tips_abun data.frame object containing taxonomic abundances across samples, with rows as taxa and columns as samples.
#' @param tree phylo object, which must contain all row names of `tips_abun` as tip labels.
#' @param ncores integer indicating number of cores to use for parallelizable steps.
#'
#' @export
fast_unweighted_UniFrac <- function(tips_abun, tree, ncores = 1) {
  
  tips_abun <- t(tips_abun)
  
  # Subset tree to be just rownames of intable.
  if (length(which(! rownames(tips_abun) %in% tree$tip.label)) > 0) {
    stop("Stopping - some features in sample are not found as tips in the tree.")
  }
  
  if (length(which(! tree$tip.label %in% rownames(tips_abun))) > 0) {
    tree <- ape::drop.tip(phy = tree,
                             tip = tree$tip.label[which(! tree$tip.label %in% rownames(tips_abun))],
                             trim.internal = TRUE)
  }
  
  # Re-order rows to match tip label order if needed.
  if (! identical(rownames(tips_abun), tree$tip.label)){
    tips_abun <- tips_abun[tree$tip.label, ]
  }
  
  tips_abun <- as.matrix(tips_abun)
  
  # Create N x 2 matrix of all pairwise combinations of samples.
  sample_combos <- combn(colnames(tips_abun), 2, simplify = FALSE)
  
  ########################################
  # Build the requisite matrices as defined 
  # in the Fast UniFrac article.
  ########################################
  ## This only needs to happen once in a call to UniFrac.
  ## Notice that A and B do not appear in this section.
  # Begin by building the edge descendants matrix (edge-by-sample)
  # `edge_array`
  
  # Create a list of descendants, starting from the first internal node (root)
  ntip <- length(tree$tip.label)
  num_samples <- ncol(tips_abun)
  
  # Create a matrix that maps each internal node to its 2 descendants
  # This matrix doesn't include the tips, so must use node#-ntip to index into it
  node.desc <- matrix(tree$edge[order(tree$edge[, 1]), 2], byrow = TRUE, ncol = 2)
  
  # Define the edge_array object
  # Right now this is a node_array object, each row is a node (including tips)
  # It will be subset and ordered to match tree$edge later
  edge_array <- matrix(0, nrow = ntip + tree$Nnode, ncol = num_samples, 
                       dimnames = list(NULL, sample_names = colnames(tips_abun)))
  
  # Load the tip counts in directly
  edge_array[1:ntip, ] <- as.matrix(tips_abun)
  
  # Get a list of internal nodes ordered by increasing depth
  ord.node <- order(ape::node.depth(tree))[(ntip + 1):(ntip + tree$Nnode)]
  
  # Loop over internal nodes, summing their descendants to get the nodes count
  for (i in ord.node) {
    edge_array[i, ] <- colSums(edge_array[node.desc[i - ntip, ], , drop = FALSE], na.rm = TRUE)
  }
  
  # Keep only those with a parental edge (i.e., drop root) and order to match tree$edge
  edge_array <- edge_array[tree$edge[, 2], , drop = FALSE]
  
  # Remove unneeded variable.
  rm(node.desc)
  
  # For unweighted UniFrac, convert the edge_array to an occurrence (presence/absence binary) array
  edge_occ <- (edge_array > 0) - 0
  
  samplesums = colSums(tips_abun)
  
  distlist <- parallel::mclapply(sample_combos,
                       function(i) {
                         A  <- i[1]
                         B  <- i[2]
                         AT <- samplesums[A]
                         BT <- samplesums[B]
                         
                         # Unweighted UniFrac
                         # Subset matrix to just columns A and B
                         edge_occ_AB <- edge_occ[, c(A, B)]
 
                         # Keep only the unique branches. Sum the lengths
                         edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[rowSums(edge_occ_AB, na.rm = TRUE) < 2, ], na.rm = TRUE)

                         # Normalize this sum to the total branches among these two samples, A and B
                         uwUFpairdist <- edge_uni_AB_sum / sum(tree$edge.length[rowSums(edge_occ_AB, na.rm = TRUE) > 0])

                         return(uwUFpairdist)
                       },
                       mc.cores = ncores)
  
  # Initialize UniFracMat with NAs
  UniFracMat <- matrix(NA_real_, num_samples, num_samples)
  rownames(UniFracMat) <- colnames(UniFracMat) <- colnames(tips_abun)
  
  # Matrix-assign upper-triangle of UniFracMat.
  matIndices <- do.call(rbind, sample_combos)
  
  # Take care of edge case where there are two samples -> 1 pair of indices -> rbind doesn't return a matrix
  if (! is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol = 2)
  
  UniFracMat[matIndices] <- unlist(distlist)
  
  return(UniFracMat)
  
}
