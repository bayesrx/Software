#' Cutting dendogram
#'
#' This helper function cuts chooses an "optimal" dendogram for the clustering pipeline
#' @param scale  the scale being considered
#'

cut_height <- function(scale,pd){
  #find distances where dendrogram nodes are formed
  node_dist <- cophenetic(result.cluster[[scale]])
  node_dist <- as.matrix(node_dist)
  min_val_cut <- max_val_cut <- rep(0, pd)
  for(i in 1:pd){
    #find the unique sorted distances
    unique_node_dist <- c(-99, sort(unique(node_dist[i,])))
    
    gap <- rep(0, pd)
    
    for(j in 1:pd){
      gap[j] <- node_dist[i,j] - unique_node_dist[which(unique_node_dist == node_dist[i,j]) - 1]
    }
    gap[gap == 99] <- 0
    dist_val <- Dist[[scale]][i, which(gap == max(gap))]
    min_val_cut[i] <- min(dist_val)
    max_val_cut[i] <- max(dist_val)
    
  }
  return(list(min_val_cut = min_val_cut, max_val_cut = max_val_cut))
}