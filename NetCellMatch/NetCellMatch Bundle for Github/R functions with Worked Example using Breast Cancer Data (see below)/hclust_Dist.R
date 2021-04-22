#'hclust distances helper function 
#'@param Dist is distance matrix in main function
#'@param dim is number of nodes in laplacian 

hclust.Dist <- function(Dist,dim){
  res <- NULL
  for(j in 1:dim){
    #distance matrix for hierarchical clustering
    Dist[[j]] <- as.dist(Dist[[j]])
    
    res[[j]] <- hclust(Dist[[j]], method = "average")
  }
  return(res)
}
