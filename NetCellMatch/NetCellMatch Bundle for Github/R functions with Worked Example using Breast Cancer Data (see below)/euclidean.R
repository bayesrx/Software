#'euclidean distance
#'helper function for euclidean distance 
#'@param x1  one of the values to plug into helper function
#'
#'

eucl <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}