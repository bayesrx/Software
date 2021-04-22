#'helper distance function
manhattan <- function(rating1, rating2){
  distance <- abs(rating1-rating2)
  distance <- sum(distance)
  return(distance)
}