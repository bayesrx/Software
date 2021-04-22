#' Helper function for graphical wavelets filter
sfunc <- function(x,alpha,beta,x1,x2){
  Xmat <- matrix(c(1,1,0,0,x1,x2,x1,x2,x1^2,x2^2,2*x1^2,2*x2^2,x1^3,x2^3,3*x1^3,3*x2^3),4,4)
  a <- solve(Xmat, c(1,1,alpha,beta))
  sfunc <- a[1] + a[2]*x + a[3]*x^2 + a[4]*x^3
  return(sfunc)
}

