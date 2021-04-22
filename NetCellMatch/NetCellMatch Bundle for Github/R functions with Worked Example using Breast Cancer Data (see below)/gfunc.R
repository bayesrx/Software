#' Graphical Wavelets helper function
#' graphcial wavelet function, the high band pass filter as detailed in literature
gfunc <- function(x,alpha, beta, x1, x2){
  if( sum(x < x1) == 1) gfunc = x1^(-alpha)*x^(alpha)
  if( sum(x1<=x)*sum(x<=x2) == 1) gfunc = sfunc(x,alpha,beta,x1,x2)
  if( sum(x > x2) == 1) gfunc = x2^(beta)*x^(-beta)      #beta gets screwed up because of how its defined
  return(gfunc) 
}