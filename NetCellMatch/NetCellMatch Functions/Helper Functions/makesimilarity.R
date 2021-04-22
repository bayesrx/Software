#' Function to make similarity matrix (can make on own)

make.similarity <- function(my.data, similarity) {
  N <- ncol(my.data) #note change here I had made column 
  S <- matrix(rep(NA,N^2), ncol=N)
  #print(dim(S))
  if(similarity=='spear'){
    S=cor(t(my.data),method='spear')
    S=abs(S)
  }
  else{
    for(i in 1:N) {
      for(j in 1:N) {
        S[i,j] <- similarity(my.data[i,], my.data[j,])
        #print(S[i,j])
        print(j)
      }
    }}
  return(S)
}