#' Make laplacian
#' Function to make laplacian takes in input is data 

make.laplacian=function(cellpatient){
  D <- diag(apply(cellpatient, 1, sum)) # sum rows
  D[1:8,1:8]
  
  #Compute Unnormalized Graph Laplacian 
  U <- D - cellpatient
  #round(U[1:12,1:12],1)
  Ue=solve(D)
  dim(Ue)
  #check normalized laplacian 
  Lp=Ue^(1/2)%*%as.matrix(U)%*%Ue^(1/2)  
  return(Lp)
}