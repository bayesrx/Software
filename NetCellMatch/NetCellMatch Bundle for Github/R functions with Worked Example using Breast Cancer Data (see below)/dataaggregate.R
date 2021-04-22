#' Useful data aggregation
#'
#' This function aggregates data in a usable form outputting a barplat and returning a Countermatrix of coclustering
#' @param pd number of nodes in laplacian 
#' @param scales number of scales to consider
#' @param clustmem is cluster assignments at each scale (a list output of clustering function )



data.aggregate=function(pd,scales,clustmem){
#Make a matrix cataloguing responses for each 
scaleclusterspatient=matrix(NA, nrow=pd,ncol=scales)
for(i in 1:scales){
  scaleclusterspatient[,i]=clustmem[[i]]
}


#Get number of clusters at different scales 2 ways of summarizing 
clustersatscalepatients=c()
for(i in 1:scales){
  clustersatscalepatients[i]=length(unique(scaleclusterspatient[,i]))
}

uniqueclusterspatients=list(0)
for(i in 1:scales){
  uniqueclusterspatients[[i]]=unique(scaleclusterspatient[,i])
}

#Name trunccells as whole thing, because no non-unique clusters 
trunccells=clustmem

#Compute Number of Times Scales Come Together 
matchpatient=list()
matchpatient_=list()
for(i in 1:scales){
  for(j in 1:pd){
    matchpatient[[j]]=as.vector(which(trunccells[[i]] %in% j))
  }
  matchpatient_[[i]]=matchpatient
}

#Make plots indicating cluster differences 
plot(clustersatscalepatients,type="l",xlab="scale",ylab="# of clusters",main="# of clusters by scale",col="blue",xaxt='n')
axis(1,at=seq(1,150,by=5))
#-------------------------------------------------------------------------------
#PART 5
#MAKE GIANT COUNTER MATRIX
Countpatient=matrix(0,nrow=pd,ncol=pd)
for(i in 1:scales){
  for(j in 1:pd){
    if(length(matchpatient_[[i]][[j]])!=0){
      for(k in 1:length(matchpatient_[[i]][[j]])){
        if(k!=length(matchpatient_[[i]][[j]])){
          Countpatient[matchpatient_[[i]][[j]][[k]],matchpatient_[[i]][[j]][[k+1]]]=Countpatient[matchpatient_[[i]][[j]][[k]],matchpatient_[[i]][[j]][[k+1]]]+1
          Countpatient[matchpatient_[[i]][[j]][[k+1]],matchpatient_[[i]][[j]][[k]]]=Countpatient[matchpatient_[[i]][[j]][[k+1]],matchpatient_[[i]][[j]][[k]]]+1        
        }
      }
    }
  }
  cat("Scale Completed:", i," ")
}
return(Countpatient)
}

