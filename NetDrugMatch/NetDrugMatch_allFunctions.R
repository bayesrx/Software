#----hclust----#
hclust.Dist = function (Dist, dim) {
  res <- NULL
  for (j in 1:dim) {
    Dist[[j]] <- as.dist(Dist[[j]])
    res[[j]] <- hclust(Dist[[j]], method = "average")
  }
  return(res)
}

#----new_clustering----#
new_clustering=function(pd,scales,Wave){
  # pd number of nodes in Laplacian
  # scales number of scales to consider
  # Wave is the wavelet parameter (run from graph.wavelet)
  Dist = vector("list",scales) #initialize the Distance matrices
  T <- 1000
  R <- matrix(rnorm(pd*T,mean=0,sd=1), pd, T)
  FWT <- vector("list",scales)
  for(i in 1:scales){
    FWT[[i]] <- t(Wave[[i]])%*%R
    Dist[[i]] <- 1 - cor(t(FWT[[i]]))
  }
  distance=c()
  result.cluster = hclust.Dist(Dist = Dist,scales)
  clustmem <- list()
  for(ii in 1:scales){
    cophenetic_dist = cophenetic(result.cluster[[ii]])
    cophenetic_dist = as.matrix(cophenetic_dist)
    y_matrix = matrix(0,nrow = 500, ncol = length(result.cluster[[ii]]$order))
    for (i in 1:length(result.cluster[[ii]]$order)) {
      unique_cophenetic <- c(sort(unique(cophenetic_dist[i,])))
      y = vector()
      x = seq(from = 0, to = max(result.cluster[[ii]]$height),length = 500)
      gap = c(unique_cophenetic,-99) - c(-99, unique_cophenetic)
      gap = gap[2:(length(gap)-1)]
      for (j in 1:(length(gap)-1)){
        y = c(y, rep(gap[j], floor((gap[j]/sum(gap))*500)))
      }
      y = c(y, rep(gap[length(gap)], 500-length(y)))
      y_matrix[,i] = y
    }
    y_average_gap = rowSums(y_matrix)/length(result.cluster[[ii]]$order)
    x_max = x[which(y_average_gap == sort(unique(y_average_gap),decreasing = TRUE)[3])]
    cut_height = median(x_max)
    clustmem[[ii]] = cutree(result.cluster[[ii]],h = cut_height)
  }
  cluster_list = list(result.cluster= result.cluster,clustmem = clustmem)
  return(cluster_list)
}

#----sfunc----#
sfunc = function (x, alpha, beta, x1, x2) {
  Xmat <- matrix(c(1, 1, 0, 0, x1, x2, x1, x2, x1^2, x2^2, 
                   2 * x1^2, 2 * x2^2, x1^3, x2^3, 3 * x1^3, 3 * x2^3), 
                 4, 4)
  a <- solve(Xmat, c(1, 1, alpha, beta))
  sfunc <- a[1] + a[2] * x + a[3] * x^2 + a[4] * x^3
  return(sfunc)
}


#----gfunc----#
gfunc = function (x, alpha, beta, x1, x2) {
  if (sum(x < x1) == 1) 
    gfunc = x1^(-alpha) * x^(alpha)
  if (sum(x1 <= x) * sum(x <= x2) == 1) 
    gfunc = sfunc(x, alpha, beta, x1, x2)
  if (sum(x > x2) == 1) 
    gfunc = x2^(beta) * x^(-beta)
  return(gfunc)
}

#----Graph.wavelet----#
Graph.wavelet = function (dimension, Lp, nscales) {
  p = dimension
  ev <- eigen(Lp)
  ev.val <- ev$value[p:1]
  ev.vec <- ev$vector[, p:1]
  Gs <- list()
  x1 <- 1
  x2 <- x1/ev.val[2]
  alpha <- 2
  beta <- 1/log10(ev.val[3]/ev.val[2])
  s.min = x1/ev.val[2]
  s.max = x1/ev.val[2]^2
  library(sfsmisc)
  s = lseq(from = s.min, to = s.max, length = nscales)
  for (j in 1:length(s)) {
    Gs[[j]] <- matrix(0, p, p)
    for (i in 1:p) {
      Gs[[j]][i, i] = gfunc(s[j] * ev.val[i], alpha, beta, 
                            x1, x2)
    }
  }
  Wave <- list()
  for (j in 1:length(s)) {
    Wave[[j]] <- ev.vec %*% Gs[[j]] %*% t(ev.vec)
  }
  res <- list()
  res <- NULL
  res$Wave <- Wave
  res$s <- s
  res$eval <- ev.val
  return(res)
}

#----get.co-clustering score.matrix----#
get.prop.matrix = function(W){
  # weighted matrix
  print("start0")
  D <- diag(apply(W, 1, sum))
  print("start1")
  L = diag(1,ncol(W)) - solve(D)^(1/2)%*%W%*%solve(D)^(1/2)
  print("start2")
  nfeatures = nrow(L)
  print("start3")
  k = 10 # number of scales = klog2(N)
  nscales = floor(k*log2(nfeatures))
  print("start4")
  wavelet = Graph.wavelet(nfeatures,L,nscales)
  wave = wavelet[[1]]
  print("start_cluster")
  cluster_list = new_clustering(nfeatures, nscales ,wave)
  aggregate_matrix = matrix(0, nrow = nscales, ncol = nfeatures)
  pb1 = txtProgressBar(style=3)
  for (i in 1:nscales) {
    aggregate_matrix[i,] = cluster_list$clustmem[[i]] # clustmem[[i]]
    setTxtProgressBar(pb1,i/nscales)
  }
  close(pb1)
  
  clust_result = diag(nfeatures)
  pb2 = txtProgressBar(style=3)
  for (scan1 in 1:nfeatures) {
    for (scan2 in 1:nfeatures) {
      # Select two col
      vec_compare = aggregate_matrix[,scan1] - aggregate_matrix[,scan2]
      conn_coef = length(which(vec_compare == 0))/length(vec_compare)
      clust_result[scan1,scan2] = conn_coef
    }
    setTxtProgressBar(pb2,scan1/nfeatures)
  }
  close(pb2)
  result_list = list(clust_result = clust_result, aggregate_matrix = aggregate_matrix, dentrograms=cluster_list$result.cluster)
  return(result_list)
}

#----get.adjusted co-clustering score.matrix----#
get.adjusted.prop.matrix = function(stability,usescales,aggregate_matrix){
  nfeatures = ncol(aggregate_matrix)
  clust_result = matrix(1,nrow = nfeatures, ncol = nfeatures)
  aggregate_matrix_use = aggregate_matrix[usescales:nrow(aggregate_matrix),]
  weighted = stability[usescales:length(stability)]/sum(stability[usescales:length(stability)])
  for (scan1 in 1:nfeatures) {
    for (scan2 in 1:nfeatures) {
      conn_coef = 0
      for (scale in 1:nrow(aggregate_matrix_use)) {
        if (aggregate_matrix_use[scale,scan1] == aggregate_matrix_use[scale,scan2]) {
          conn_coef = conn_coef + weighted[scale]
        }
      }
      clust_result[scan1,scan2] = conn_coef
    }
  }
  clust_result = (clust_result-diag(1,nfeatures))*(1/max(clust_result-diag(1,nfeatures)))+diag(1,nfeatures)
  return(clust_result)
}


#----computeR----#
computeR = function(vector1,vector2) {
  a = 0
  b = 0
  Total = length(vector1)*(length(vector2)-1)/2
  count =0
  for (i in 1:(length(vector1)-1)){
    for (j in (i+1):(length(vector2))){
      if (vector1[i]==vector1[j] && vector2[i]==vector2[j]){
        a = a + 1
      }
      if (vector1[i]!=vector1[j] && vector2[i]!=vector2[j]){
        b = b + 1
      }
    }
  }
  (a + b)/Total
}

#----compute_stability----#
compute_stability = function(scalemm){
  vector_rand_index = vector()
  for (i in 1:(nrow(scalemm)-1)){
    for (j in (i+1):nrow(scalemm)){
      vector1 = scalemm[i,]
      vector2 = scalemm[j,]
      Rand_index = computeR(vector1,vector2)
      vector_rand_index = c(vector_rand_index,Rand_index)
    }
  }
  ExpectedR = computeR(scalemm[1,][sample(1:nrow(scalemm))],scalemm[2,][sample(1:nrow(scalemm))])
  MaxIndex = max(vector_rand_index)
  vector_rand_index_adjusted = vector()
  for (i in 1:(nrow(scalemm)-1)){
    for (j in (i+1):nrow(scalemm)){
      vector1 = scalemm[i,]
      vector2 = scalemm[j,]
      Rand_index = computeR(vector1,vector2)
      a = Rand_index - ExpectedR
      b = MaxIndex - ExpectedR
      if (b == 0){
        rand_adjusted = 1
      } else {
        rand_adjusted = (Rand_index - ExpectedR)/(MaxIndex - ExpectedR)
      }
      vector_rand_index_adjusted = c(vector_rand_index_adjusted,rand_adjusted)
    }
  }
  result = sum(vector_rand_index_adjusted)/((nrow(scalemm)-1)*nrow(scalemm)/2)
  if (result>1) {
    result = 1
  }
  return(result)
}