# This file includes the sample code to compute PIDS-Score for lung cancer (two subtypes: LUAD and LUSC).
# Run "NetDrugMatch_allFunctions.R", a file includes all useful functions, before running this file.  
# Inpute datasets for this file include: patient_cellLine_corrMatrix, cellLine_drug_LnIC50. 
# patient_cellLine_corrMatrix: Correlation matrix for patient tumors and cell lines.
# cellLine_drug_LnIC50: Experiment-based drug sensitivity data for cell lines.
# Output dataset is patient_drug_PIDS.
# patient_drug_PIDS: PIDS-Scores for individual patients and drugs.
# Load these datasets from "NetDrugMatch.RData".
# In practice, the correlation matrix of patient-cell line network has multiple choices. We use distance correlation matrix (patient_cellLine_corrMatrix).
# Similarly, the measure of the drug sensitivity for cell lines also has multiple choices. We use LogIC50 (cellLine_drug_logIC50), and the data is from GDSC.
# If you meet any questions or bugs, contact Qingzhi Liu for help (email: qingzliu@umich.edu). 

rm(list=ls(all=TRUE))
load("NetDrugMatch.RData")
# Get Co-Clustering Score Martrix
result_LU_1 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_2 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_3 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_4 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_5 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_6 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_7 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_8 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_9 = get.prop.matrix(patient_cellLine_corrMatrix)
result_LU_10 = get.prop.matrix(patient_cellLine_corrMatrix)

clust_aggragate_list = list()
clust_aggragate_list[[1]] = clust_result1$aggregate_matrix
clust_aggragate_list[[2]] = clust_result2$aggregate_matrix
clust_aggragate_list[[3]] = clust_result3$aggregate_matrix
clust_aggragate_list[[4]] = clust_result4$aggregate_matrix
clust_aggragate_list[[5]] = clust_result5$aggregate_matrix
clust_aggragate_list[[6]] = clust_result6$aggregate_matrix
clust_aggragate_list[[7]] = clust_result7$aggregate_matrix
clust_aggragate_list[[8]] = clust_result8$aggregate_matrix
clust_aggragate_list[[9]] = clust_result9$aggregate_matrix
clust_aggragate_list[[10]] = clust_result10$aggregate_matrix

# Compute Stability
clust_scale_list = list()
for (i in 1:103) {
  clust_scale = matrix(nrow = 10, ncol = 2500)
  for (j in 1:10) {
    clust_scale[j,] = clust_aggragate_list[[j]][i,]
  }
  clust_scale_list[[i]] = clust_scale
}

pb = txtProgressBar(style=3)
adjusted_rand_index = vector()
for (i in 1:103) {
  ari = compute_stability(clust_scale_list[[i]])
  adjusted_rand_index = c(adjusted_rand_index,ari)
  setTxtProgressBar(pb,i/103)
}
close(pb)

# Get Adjusted Co-Clustering Score Matrix
clust_result_adjusted1 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_1$aggregate_matrix)
clust_result_adjusted2 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_2$aggregate_matrix)
clust_result_adjusted3 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_3$aggregate_matrix)
clust_result_adjusted4 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_4$aggregate_matrix)
clust_result_adjusted5 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_5$aggregate_matrix)
clust_result_adjusted6 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_6$aggregate_matrix)
clust_result_adjusted7 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_7$aggregate_matrix)
clust_result_adjusted8 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_8$aggregate_matrix)
clust_result_adjusted9 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_9$aggregate_matrix)
clust_result_adjusted10 = get.adjusted.prop.matrix(adjusted_rand_index,1,result_LU_10$aggregate_matrix)

clust_result_adjusted_final = (clust_result_adjusted1+clust_result_adjusted2+clust_result_adjusted3+clust_result_adjusted4+clust_result_adjusted5+
                                 clust_result_adjusted6+clust_result_adjusted7+clust_result_adjusted8+clust_result_adjusted9+clust_result_adjusted10)/10
clust_result_adjusted_final = (clust_result_adjusted_final-diag(1,1335))*(1/max(max(clust_result_adjusted_final -diag(1,1335))))+diag(1,1335)

# Predict Links between Patients and Drugs (Choose cell lines overlapped in patient-cell line network and cell line-drug network)
patients_cellLines = clust_result_adjusted_final[649:1335,1:648]
cellLines_inter1 = intersect(rownames(patient_cellLine_corrMatrix)[1:648], rownames(cellLine_drug_logIC50))
DrugResponse_IC50 = cellLine_drug_logIC50[cellLines_inter1,] #(Use LN_IC50 from GDSC)
patients_cellLines_use1 = patients_cellLines[,match(cellLines_inter1, rownames(patient_cellLine_corrMatrix)[1:648])]
colnames(patients_cellLines_use1) = cellLines_inter1
rownames(patients_cellLines_use1)[1:362] = paste(substr(rownames(patient_cellLine_corrMatrix)[649:1010],1,12),"LUAD",sep="-")
rownames(patients_cellLines_use1)[363:687] = paste(substr(rownames(patient_cellLine_corrMatrix)[1011:1335],1,12),"LUSC",sep="-")

# Predict Links between Patients and Drugs (Determine Rate Parameters)
cellLines_network = clust_result_adjusted_final[1:648,1:648]
cellLines_network_use = cellLines_network[match(cellLines_inter1, rownames(patient_cellLine_corrMatrix)[1:648]),match(cellLines_inter1, rownames(patient_cellLine_corrMatrix)[1:648])]
sigma_vec = seq(from=0.001,to=0.1,by=0.025)
error_vec = vector()
pb = txtProgressBar(style=3)
for (i in 1:length(sigma_vec)) {
  error_vec = c(error_vec, cellLines_predict_error(sigma_vec[i],cellLines_network_use,DrugResponse_IC50))
  setTxtProgressBar(pb,i/length(sigma_vec))
}
close(pb)
sigma = sigma_vec[which(error_vec==min(error_vec))]

# Predict Links between Patients and Drugs (Locally weighted linear model)
patient_drug_IC50 = matrix(nrow=687,ncol=251)
pb = txtProgressBar(style=3)
for (p in 1:687) {
  for (d in 1:251) {
    sens = 0
    wc_all = vector()
    for (c in 1:293) {
      if (!is.na(DrugResponse_IC50[c,d])) { # Available CellLine-Drug Link
        wc = exp(-(1-patients_cellLines_use1[p,c])^2/(2*sigma))
        wc_all = c(wc_all,wc)
        sens = sens + DrugResponse_IC50[c,d]*wc
      }
    }
    sum_wc_all = sum(wc_all)
    patient_drug_IC50[p,d] = sens/sum_wc_all
  }
  setTxtProgressBar(pb,p/687)
}
close(pb)

# Predict Links between Patients and Drugs (Compute PIDS-Score)
patient_drug_Z = matrix(nrow = 687,ncol = 251)
for (i in 1:251) {
  patient_drug_Z[,i] = (patient_drug_IC50[,i]-mean(DrugResponse_IC50[,i],na.rm = TRUE))/sd(DrugResponse_IC50[,i],na.rm = TRUE)
}
colnames(patient_drug_Z) = colnames(DrugResponse_IC50)
rownames(patient_drug_Z)[1:362] = paste(substr(rownames(patient_cellLine_corrMatrix)[649:1010],1,12),"LUAD",sep="-")
rownames(patient_drug_Z)[363:687] = paste(substr(rownames(patient_cellLine_corrMatrix)[1011:1335],1,12),"LUSC",sep="-")
names(sort(apply(patient_drug_Z[1:362,],2,mean)))
names(sort(apply(patient_drug_Z[363:687,],2,mean)))
patient_drug_PIDS  = -patient_drug_Z # PIDS-Score (Final Result)
rownames(patient_drug_PIDS)[1:362] = paste(substr(rownames(patient_cellLine_corrMatrix)[649:1010],1,12),"LUAD",sep="-")
rownames(patient_drug_PIDS)[363:687] = paste(substr(rownames(patient_cellLine_corrMatrix)[1011:1335],1,12),"LUSC",sep="-")

# Visualization of PIDS-Score
library(gplots)
cancertype = c(rep("green",362),rep("darkorchid",325))
heatmap.2(patient_drug_PIDS , hclustfun=myclust, distfun=mydist, scale="none", dendrogram="both", margins=c(6,12),
          Rowv=TRUE, Colv=TRUE, RowSideColors=cancertype, col=jet.colors(100), trace="none",labRow = FALSE,labCol = FALSE, 
          key.title="Color Key",keysize=0.9, key.par = list(cex=0.5),key.xlab = "PIDS")
legend("topright",legend=c("LUAD","LUSC","BRCA"),title="Cancer Types",
       fill=c("green","darkorchid","brown"), border=FALSE, bty="n", y.intersp =1, cex=0.9)
