#NetCellMatch Example 

#----------------------------------------------------
#STEP 1: LOAD DATA, FORM SIMILARITY MATRIX AND LAPLACIAN
#----------------------------------------------------
Breast=read.csv('Breast_Cancer_BatchCorrected.csv')
S=make.similarity(Breast[,5:223],'spear')
Lap=make.laplacian(S)
#-----------------------------------------------------

#-----------------------------------------------------
#STEP 2:RUN GRAPHICAL WAVELETS ACROSS DESIRED SCALES 
#-----------------------------------------------------
#RUNNING ACROSS 50 SCALES 
Waves=Graph.wavelet(ncol(Lap),Lap,50)
#-----------------------------------------------------

#-----------------------------------------------------
#STEP 3: COMPUTE DISTANCE CORRELATION MATRIX AND CLUSTER RESULTS AT EACH SCALE
#-----------------------------------------------------
Distance=clustering(ncol(Lap),50,Waves$Wave)
#-----------------------------------------------------

#----------------------------------------------------
#STEP 4: AGGREGATE CLUSTER RESULTS ACROSS SCALES 
#-----------------------------------------------------
Aggregation=data.aggregate(ncol(Lap),50,Distance)
#-----------------------------------------------------

#----------------------------------------------------
#STEP 5: GIVEN GIANT COUNTER MATRIX, IDENTFY AVATAR CELL LINES
#----------------------------------------------------
Avatars=Avatar(Aggregation,which(Breast$Ind==2),which(Breast$Ind==1),
               scales=50,cutoff=0.75)

#RETURNS A LIST WITH (IN ORDER): AVATAR CELL LINE INDEX, Z SCORE, R SCORE MATRIX











