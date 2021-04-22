#'Circos Plot
#'
#'This function creates a circos plot visualizing results
#' @param  cutoff point for connections to keep
#' @param  Hetpatient matrix of connections to assess 
#'
#'

circos=function(cutoff,Hetpatient){
d=rowSums(Hetpatient)
index=which(d>=cutoff)
mat3=Hetpatient[c(index),]
c=colSums(mat3)
index2=which(c!=0)
mat3=mat3[,c(index2)]
rownames(mat3) = paste0("CL", 1:nrow(mat3))
colnames(mat3) = paste0("P", 1:ncol(mat3))
mycols=c('red','green','black','blue4','chartreuse')

circos.par(gap.after = 2) #rep(c(rep(1, 4), 8), 22)) 
col_fun = colorRamp2(c(0,1,2,5,9), mycols, transparency = 0.5)
chordDiagram(mat3, col=col_fun, annotationTrack = c("grid", "axis"), preAllocateTracks = list( track.height = uh(4, "mm"), 
                                                                                               track.margin = c(uh(4, "mm"), 0)
))

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, niceFacing = TRUE)
}, bg.border = NA)

highlight.sector(rownames(mat3), track.index = 1, col = "red", text = "Cell Lines", cex = 0.8, text.col = "white", niceFacing = TRUE)
highlight.sector(colnames(mat3), track.index = 1, col = "green", text = "Patients", cex = 0.8, text.col = "white", niceFacing = TRUE)
}

#Het=matrix(runif(100,0,1),10,10)
#circos(5,Het)

