#IDENTIFY AVATAR CELL LINES AND R SCORES 

Avatar=function(Counter,cindex,pindex,scales,cutoff){
  Hetpatient=Counter[cindex,pindex]
  #Divide By Scales 
  Hetpatient=Hetpatient/scales
  Hetpatient=(1/max(Hetpatient))*Hetpatient
  Z=rowSums(Hetpatient)
  Av=which(Z>=cutoff)
  return(list(Av,Z,Hetpatient))
}
