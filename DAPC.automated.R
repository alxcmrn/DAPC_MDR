#Load required libraries:
library(adegenet)
library(randomForest)

#Set working directory
setwd("/home/TurnerLab/Pop_Genetics_software/fsc27_linux64/missing_data/population_MD/")
dirs <- list.dirs()
for (folder in dirs[c(2,3,4,5,6,7,8,9,10)]){
  setwd(folder)
  files<-Sys.glob("*LEA_imputed.mat")
  output1<-matrix(NA,100,length(files))
  output2<-matrix(NA,2,length(files))
  output3<-matrix(NA,1,length(files))
  for (mismat in 1:length(files)){
    df<-read.table(files[mismat],header=FALSE, colClasses=c('character')) # Import converted genotype table
    genind<-df2genind(df,ploidy=2,ncode=1,type="codom") # Convert to genind.
    pops<-c("Sample 1","Sample 2") # Define population names.
    genind@pop<-as.factor(rep(pops,c(50,50))) # Assuming that individuals are sorted by population, assign each group of individual to a population.
    scale<-scaleGen(genind,center=TRUE,scale=FALSE,NA.method="mean")
    dapc1<-dapc(scale,center=FALSE,scale=FALSE,n.pca=100,grp=genind$pop,n.da=2) # 1st DAPC (using all PCs based on preliminary runs and DFs).
    a_score<-optim.a.score(dapc1,n.pca=1:ncol(dapc1$tab),smart=TRUE,n=10,plot=FALSE,
                         n.sim=100,n.da=length(levels(dapc1$grp))) # Find optimal number of principal components (PCs) to retain.
    dapc2<-dapc(scale,center=FALSE,scale=FALSE,n.pca=as.numeric(paste0(a_score$best)),grp=genind$pop,n.da=2) # 2nd DAPC using optimal nr of PC found with 'optim.a.score' and the highest DA eigenvalues.
    output1[,mismat]<-dapc2$ind.coord[,1]
    output2[,mismat]<-dapc2$grp.coord[,1]
    output3[,mismat]<-as.numeric(paste0(a_score$best))
  }
  #Add row and column names:
  rownames(output1)<-rownames(dapc2$ind.coord)
  colnames(output1)<-c(1:100)
  rownames(output2)<-rownames(dapc2$grp.coord)
  colnames(output2)<-c(1:100)
  rownames(output3)<-rownames("bestPCs")
  colnames(output3)<-c(1:100)
  # Write table with individual and group PC coordinates:
  write.table(output1,file="ind_coord_mean_LEA.txt",quote=FALSE,sep="\t",col.names=T,row.names=T)
  write.table(output2,file="grp_coord_mean_LEA.txt",quote=FALSE,sep="\t",col.names=T,row.names=T)
  write.table(output3,file="bestPCs_mean_LEA.txt",quote=FALSE,sep="\t",col.names=T,row.names=T)
  setwd("..")
}





# Scatter plot
scatter(dapc2, scree.pca=F,posi.pca="bottomleft",scree.da=F,cstar=0, clab=0, cex=2, cleg=0.9, cex.lab=1.5,
        leg=T, posi.leg="topright",posi.da="bottomright")
