##Author: Alex Cameron

library(vcfR)
library(hierfstat)
library(adegenet)
library(parallel)

popmap<- read.table("./popmap.txt")
popmap$V2<- as.factor(popmap$V2)

vcf.files<- list.files(path = ".",
                       pattern = ".vcf.gz",
                       full.names = F)


pFST <- mclapply(vcf.files, function(x){
  temp<- vcfR::read.vcfR(file = x)
  temp.genind<- vcfR2genind(temp, sep = "[|/]")
  temp.genind@pop <- popmap$V2
  temp.heirfstat <- genind2hierfstat(temp.genind)
  temp.heirfstat[,1]<- as.integer(temp.heirfstat[,1])
  mtemp.heirfstat<- temp.heirfstat[order(temp.heirfstat[,1]),]
  pairfst<- pairwise.WCfst(temp.heirfstat)
  matNames<- levels(pop(temp.genind))
  dimnames(pairfst)<- list(matNames, matNames)
  pairfst[upper.tri(pairfst)] <- NA
  pfst_melt<- reshape2::melt(pairfst)
  pfst_melt$V4<- paste(sep = '_', pfst_melt$Var1, pfst_melt$Var2)
  pfst_melt<- pfst_melt[ ,c(1,2,4,3)]
  pfst_melt<-pfst_melt[ ,-c(1:2)]
  colnames(pfst_melt)<- c("pop_comp", "fst")
  pfst_melt<- na.omit(pfst_melt)
  return(pfst_melt)
}, mc.cores = 16L)
rm(popmap)

vcf.names<- gsub("\\..*","",vcf.files)
pFST<- setNames(pFST, vcf.names)
list2env(pFST, envir = .GlobalEnv)
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
concat_data<- as.data.frame(do.call(rbind, dfs))
concat_data$md<- row.names(concat_data); row.names(concat_data) <- NULL
##output file should have 3 X k(k-1)/2 number of lines +1 for header line; where k is number of populations
write.table(concat_data, file = "./output_pfst.txt", quote = F, sep = '\t', row.names = F)















