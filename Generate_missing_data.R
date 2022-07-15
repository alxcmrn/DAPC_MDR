library(missMethods)
library(adegenet)
library(zvau)
library(stringr)

#Set working directory
setwd("/home/TurnerLab/Pop_Genetics_software/fsc27_linux64/missing_data/")

# Create a folder for each scenario of missing data distribution with XXX subfolders
dir.create(file.path(getwd(), "./random_MD"))
sapply(c("05","20","10"), function(x) dir.create(file.path(getwd(), "population_MD", paste0(x, "_deep_div_no_gene_flow"))))

files<-Sys.glob("/home/TurnerLab/Pop_Genetics_software/fsc27_linux64/no_missing_data/shallow_div_high_gene_flow_50gen/*.no.missing.txt")

# Generate random missing data across full data set
for (mismat in files){
  #df<-read.table(mismat,header=FALSE,colClasses=c('character')) #Import genotype_table to data.frame
  #df<-t(df)
  df_miss_05<-delete_MCAR(df, cols_mis=seq(1,100),p=0.05) #Generates 5% of missing data randomly across the data set.
  df_miss_05<-t(df_miss_05)
  df_miss_20<-delete_MCAR(df, cols_mis=seq(1,100),p=0.2) #Generates 20% of missing data randomly across the data set.
  df_miss_20<-t(df_miss_20)
  df_miss_10<-delete_MCAR(df, cols_mis=seq(1,100),p=0.1) #Generates 10% of missing data randomly across the data set.
  df_miss_10<-t(df_miss_10)
  path_05<-str_replace(mismat,"no_missing_data/shallow_div_high_gene_flow_50gen/","missing_data/random_MD/05_shallow_div_high_gene_flow/")
  path_05<-str_replace(path_05,"no.missing","mismat")
  write.table(df_miss_05,path_05,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  path_20<-str_replace(mismat,"no_missing_data/shallow_div_high_gene_flow_50gen/","missing_data/random_MD/20_shallow_div_high_gene_flow/")
  path_20<-str_replace(path_20,"no.missing","mismat")
  write.table(df_miss_20,path_20,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  path_10<-str_replace(mismat,"no_missing_data/shallow_div_high_gene_flow_50gen/","missing_data/random_MD/10_shallow_div_high_gene_flow/")
  path_10<-str_replace(path_10,"no.missing","mismat")
  write.table(df_miss_10,path_10,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

# Generate missing data biased towards half of the individuals of both populations
files<-Sys.glob("/home/TurnerLab/Pop_Genetics_software/fsc27_linux64/no_missing_data/shallow_div_no_gene_flow/*.no.missing.txt")
for (mismat in files){
  df<-read.table(mismat,header=FALSE,colClasses=c('character')) #Import genotype_table to data.frame
  df<-t(df)
  df_ind_miss_05<-delete_MCAR(df, cols_mis=seq(1,25),p=0.075) #Generates 15% of missing data randomly across each of the 1st 25 individuals in sample 1.
  df_ind_miss_05<-delete_MCAR(df, cols_mis=seq(26,50),p=0.025) #Generates 5% of missing data randomly across each of the 2nd 25 individuals in sample 1.
  df_ind_miss_05<-delete_MCAR(df_ind_miss_05, cols_mis=seq(51,75),p=0.075) #Generates 15% of missing data randomly across each of the 1st 25 individuals in sample 2.
  df_ind_miss_05<-delete_MCAR(df_ind_miss_05, cols_mis=seq(76,100),p=0.025) #Generates 5% of missing data randomly across each of the 2nd 25 individuals in sample 2.
  df_ind_miss_05<-t(df_ind_miss_05)
  path_05<-str_replace(mismat,"no_missing_data/shallow_div_no_gene_flow/","missing_data/individual_MD/05_shallow_div_no_gene_flow/")
  path_05<-str_replace(path_05,"no.missing","mismat")
  write.table(df_ind_miss_05,path_05,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  df_ind_miss_20<-delete_MCAR(df, cols_mis=seq(1,25),p=0.30) #Generates 15% of missing data randomly across each of the 1st 25 individuals in sample 1.
  df_ind_miss_20<-delete_MCAR(df, cols_mis=seq(26,50),p=0.10) #Generates 20% of missing data randomly across each of the 2nd 25 individuals in sample 1.
  df_ind_miss_20<-delete_MCAR(df_ind_miss_20, cols_mis=seq(51,75),p=0.30) #Generates 20% of missing data randomly across each of the 1st 25 individuals in sample 2.
  df_ind_miss_20<-delete_MCAR(df_ind_miss_20, cols_mis=seq(76,100),p=0.10) #Generates 20% of missing data randomly across each of the 2nd 25 individuals in sample 2.
  df_ind_miss_20<-t(df_ind_miss_20)
  path_20<-str_replace(mismat,"no_missing_data/shallow_div_no_gene_flow/","missing_data/individual_MD/20_shallow_div_no_gene_flow/")
  path_20<-str_replace(path_20,"no.missing","mismat")
  write.table(df_ind_miss_20,path_20,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  df_ind_miss_10<-delete_MCAR(df, cols_mis=seq(1,25),p=0.15) #Generates 15% of missing data randomly across each of the 1st 25 individuals in sample 1.
  df_ind_miss_10<-delete_MCAR(df, cols_mis=seq(26,50),p=0.05) #Generates 5% of missing data randomly across each of the 2nd 25 individuals in sample 1.
  df_ind_miss_10<-delete_MCAR(df_ind_miss_10, cols_mis=seq(51,75),p=0.15) ##Generates 15% of missing data randomly across each of the 1st 25 individuals in sample 1.
  df_ind_miss_10<-delete_MCAR(df_ind_miss_10, cols_mis=seq(76,100),p=0.05) ##Generates 5% of missing data randomly across each of the 2nd 25 individuals in sample 1.
  df_ind_miss_10<-t(df_ind_miss_10)
  path_10<-str_replace(mismat,"no_missing_data/shallow_div_no_gene_flow/","missing_data/individual_MD/10_shallow_div_no_gene_flow/")
  path_10<-str_replace(path_10,"no.missing","mismat")
  write.table(df_ind_miss_10,path_10,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

# Generate missing data biased towards one of the populations:
files<-Sys.glob("/home/TurnerLab/Pop_Genetics_software/fsc27_linux64/no_missing_data/shallow_div_no_gene_flow/*.no.missing.txt")
for (mismat in files){
  df<-read.table(mismat,header=FALSE,colClasses=c('character')) #Import genotype_table to data.frame
  df<-t(df)
  df_pop_miss_05<-delete_MCAR(df, cols_mis=seq(1,50),p=0.075) #Generates 5% of missing data randomly across all individuals in sample 1.
  df_pop_miss_05<-delete_MCAR(df_pop_miss_05, cols_mis=seq(51,100),p=0.025) #Generates 5% of missing data randomly across the first 10 individuals in sample 2.
  df_pop_miss_05<-t(df_pop_miss_05)
  path_05<-str_replace(mismat,"no_missing_data/shallow_div_no_gene_flow/","missing_data/population_MD/05_shallow_div_no_gene_flow/")
  path_05<-str_replace(path_05,"no.missing","mismat")
  write.table(df_pop_miss_05,path_05,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  df_pop_miss_20<-delete_MCAR(df, cols_mis=seq(1,50),p=0.3) #Generates 20% of missing data randomly across all individuals in sample 1.
  df_pop_miss_20<-delete_MCAR(df_pop_miss_20, cols_mis=seq(51,100),p=0.1) #Generates 20% of missing data randomly across the first 10 individuals in sample 2.
  df_pop_miss_20<-t(df_pop_miss_20)
  path_20<-str_replace(mismat,"no_missing_data/shallow_div_no_gene_flow/","missing_data/population_MD/20_shallow_div_no_gene_flow/")
  path_20<-str_replace(path_20,"no.missing","mismat")
  write.table(df_pop_miss_20,path_20,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
  df_pop_miss_10<-delete_MCAR(df, cols_mis=seq(1,50),p=0.15)  #Generates 30% of missing data randomly across all individuals in sample 1.
  df_pop_miss_10<-delete_MCAR(df_pop_miss_10, cols_mis=seq(51,100),p=0.05) #Generates 30% of missing data randomly across the first 10 individuals in sample 2.
  df_pop_miss_10<-t(df_pop_miss_10)
  path_10<-str_replace(mismat,"no_missing_data/shallow_div_no_gene_flow/","missing_data/population_MD/10_shallow_div_no_gene_flow/")
  path_10<-str_replace(path_10,"no.missing","mismat")
  write.table(df_pop_miss_10,path_10,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

#Convert files with missing data to matrix with genotypes at a locus coded as 00, 01, 11
setwd("/home/TurnerLab/Pop_Genetics_software/fsc27_linux64/missing_data/individual_biased_MD/")
dirs <- list.dirs()
for (folder in dirs[2:length(dirs)]){
  setwd(folder)
  files<-Sys.glob("*LEA_imputed.mat")
  for (mismat in files){
    df<-read.table(mismat,header=FALSE,na.strings="NA") # Import
    df[df==01]<-12
    df[df==00]<-11
    #df[df==2]<-22
    #df<-t(df)
    path<-str_replace(mismat,"mat","mat.txt")
    write.table(df,file=path,quote=FALSE,sep="\t",na="NA",row.names=FALSE,col.names=FALSE)
  }
  setwd("..")
}


for (gen in 1:length(files)){
  genind<-import2genind(files[gen], ncode=3) # Import genepop file to genind object.
  df<-genind2df(genind)
  pop<-df$pop
  df$pop<-NULL
  df_t<-as.data.frame(t(df))
  df_t_miss<-delete_MCAR(df_t, cols_mis=seq(1,100),p=0.2) #Generates 30% of missing data randomly across both samples.
  #df_t_miss<-delete_MCAR(df_t, cols_mis=seq(1,25),p=0.2) #Generates 30% of missing data randomly across 1st 25 individuals in sample 1.
  #df_t_miss<-delete_MCAR(df_t_miss, cols_mis=seq(51,75),p=0.2) #Generates 30% of missing data randomly across 1st 25 individuals in sample 2.
  df_miss<-as.data.frame(t(df_t_miss))
  df_miss2<-cbind(pop,df_miss)
  genind_miss<-df2genind(df_miss2[,2:2001],ncode=1,type=c("codom"),pop=df_miss2$pop)
  path<-str_replace(files[gen],"deep_div_high_gene_flow","missing_data/deep_div_high_gene_flow_0.2")
  writeGenPop(genind_miss,path,"missing data 0.2")
}

