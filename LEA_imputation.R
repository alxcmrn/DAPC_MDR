######################
### Empirical Data ###
######################

library("LEA")

### 5% Missing Data ###

gila1<- geno2lfmm(input.file = "GilaTrout_05per.geno")
gila1.snmf<- snmf("GilaTrout_05per.lfmm", 
     K = 1:5, 
     entropy = TRUE,
     percentage = 0.05,
     repetitions = 50,
     project = "new",
     CPU=8
)

best<- which.min(cross.entropy(gila1.snmf, K = 5))

impute(gila1.snmf, "GilaTrout_05per.geno", method = 'mode', K = 5, run = best)

### 10% Missing Data ###

gila2<- geno2lfmm(input.file = "GilaTrout_10per.geno")
gila2.snmf<- snmf("GilaTrout_10per.lfmm", 
                  K = 1:5, 
                  entropy = TRUE,
                  percentage = 0.05,
                  repetitions = 50,
                  project = "new",
                  CPU=8
)

best<- which.min(cross.entropy(gila2.snmf, K = 5))
impute(gila2.snmf, "GilaTrout_10per.geno", method = 'mode', K = 5, run = best)

### 20% Missing Data ###

gila3<- geno2lfmm(input.file = "GilaTrout_20per.geno")
gila3.snmf<- snmf("GilaTrout_20per.lfmm", 
                  K = 1:5, 
                  entropy = TRUE,
                  percentage = 0.05,
                  repetitions = 50,
                  project = "new",
                  CPU=8
)

best<- which.min(cross.entropy(gila3.snmf, K = 5))
impute(gila3.snmf, "GilaTrout_20per.geno", method = 'mode', K = 5, run = best)

######################
### Simulated Data ###
######################

genofiles<- as.list(list.files(path = ".", pattern = ".geno", full.names = T))

convert<- mclapply(genofiles, function(x){geno2lfmm(input.file = x)}, mc.cores = 4L)

lfmm.files<- as.list(list.files(path = ".", pattern = ".lfmm", full.names = T))

impute_sim<- function (x, y){ 
  xsnmf<- snmf(input.file = x, K = 2, entropy = TRUE, percentage = 0.05, repetitions = 10, project = "new")
  best <- which.min(cross.entropy(xsnmf, K = 2))
  impute(xsnmf, input.file = y, method = 'mode', run = best)
}

mapply(impute_sim, x = genofiles, y = lfmm.files)