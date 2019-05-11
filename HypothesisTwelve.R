rm(list = objects()); ls()

library(MatrixModels)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(asreml)
library(beepr)
library(rrBLUP)
library(Hmisc)
library(broom)
library(MVN)
library(MASS)
library(car)
library(ape)
library(RColorBrewer)
library(caret)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

snpChip <- read_delim(
  "./Genotype_Database/90KsnpChipHapMap/AMsnpChipImputed.hmp.txt", 
  "\t", escape_double = FALSE, trim_ws = TRUE)
snpChip<- snpChip %>% 
  clean_names()

missAmbiguous = c('0', '+', '-')
hetCodes = c('R','Y','S','W','K','M','B','D','H','V')
hapgeno=as.matrix(snpChip[,13:ncol(snpChip)])
hapgeno[hapgeno %in% missAmbiguous]=NA
hapgeno[hapgeno=='N']=NA
hapgeno[hapgeno %in% hetCodes]='H'
snpChip=cbind(snpChip[,1:12], hapgeno)
rm(hapgeno)

write.table(snpChip, file="./Genotype_Database/SelectedImputedBeagle.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./Genotype_Database/SelectedImputedBeagle.txt", 
                header=TRUE, check.names=F, sep = "\t")

snpChip[snpChip == snpChip$allele_a] = 1
snpChip[snpChip == snpChip$allele_b] = -1
snpChip[snpChip == "H"] = 0
snpChip[snpChip == "C"] = NA
snpChip[snpChip == "A"] = NA
snpChip[snpChip == "T"] = NA
snpChip[snpChip == "G"] = NA
snpChip[snpChip == "-"] = NA
snpChip[snpChip == "."] = NA

snpChip<- snpChip[ ,c(1,4,5,13:311)]

write.table(snpChip, 
            file="./Genotype_Database/SelectedImputedBeagleNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./Genotype_Database/SelectedImputedBeagleNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")
snpMatrix<- t(snpChip[ , c(-1, -2, -3)])
snpMatrix %>% glimpse()
snpMatrix[1:10,1:10]

snpMatrix<- setDT(as.data.frame(snpMatrix),keep.rownames = T)
snpMatrix[1:10,1:10]

##### Phenotypes ####

pheno17<- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18<- fread("./Phenotype_Database/pheno18_htpLong.txt")
phenoLong<- fread("./Phenotype_Database/Pheno_Long1718.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(phenoLong)

mean(pheno17$GRYLD)
Matrix::mean(pheno18$GRYLD)

phenoLong<- phenoLong %>% 
  dplyr::rename(Plot_ID = entity_id) %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column)

pheno17$Date<- as.Date(pheno17$Date,format = "%Y-%m-%d")
pheno17$Date<- format(pheno17$Date, "%Y%m%d")
pheno18$Date<- as.Date(pheno18$Date,format = "%Y-%m-%d")
pheno18$Date<- format(pheno18$Date, "%Y%m%d")

pheno17<- pheno17 %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20170331:RedEdge_20170609) %>% 
  tidylog::inner_join(phenoLong) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse()

pheno18<- pheno18 %>% 
  dplyr::rename(Plot_ID = entity_id)  %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20171120:RE_20180613) %>% 
  tidylog::inner_join(phenoLong) %>% 
  tidylog::select(-starts_with("height_")) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse()

pheno17$Plot_ID<- as.factor(pheno17$Plot_ID)
pheno17$Variety<- as.factor(pheno17$Variety)
pheno17$block<- as.factor(pheno17$block)
pheno17$rep<- as.factor(pheno17$rep)
pheno17$range<- as.factor(pheno17$range)
pheno17$column<- as.factor(pheno17$column)

pheno18$Plot_ID<- as.factor(pheno18$Plot_ID)
pheno18$Variety<- as.factor(pheno18$Variety)
pheno18$block<- as.factor(pheno18$block)
pheno18$rep<- as.factor(pheno18$rep)
pheno18$range<- as.factor(pheno18$range)
pheno18$column<- as.factor(pheno18$column)

asreml.license.status()

##### Trialing something ####

set.seed(1962)

par(mar=c(1,1,1,1))
mean(pheno17$GRYLD)

t17<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,#+ vm(Variety,relMat,singG="PSD")
             #mef = list(relMat = snpMatrix17),
             #na.method(y = "include",x = "include"),
             data = pheno17)
summary(t17)
plot(t17)
blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat17<- blues %>% 
  rename(GRYLD = effect) %>% 
  glimpse() 

##### Making it a function for all of the VI's 2017

effectvars <- names(pheno17) %in% c("block", "rep", "Variety", "year", 
                                    "column","range", "Plot_ID","GRYLD")
traits <- colnames(pheno17[ , !effectvars])
traits
fieldInfo<- pheno17 %>% 
  tidylog::select(Variety, rep, block, column, range)

for (i in traits) {
  print(paste("Working on trait", i))
  j<- i
  
  data<- cbind(fieldInfo, pheno17[,paste(i)])
  names(data)<- c("Variety","rep","block","column","range","Trait")
  print(colnames(data))
  
  t17<- asreml(fixed = Trait ~ 0 + Variety,
               random = ~ rep + rep:block,
               data = data)
  pdf(paste0("./Figures/AsremlPlots/ASREML_Blues17_",
             i,".pdf"))
  plot(t17)
  
  blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat17<- blues %>% 
    inner_join(dat17)
  
}

#dat17[1:5,1:15]

dev.off()
graphics.off()

beep()

##### Making it a function for all of the VI's 2018

par(mar=c(1,1,1,1))

t18<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,
             data = pheno18)
plot(t18)
blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat18<- blues %>% 
  rename(GRYLD = effect) %>% 
  glimpse() 

effectvars <- names(pheno18) %in% c("block", "rep", "Variety", "year", 
                                    "column","range", "Plot_ID","GRYLD")
traits <- colnames(pheno18[ , !effectvars])
traits
fieldInfo<- pheno18 %>% 
  tidylog::select(Variety, rep, block, column, range)

for (i in traits) {
  print(paste("Working on trait", i))
  j<- i
  
  data<- cbind(fieldInfo, pheno18[,paste(i)])
  names(data)<- c("Variety","rep","block","column","range","Trait")
  print(colnames(data))
  
  t18<- asreml(fixed = Trait ~ 0 + Variety,
               random = ~ rep + rep:block,
               data = data)
  pdf(paste0("./Figures/AsremlPlots/ASREML_Blues18_",
             i,".pdf"))
  plot(t17)
  
  blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat18<- blues %>% 
    inner_join(dat18)
  dev.off()
}

dat18[1:5,1:15]

beep(2)


###############################################################
####                    rrBlup trial                      ####

snpMatrix[1:5,1:5]

dat17<- dat17 %>% 
  semi_join(snpMatrix, by = "rn")

snpMatrix17<- snpMatrix %>% 
  semi_join(dat17,by = "rn") 
rownames(snpMatrix17) <- snpMatrix17[,1]
snpMatrix17[,1] <- NULL

snpMatrix17<- as.matrix(snpMatrix17)

dat18<- dat18 %>% 
  semi_join(snpMatrix, by = "rn")

snpMatrix18<- snpMatrix %>% 
  semi_join(dat18,by = "rn") 
rownames(snpMatrix18) <- snpMatrix18[,1]
snpMatrix18[,1] <- NULL

snpMatrix18<- as.matrix(snpMatrix18)

##### Defining the training and test populations ####
#define the training and test populations
#training-80% validation-20%
train17= as.matrix(sample(1:nrow(dat17), (0.4 * nrow(dat17))))
test17= as.matrix(setdiff(1:nrow(dat17),train17))
train17 = dat17[train17]
train17Lines<- as.data.frame(train17)
test17 = dat17[test17]
Pheno_train17=dat17 %>% 
  semi_join(train17Lines, by = c("rn" = "train17"))
m_train17=snpMatrix17[train17,]
Pheno_valid17=dat17 %>% 
  anti_join(train17Lines, by = c("rn" = "train17"))
m_valid17=snpMatrix17[test17,]

##### Predicting Phenotypes ####
## GRYLD
yield=(Pheno_train17[,"GRYLD"])
covar<- (Pheno_train17[,"NDRE_20170505"])
yield_answer<-mixed.solve(yield, Z = m_train17, #X = covar, 
                          SE = TRUE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid17 %*% e
pred_yield=(pred_yield_valid[,1]) + yield_answer$beta
pred_yield
yield_valid = Pheno_valid17[,"GRYLD"]
YLD_accuracy <-cor(pred_yield_valid, yield_valid )
YLD_accuracy

## VI
ndvi0512=(Pheno_train17[,"NDVI_20170512"])
ndvi0512_answer<-mixed.solve(ndvi0512, Z = m_train17, SE = TRUE)
NDVI0512 = ndvi0512_answer$u
e = as.matrix(NDVI0512)
pred_ndvi0512_valid =  m_valid17 %*% e
pred_ndvi0512=(pred_ndvi0512_valid[,1]) + ndvi0512_answer$beta
pred_ndvi0512
ndvi0512_valid = Pheno_valid17[,"NDVI_20170512"]
NDVI_20170512_accuracy <-cor(pred_ndvi0512_valid, ndvi0512_valid )
NDVI_20170512_accuracy

re0512=(Pheno_train17[,"RedEdge_20170512"])
re0512_answer<-mixed.solve(re0512, Z = m_train17, SE = TRUE)
RE0512 = re0512_answer$u
e = as.matrix(RE0512)
pred_re0512_valid =  m_valid17 %*% e
pred_re0512=(pred_re0512_valid[,1]) + re0512_answer$beta
pred_re0512
re0512_valid = Pheno_valid17[,"RedEdge_20170512"]
RedEdge_20170512_accuracy <-cor(pred_re0512_valid, re0512_valid )
RedEdge_20170512_accuracy

##### Cross-Validation ####

effectvars <- names(dat17) %in% c("rn")
traits <- colnames(dat17[ , !effectvars])
traits
cycles=2000
accuracy = matrix(nrow=cycles, ncol=length(traits))
colnames(accuracy)<- traits

for(r in 1:cycles) {
  print(paste("Rep cycle: ",r))
  train= as.matrix(sample(1:nrow(dat17), (0.4 * nrow(dat17))))
  test= as.matrix(setdiff(1:nrow(dat17),train17))
  train = dat17[train]
  trainLines<- as.data.frame(train)
  test = dat17[test]
  Pheno_train=dat17 %>% 
    semi_join(trainLines, by = c("rn" = "train"))
  m_train=snpMatrix17[train,]
  Pheno_valid=dat17 %>% 
    anti_join(trainLines, by = c("rn" = "train"))
  m_valid=snpMatrix17[test17,]
  
  for (i in traits) {
    print(paste(i))
    trait=(Pheno_train[,paste(i)])
    trait_answer<-mixed.solve(trait, Z=m_train, SE = F)
    TRT = trait_answer$u
    e = as.matrix(TRT)
    pred_trait_valid =  m_valid %*% e
    pred_trait = (pred_trait_valid[,1]) + trait_answer$beta
    pred_trait
    trait_valid = Pheno_valid[,paste(i)]
    accuracy[r,paste(i)] <-cor(pred_trait_valid, trait_valid, use="complete" )
  }
  
}

write.table(accuracy, "./R/rrBlup/HypothesisTwelve/GenomicSelection_40_2000_accuracy17.txt", quote = F, 
            sep = "\t",row.names = F,col.names = T)
#accuracy<- as.data.frame(accuracy[1:129,])
mean(accuracy$GRYLD, na.rm = T)
sd(accuracy$GRYLD, na.rm = T)



## 2018
train18= as.matrix(sample(1:nrow(dat18), (0.4 * nrow(dat18))))
test18= as.matrix(setdiff(1:nrow(dat18),train18))
train18 = dat18[train18]
train18Lines<- as.data.frame(train18)
test18 = dat18[test18]
Pheno_train18=dat18 %>% 
  semi_join(train18Lines, by = c("rn" = "train18"))
m_train18=snpMatrix18[train18,]
Pheno_valid18=dat18 %>% 
  anti_join(train18Lines, by = c("rn" = "train18"))
m_valid18=snpMatrix18[test18,]

##### Predicting Phenotypes ####
## GRYLD
yield=(Pheno_train18[,"GRYLD"])
yield_answer<-mixed.solve(yield, Z = m_train18, SE = TRUE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid18 %*% e
pred_yield=(pred_yield_valid[,1]) + yield_answer$beta
pred_yield
yield_valid = Pheno_valid18[,"GRYLD"]
YLD_accuracy <-cor(pred_yield_valid, yield_valid )
YLD_accuracy

##### Cross-Validation ####

effectvars <- names(dat18) %in% c("rn")
traits <- colnames(dat18[ , !effectvars])
traits
cycles=2000
accuracy = matrix(nrow=cycles, ncol=length(traits))
colnames(accuracy)<- traits

for(r in 1:cycles) {
  print(paste("Rep cycle: ",r))
  train= as.matrix(sample(1:nrow(dat18), (0.4 * nrow(dat18))))
  test= as.matrix(setdiff(1:nrow(dat18),train))
  train = dat18[train]
  trainLines<- as.data.frame(train)
  test = dat18[test]
  Pheno_train=dat18 %>% 
    semi_join(trainLines, by = c("rn" = "train"))
  m_train=snpMatrix18[train,]
  Pheno_valid=dat18 %>% 
    anti_join(trainLines, by = c("rn" = "train"))
  m_valid=snpMatrix18[test,]
  
  for (i in traits) {
    print(paste(i))
    trait=(Pheno_train[,paste(i)])
    trait_answer<-mixed.solve(trait, Z=m_train, SE = F)
    TRT = trait_answer$u
    e = as.matrix(TRT)
    pred_trait_valid =  m_valid %*% e
    pred_trait = (pred_trait_valid[,1]) + trait_answer$beta
    print(pred_trait)
    trait_valid = Pheno_valid[,paste(i)]
    accuracy[r,paste(i)] <-cor(pred_trait_valid, trait_valid, use="complete" )
  }
  
}

write.table(accuracy, 
            "./R/rrBlup/HypothesisTwelve/GenomicSelection_40_2000_accuracy18.txt", 
            quote = F,
            sep = "\t",row.names = F,col.names = T)

###### Covariate added in? #####
traits = c("GRYLD","GRYLDCovar")
cycles = 2000
accuracyCovar<- matrix(nrow = cycles, ncol = length(traits))
colnames(accuracyCovar) <- traits

for (i in 1:cycles) {
  print(paste("number of cycles:", i))
  train= as.matrix(sample(1:nrow(dat17), (0.4 * nrow(dat17))))
  test= as.matrix(setdiff(1:nrow(dat17),train))
  train = dat17[train]
  trainLines<- as.data.frame(train)
  test = dat17[test]
  Pheno_train=dat17 %>% 
    semi_join(trainLines, by = c("rn" = "train"))
  m_train=snpMatrix17[train,]
  Pheno_valid=dat17 %>% 
    anti_join(trainLines, by = c("rn" = "train"))
  m_valid=snpMatrix17[test,]
  
  covar<- (Pheno_train17[,"NDRE_20170505"])
  yield=(Pheno_train17[,"GRYLD"])
  
  yield=(Pheno_train17[,"GRYLD"])
  yield_answer<-mixed.solve(yield, Z=m_train17,
                            K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid17[,"GRYLD"]
  accuracyCovar[i,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
  yield_answer<-mixed.solve(yield, Z=m_train17, X = covar,
                            K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid17[,"GRYLD"]
  accuracyCovar[i,2] <-cor(pred_yield_valid, yield_valid, use="complete" )
  print(accuracyCovar[i,])
  
}
colMeans(accuracyCovar)
sd(accuracyCovar[,1])
sd(accuracyCovar[,2])
cor(accuracyCovar)
t.test(accuracyCovar)

write.table(accuracyCovar,"./R/rrBlup/HypothesisTwelve/GSGRYLD17_Covar_40_1000_accuracy.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

effectvars <- names(dat17) %in% c("rn")
traits <- colnames(dat17[ , !effectvars])
cycles = 2000
accuracyCovar<- matrix(nrow = cycles, ncol = 1)
colnames(accuracyCovar) <- "Cycles"
accuracyCovar<- as.data.frame(accuracyCovar)
accuracyCovar$Cycles<- 1:cycles
CovariateFile<- matrix(nrow = cycles, ncol = length(traits))
colnames(CovariateFile) <- traits

for (r in 1:cycles) {
  print(paste("number of cycles:", r))
  
  
  train= as.matrix(sample(1:nrow(dat17), (0.4 * nrow(dat18))))
  test= as.matrix(setdiff(1:nrow(dat17),train))
  train = dat17[train]
  trainLines<- as.data.frame(train)
  test = dat17[test]
  Pheno_train=dat17 %>% 
    semi_join(trainLines, by = c("rn" = "train"))
  m_train=snpMatrix17[train,]
  Pheno_valid=dat17 %>% 
    anti_join(trainLines, by = c("rn" = "train"))
  m_valid=snpMatrix17[test,]
  
  for (i in traits) {
    print(paste(i))
    corDat17<- rcorr(as.matrix(dat17[,-1]))
    corDat17<- flattenCorrMatrix(corDat17$r,corDat17$P)
    covar<- corDat17 %>%
      tidylog::filter(row == paste(i) | column == paste(i)) %>%
      tidylog::filter(cor == max(cor))
    covariate<- if_else(covar$row == paste(i), paste(covar$column),
                        paste(covar$row))
    print(covariate)
    CovariateFile[r,c(paste(i))]<- covariate
    
    trait=(Pheno_train[,paste(i)])
    covar=(Pheno_train[,covariate])
    
    trait_answer<-mixed.solve(trait, Z=m_train, X = covar, SE = F)
    TRT = trait_answer$u
    e = as.matrix(TRT)
    pred_trait_valid =  m_valid %*% e
    pred_trait = (pred_trait_valid[,1]) + trait_answer$beta
    print(pred_trait)
    trait_valid = Pheno_valid[,paste(i)]
    accuracyCovar[r,paste(i)]<- cor(pred_trait_valid, trait_valid,
                                    use="complete" )
    print(accuracyCovar[r,paste(i)])
  }
}

Means<-as.data.frame(colMeans(accuracyCovar))
cor(accuracyCovar)
t.test(accuracyCovar)

write.table(accuracyCovar,"./R/rrBlup/HypothesisTwelve/GenomicSelection_Covar_40_2000_accuracy17.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
write.ftable(CovariateFile, "./R/rrBlup/HypothesisTwelve/CovarCorrFile17_40_2000.txt", quote = F,
             sep = "\t", col.names = T, row.names = F)


effectvars <- names(dat18) %in% c("rn")
traits <- colnames(dat18[ , !effectvars])
cycles = 2000
accuracyCovar<- matrix(nrow = cycles, ncol = 1)
colnames(accuracyCovar) <- "Cycles"
accuracyCovar<- as.data.frame(accuracyCovar)
accuracyCovar$Cycles<- 1:cycles
CovariateFile<- matrix(nrow = cycles, ncol = length(traits))
colnames(CovariateFile) <- traits

for (r in 1:cycles) {
  print(paste("number of cycles:", r))
  
  
  train= as.matrix(sample(1:nrow(dat18), (0.4 * nrow(dat18))))
  test= as.matrix(setdiff(1:nrow(dat18),train))
  train = dat18[train]
  trainLines<- as.data.frame(train)
  test = dat18[test]
  Pheno_train=dat18 %>% 
    semi_join(trainLines, by = c("rn" = "train"))
  m_train=snpMatrix18[train,]
  Pheno_valid=dat18 %>% 
    anti_join(trainLines, by = c("rn" = "train"))
  m_valid=snpMatrix18[test,]
  
  for (i in traits) {
    print(paste(i))
    corDat18<- rcorr(as.matrix(dat18[,-1]))
    corDat18<- flattenCorrMatrix(corDat18$r,corDat18$P)
    covar<- corDat18 %>%
      tidylog::filter(row == paste(i) | column == paste(i)) %>%
      tidylog::filter(cor == max(cor))
    covariate<- if_else(covar$row == paste(i), paste(covar$column),
                        paste(covar$row))
    print(covariate)
    CovariateFile[r,c(paste(i))]<- covariate
    
    trait=(Pheno_train[,paste(i)])
    covar=(Pheno_train[,covariate])
    
    trait_answer<-mixed.solve(trait, Z=m_train, X = covar, SE = F)
    TRT = trait_answer$u
    e = as.matrix(TRT)
    pred_trait_valid =  m_valid %*% e
    pred_trait = (pred_trait_valid[,1]) + trait_answer$beta
    print(pred_trait)
    trait_valid = Pheno_valid[,paste(i)]
    accuracyCovar[r,paste(i)]<- cor(pred_trait_valid, trait_valid,
                                    use="complete" )
    print(accuracyCovar[r,paste(i)])
  }
}

Means<-as.data.frame(colMeans(accuracyCovar))
cor(accuracyCovar)
t.test(accuracyCovar)

write.table(accuracyCovar,"./R/rrBlup/HypothesisTwelve/GenomicSelection_Covar_40_2000_accuracy18.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)
write.ftable(CovariateFile, "./R/rrBlup/HypothesisTwelve/CovarCorrFile18_40_2000.txt", quote = F,
             sep = "\t", col.names = T, row.names = F)



###############################################################################
##### Examining results ####
getwd()

fileNames<- list.files(path = "./R/rrBlup/HypothesisTwelve",
                       full.names = T,
                       pattern = "accuracy17.txt$")
traitNames<- basename(fileNames) %>%
  str_remove_all("_accuracy17.txt")

load.file<- function (filename) {
  d<- fread(file = filename,header = TRUE,check.names = F,data.table = F)
  d
}

data17<- lapply(fileNames, load.file)

names(data17)<- traitNames

data17<- data17 %>% 
  map(function(x) setDT(x, keep.rownames = T)) %>% 
  map(function(x) gather(x,key = "Phenotype", 
                         value = "correlation",
                         RedEdge_20170609:GRYLD)) 
data17<- data17 %>% 
  map(function(x) group_by(x, Phenotype)) %>% 
  map(function(x) dplyr::summarise(x, average = mean(correlation),
                                   standardDev = sd(correlation),
                                   standardError = sd(correlation)/sqrt(length(correlation)))) 
dats17<- data.frame()

for (i in traitNames) {
  dats17<- data17[[paste(i)]] %>% 
    mutate(Trial = paste(i)) %>% 
    bind_rows(dats17)
}

fileNames<- list.files(path = "./R/rrBlup/HypothesisTwelve",
                       full.names = T,
                       pattern = "accuracy18.txt$")
traitNames<- basename(fileNames) %>%
  str_remove_all("_accuracy18.txt")

data18<- lapply(fileNames, load.file)

names(data18)<- traitNames

data18<- data18 %>% 
  map(function(x) setDT(x, keep.rownames = T)) %>% 
  map(function(x) gather(x,key = "Phenotype", 
                         value = "correlation",
                         RE_20180613:GRYLD)) 
data18<- data18 %>% 
  map(function(x) group_by(x, Phenotype)) %>% 
  map(function(x) 
    dplyr::summarise(x, average = mean(correlation),
                     standardDev = sd(correlation),
                     standardError = sd(correlation)/sqrt(length(correlation)))) 
dats18<- data.frame()

for (i in traitNames) {
  dats18<- data18[[paste(i)]] %>% 
    mutate(Trial = paste(i)) %>% 
    bind_rows(dats18)
}

dats17<- dats17 %>% 
  separate(Phenotype, c("Trait","date"), sep = "_") 
dats17$date<- as.Date(dats17$date, format = "%Y%m%d")
dats18<- dats18 %>% 
  separate(Phenotype, c("Trait","date"), sep = "_") 
dats18$date<- as.Date(dats18$date, format = "%Y%m%d")

dats17$Trait<- as.factor(dats17$Trait)
dats17$Trial<- as.factor(dats17$Trial)

dats17<- dats17  %>% 
  mutate(Covariate = str_detect(Trial,"Covar")) %>% 
  separate(Trial,c("GS","Trial"), sep = "n_") %>% 
  tidylog::select(-GS) 
dats17$Trial<- str_replace(dats17$Trial,"Covar_","")
dats17$Split<- dats17$Trial
dats17<- dats17 %>% 
  separate(Split,c("Split","Rep"), sep = "_")

dats18<- dats18  %>% 
  mutate(Covariate = str_detect(Trial,"Covar")) %>% 
  separate(Trial,c("GS","Trial"), sep = "n_") %>% 
  tidylog::select(-GS) 
dats18$Trial<- str_replace(dats18$Trial,"Covar_","")
dats18$Split<- dats18$Trial
dats18<- dats18 %>% 
  separate(Split,c("Split","Rep"), sep = "_")

dats17 %>% 
  filter(Trait =="GRVI") %>% 
  filter(Split == 80) %>% 
  ggplot(aes(x = date, y = average, 
             colour = Trial, 
             shape = Covariate)) +
  geom_point() +
  theme_bw() +
  scale_x_date(breaks = "1 week",date_labels = "%m/%d") +
  geom_errorbar(aes(ymin = average - standardError, 
                    ymax = average + standardError), 
                alpha = 0.75) +
  facet_wrap(~Rep,ncol = 2) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.y = 
          element_line(colour = "#bdbdbd",
                       linetype = 2)) +
  scale_colour_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                 '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                 '#cab2d6','#6a3d9a','#b15928')) +
  labs(title = "Accuracy of Prediction",
       subtitle = "2016/2017 GRVI 80:20")

dats17 %>% 
  filter(Trait == "GRYLD") %>% 
  ggplot(aes(x = Trial, y = average, shape = Covariate)) +
  geom_point() +
  geom_errorbar(aes(ymin = average - standardError, 
                    ymax = average + standardError)) +
  theme(axis.text = element_text(colour = "black")) +
  labs(title = "Accuracy of Prediction",
       subtitle = "2016/2017 GRYLD") +
  theme_bw()

dats18$Trait<- as.factor(dats18$Trait)
dats18$Trial<- as.factor(dats18$Trial)

dats18 %>% 
  filter(Trait =="GNDVI") %>% 
  filter(Split == 80) %>% 
  ggplot(aes(x = date, y = average, 
             colour = Trial, 
             shape = Covariate)) +
  geom_point() +
  theme_bw() +
  scale_x_date(breaks = "1 week",date_labels = "%m/%d") +
  geom_errorbar(aes(ymin = average - standardError, 
                    ymax = average + standardError), 
                alpha = 0.75) +
  facet_wrap(~Rep,ncol = 2) +  
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.y = 
          element_line(colour = "#bdbdbd",
                       linetype = 2)) +
  scale_colour_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                 '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                 '#cab2d6','#6a3d9a','#ffff66','#b15928')) +
  labs(title = "Accuracy of Prediction",
       subtitle = "2017/2018 GNDVI 80:20")

dats18 %>% 
  filter(Trait == "GRYLD") %>% 
  ggplot(aes(x = Trial, y = average, shape = Covariate)) +
  geom_point() +
  geom_errorbar(aes(ymin = average - standardError, 
                    ymax = average + standardError)) +
  theme(axis.text = element_text(colour = "black")) +
  labs(title = "Accuracy of Prediction",
       subtitle = "2017/2018 GRYLD") +
  theme_bw()
###############################################################################
