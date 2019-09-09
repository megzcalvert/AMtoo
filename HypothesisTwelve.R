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
set.seed(1642)

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

dat17<- fread("./Phenotype_Database/Hyp10BLUEs_17.txt")
dat18<- fread("./Phenotype_Database/Hyp10BLUEs_18.txt")
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
Pheno_train17<- dat17 %>% 
  dplyr::sample_frac(0.8)
Pheno_valid17<- dat17 %>% 
  anti_join(Pheno_train17, by = "rn")

m_train17=snpMatrix17[Pheno_train17$rn,]
m_valid17=snpMatrix17[Pheno_valid17$rn,]

##### Predicting Phenotypes ####
## GRYLD
yield=(Pheno_train17[,"GRYLD"])
covar<- (Pheno_train17[,"NDRE_20170505"])
yield_answer<-mixed.solve(yield, Z = m_train17, X = covar,
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
cycles=100
accuracy = matrix(nrow=cycles, ncol=length(traits))
colnames(accuracy)<- traits

for(r in 1:cycles) {
  print(paste("Rep cycle: ",r))
  
  Pheno_train<- dat17 %>%
    dplyr::sample_frac(0.8)
  Pheno_valid<- dat17 %>%
    anti_join(Pheno_train, by = "rn")
  
  m_train=snpMatrix17[Pheno_train$rn,]
  m_valid=snpMatrix17[Pheno_valid$rn,]
  
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

write.csv(accuracy, 
          "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy17.txt", 
          quote = F,row.names = F)
colMeans(accuracy)
mean(accuracy$GRYLD, na.rm = T)
sd(accuracy$GRYLD, na.rm = T)


# 2018
Pheno_train18<- dat18 %>% 
  dplyr::sample_frac(0.8)
Pheno_valid18<- dat18 %>% 
  anti_join(Pheno_train18, by = "rn")

m_train18=snpMatrix18[Pheno_train18$rn,]
m_valid18=snpMatrix18[Pheno_valid18$rn,]

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
cycles=100
accuracy = matrix(nrow=cycles, ncol=length(traits))
colnames(accuracy)<- traits

for(r in 1:cycles) {
  print(paste("Rep cycle: ",r))
  
  Pheno_train<- dat18 %>%
    dplyr::sample_frac(0.8)
  Pheno_valid<- dat18 %>%
    anti_join(Pheno_train, by = "rn")
  
  m_train=snpMatrix18[Pheno_train$rn,]
  m_valid=snpMatrix18[Pheno_valid$rn,]
  
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

write.csv(accuracy, 
          "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy18.txt", 
          quote = F,row.names = F)
colMeans(accuracy)
mean(accuracy$GRYLD, na.rm = T)
sd(accuracy$GRYLD, na.rm = T)


##### Predicting across years #####

Pheno_train17<- dat17 %>% 
  tidylog::semi_join(dat18, by = "rn")
Pheno_valid18<- dat18 %>% 
  tidylog::semi_join(dat17, by = "rn")
Pheno_train18<- dat18 %>% 
  tidylog::semi_join(dat17, by = "rn")
Pheno_valid17<- dat17 %>% 
  tidylog::semi_join(dat18, by = "rn")

m_train17=snpMatrix17[Pheno_train17$rn,]
m_valid18=snpMatrix18[Pheno_valid18$rn,]
m_train18=snpMatrix18[Pheno_train18$rn,]
m_valid17=snpMatrix17[Pheno_valid17$rn,]

yield=(Pheno_train17[,"GRYLD"])
yield_answer<-mixed.solve(yield, Z=m_train17,
                          K=NULL, SE = FALSE, return.Hinv=FALSE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid18 %*% e
pred_yield=(pred_yield_valid[,1])+yield_answer$beta
pred_yield

yield_valid = Pheno_valid18[,"GRYLD"]
cor(pred_yield_valid, yield_valid, use="complete" )

yield=(Pheno_train18[,"GRYLD"])
yield_answer<-mixed.solve(yield, Z=m_train18,
                          K=NULL, SE = FALSE, return.Hinv=FALSE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid17 %*% e
pred_yield=(pred_yield_valid[,1])+yield_answer$beta
pred_yield
yield_valid = Pheno_valid17[,"GRYLD"]
cor(pred_yield_valid, yield_valid, use="complete" )

## Flatten a correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# effectvars <- names(dat17) %in% c("rn")
# traits <- colnames(dat17[ , !effectvars])
# cycles = 1000
# accuracyCovar<- matrix(nrow = cycles, ncol = 1)
# colnames(accuracyCovar) <- "Cycles"
# accuracyCovar<- as.data.frame(accuracyCovar)
# accuracyCovar$Cycles<- 1:cycles
# CovariateFile<- matrix(nrow = cycles, ncol = length(traits))
# colnames(CovariateFile) <- traits
# 
# for (r in 1:cycles) {
#   print(paste("number of cycles:", r))
# 
# 
# dplyr::sample_frac(0.8)
# Pheno_valid17<- dat17 %>% 
#   anti_join(train, by = "rn")

# m_train17=snpMatrix17[Pheno_train17$rn,]
# m_valid17=snpMatrix17[Pheno_valid17$rn,]
# 
#   for (i in traits) {
#     print(paste(i))
#     corDat17<- rcorr(as.matrix(dat17[,-1]))
#     corDat17<- flattenCorrMatrix(corDat17$r,corDat17$P)
#     covar<- corDat17 %>%
#       tidylog::filter(row == paste(i) | column == paste(i)) %>%
#       tidylog::filter(cor == max(cor))
#     covariate<- if_else(covar$row == paste(i), paste(covar$column),
#                         paste(covar$row))
#     print(covariate)
#     CovariateFile[r,c(paste(i))]<- covariate
# 
#     trait=(Pheno_train[,paste(i)])
#     covar=(Pheno_train[,covariate])
# 
#     trait_answer<-mixed.solve(trait, Z=m_train, X = covar, SE = F)
#     TRT = trait_answer$u
#     e = as.matrix(TRT)
#     pred_trait_valid =  m_valid %*% e
#     pred_trait = (pred_trait_valid[,1]) + trait_answer$beta
#     pred_trait
#     trait_valid = Pheno_valid[,paste(i)]
#     accuracyCovar[r,paste(i)]<- cor(pred_trait_valid, trait_valid,
#                                     use="complete" )
#     print(accuracyCovar[r,paste(i)])
#   }
# }
# 
# Means<-as.data.frame(colMeans(accuracyCovar))
# cor(accuracyCovar)
# t.test(accuracyCovar)
# 
# write.table(accuracyCovar,"./R/rrBlup/HypothesisTwelve/GenomicSelection_Covar_60_2000_accuracy17.txt",
#             quote = F, sep = "\t", col.names = T, row.names = F)
# write.table(CovariateFile, "./R/rrBlup/HypothesisTwelve/CovarCorrFile17_60_2000.txt", quote = F,
#              sep = "\t", col.names = T, row.names = F)


# effectvars <- names(dat18) %in% c("rn")
# traits <- colnames(dat18[ , !effectvars])
# cycles = 1000
# accuracyCovar<- matrix(nrow = cycles, ncol = 1)
# colnames(accuracyCovar) <- "Cycles"
# accuracyCovar<- as.data.frame(accuracyCovar)
# accuracyCovar$Cycles<- 1:cycles
# CovariateFile<- matrix(nrow = cycles, ncol = length(traits))
# colnames(CovariateFile) <- traits

# for (r in 1:cycles) {
#   print(paste("number of cycles:", r))
#   
#   
#   Pheno_train18<- dat18 %>% 
#     dplyr::sample_frac(0.8)
#   Pheno_valid18<- dat18 %>% 
#     anti_join(train, by = "rn")
#   
#   m_train18=snpMatrix18[Pheno_train18$rn,]
#   m_valid18=snpMatrix18[Pheno_valid18$rn,]
#   
#   for (i in traits) {
#     print(paste(i))
#     corDat18<- rcorr(as.matrix(dat18[,-1]))
#     corDat18<- flattenCorrMatrix(corDat18$r,corDat18$P)
#     covar<- corDat18 %>%
#       tidylog::filter(row == paste(i) | column == paste(i)) %>%
#       tidylog::filter(cor == max(cor))
#     covariate<- if_else(covar$row == paste(i), paste(covar$column),
#                         paste(covar$row))
#     print(covariate)
#     CovariateFile[r,c(paste(i))]<- covariate
#     
#     trait=(Pheno_train[,paste(i)])
#     covar=(Pheno_train[,covariate])
#     
#     trait_answer<-mixed.solve(trait, Z=m_train, X = covar, SE = F)
#     TRT = trait_answer$u
#     e = as.matrix(TRT)
#     pred_trait_valid =  m_valid %*% e
#     pred_trait = (pred_trait_valid[,1]) + trait_answer$beta
#     pred_trait
#     trait_valid = Pheno_valid[,paste(i)]
#     accuracyCovar[r,paste(i)]<- cor(pred_trait_valid, trait_valid,
#                                     use="complete" )
#     print(accuracyCovar[r,paste(i)])
#   }
# }

# Means<-as.data.frame(colMeans(accuracyCovar))
# cor(accuracyCovar)
# t.test(accuracyCovar)
# 
# write.table(accuracyCovar,"./R/rrBlup/HypothesisTwelve/GenomicSelection_Covar_40_1000_accuracy18.txt",
#             quote = F, sep = "\t", col.names = T, row.names = F)
# write.table(CovariateFile, "./R/rrBlup/HypothesisTwelve/CovarCorrFile18_40_1000.txt", quote = F,
#             sep = "\t", col.names = T, row.names = F)



###############################################################################
#### Examining results ####
getwd()

gs17_50<- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_50_accuracy17.txt",
  header = T)
gs17_50<- gs17_50 %>% 
  mutate(Rep = row.names(gs17_50)) %>% 
  gather(key = "Trait", value = "Accuracy",RedEdge_20170609:GRYLD) %>% 
  # group_by(Trait) %>% 
  # summarise(average = mean(Accuracy), standardDeviation = sd(Accuracy)) %>% 
  separate(Trait, c("Trait","Date"), sep = "_") %>% 
  filter(Trait != "GRYLD")

gs18_50<- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_50_accuracy18.txt",
  header = T)
gs18_50<- gs18_50 %>% 
  mutate(Rep = row.names(gs18_50)) %>% 
  gather(key = "Trait", value = "Accuracy",RE_20180613:GRYLD)  %>% 
  # group_by(Trait) %>% 
  # summarise(average = mean(Accuracy), standardDeviation = sd(Accuracy)) %>% 
  separate(Trait, c("Trait","Date"), sep = "_") %>% 
  filter(Trait != "GRYLD")

gs17_100<- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy17.txt",
  header = T)
gs17_100<- gs17_100 %>% 
  mutate(Rep = row.names(gs17_100)) %>% 
  gather(key = "Trait", value = "Accuracy",RedEdge_20170609:GRYLD) %>% 
  # group_by(Trait) %>% 
  # summarise(average = mean(Accuracy), standardDeviation = sd(Accuracy)) %>% 
  separate(Trait, c("Trait","Date"), sep = "_") %>% 
  filter(Trait != "GRYLD")

gs18_100<- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy18.txt",
  header = T)
gs18_100<- gs18_100 %>% 
  mutate(Rep = row.names(gs18_100)) %>% 
  gather(key = "Trait", value = "Accuracy",RE_20180613:GRYLD) %>% 
  # group_by(Trait) %>% 
  # summarise(average = mean(Accuracy), standardDeviation = sd(Accuracy)) %>% 
  separate(Trait, c("Trait","Date"), sep = "_") %>% 
  filter(Trait != "GRYLD")

gs17_100$Date<- as.Date(gs17_100$Date, format = "%Y%m%d")
gs18_100$Date<- as.Date(gs18_100$Date, format = "%Y%m%d")

ggplot(gs18_100, aes(x = Date, y = Accuracy)) + 
  geom_boxplot(aes(group = Date)) +
  facet_wrap(~Trait, ncol = 2, scales = "fixed") + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  #coord_cartesian(ylim = c(0,1)) +
  #scale_x_date(labels = "%m%d") +
  labs(title = "Genomic Prediction Accuracies",
       subtitle = "2017/2018 season",
       y = "Average prediction accuracy")


