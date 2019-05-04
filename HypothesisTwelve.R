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
library(ClassDiscovery)
library(ggdendro)
library(pvclust)
library(viridis)
library(RColorBrewer)

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

write.table(snpChip, file="./Genotype_Database/SelectedImputedBeagleNumeric.txt",
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
  semi_join(dat17,by = "rn") %>% 
  tidylog::select(-rn)

snpMatrix17<- as.matrix(snpMatrix17)

dat18<- dat18 %>% 
  semi_join(snpMatrix, by = "rn")

snpMatrix18<- snpMatrix %>% 
  semi_join(dat18,by = "rn") %>% 
  tidylog::select(-rn)

snpMatrix18<- as.matrix(snpMatrix18)

##### Defining the training and test populations ####
#define the training and test populations
#training-80% validation-20%
train17= as.matrix(sample(1:nrow(dat17), (0.8 * nrow(dat17))))
test17<-setdiff(1:96,train17)
Pheno_train17=dat17[train17,]
m_train17=snpMatrix17[train17,]
Pheno_valid17=dat17[test17,]
m_valid17=snpMatrix17[test17,]

##### Predicting Phenotypes ####
## GRYLD
yield=(Pheno_train17[,"GRYLD"])
yield_answer<-mixed.solve(yield, Z = m_train17, SE = TRUE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid17 %*% e
pred_yield=(pred_yield_valid[,1]) + yield_answer$beta
pred_yield
yield_valid = Pheno_valid17[,"GRYLD"]
YLD_accuracy <-rcorr(pred_yield_valid, yield_valid )
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
NDVI_20170512_accuracy <-cor(pred_yield_valid, yield_valid )
NDVI_20170512_accuracy

re0512=(Pheno_train17[,"RedEdge_20170512"])
re0512_answer<-mixed.solve(re0512, Z = m_train17, SE = TRUE)
RE0512 = re0512_answer$u
e = as.matrix(RE0512)
pred_re0512_valid =  m_valid17 %*% e
pred_re0512=(pred_re0512_valid[,1]) + re0512_answer$beta
pred_re0512
re0512_valid = Pheno_valid17[,"RedEdge_20170512"]
RedEdge_20170512_accuracy <-cor(pred_yield_valid, yield_valid )
RedEdge_20170512_accuracy

##### Cross-Validation ####
