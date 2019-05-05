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
#pheno17$GRYLD<- as.numeric(pheno17$GRYLD)

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

# normalisedPheno17<-mvn(pheno17[,3:51], univariateTest = "SW", desc = T)
# normalisedPheno17$univariateNormality
# normalisedPheno17$multivariateNormality
# normalisedPheno17$Descriptives
# normalisedPheno18<-mvn(pheno18[,3:93], univariateTest = "SW", desc = T) 
# normalisedPheno18$univariateNormality
# normalisedPheno18$multivariateNormality
# normalisedPheno18$Descriptives

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

pheno17<- pheno17 %>% 
  semi_join(snpMatrix, by = c("Variety" = "rn"))

snpMatrix17<- snpMatrix %>% 
  semi_join(pheno17,by = c("rn" = "Variety")) 
snpMatrix17[1:5,1:5]

snpLines<- snpMatrix17$rn
snpMatrix17<- snpMatrix17 %>% 
  tidylog::select(-rn)
relMat<- A.mat(snpMatrix17)
colnames(relMat)<- snpLines
rownames(relMat)<- snpLines
relMat[1:5,1:5]
str(relMat)
relMat<- ginv(relMat)
relMat[1:5,1:5]
str(relMat)
rownames(relMat)<- snpLines
colnames(relMat)<- snpLines

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

par(mar=c(1,1,1,1))


##### Select in 2016/2017 and see what happens in 2017/2018 ####
# Using NDVI_20170512 for the selection
dat17 %>% 
  ggplot(aes(x = NDVI_20170512)) +
  geom_density() +
  geom_vline(xintercept = mean(dat17$NDVI_20170512), linetype = 2) +
  theme_bw()

ndvi0512_top<- dat17 %>% 
  dplyr::arrange(desc(NDVI_20170512)) %>% #use desc() if you want highest to lowest
  tidylog::select(rn,GRYLD,NDVI_20170512) %>% 
  glimpse()

0.05*nrow(dat17)
ndvi0512_top<- ndvi0512_top[1:round(0.05*nrow(dat17)),]

dat17 %>% 
  ggplot(aes(x = NDVI_20170512)) +
  geom_density() +
  geom_vline(xintercept = mean(dat17$NDVI_20170512), linetype = 2) +
  geom_vline(xintercept = mean(ndvi0512_top$NDVI_20170512), colour = "blue",
             linetype = 2) +
  geom_vline(xintercept = min(ndvi0512_top$NDVI_20170512), colour = "blue") +
  theme_bw() +
  labs(title = "Distribution of NDVI_20170512", 
       subtitle = "All lines - black, top 5% of NDVI_20170512 - blue")

dat17 %>% 
  ggplot(aes(x = GRYLD)) +
  geom_density() +
  geom_density(data = ndvi0512_top,aes(x = GRYLD), colour = "blue") +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2) +
  geom_vline(xintercept = mean(ndvi0512_top$GRYLD), linetype = 2, 
             colour = "blue") +
  #geom_vline(xintercept = min(RE0512_top$GRYLD), colour = "blue") +
  #geom_vline(xintercept = max(RE0512_top$GRYLD), colour = "blue") +
  theme_bw() +
  labs(title = "Distribution of GRYLD in 2016/2017", 
       subtitle = "All lines - black, lines selected from NDVI_20170512 - blue")

ndviSelect<- ndvi0512_top %>% 
  left_join(dat18, by = "rn") %>% 
  tidylog::select(rn,GRYLD.x,NDVI_20170512,GRYLD.y) %>% 
  glimpse

dat18 %>% 
  ggplot(aes(x = GRYLD)) +
  geom_density() +
  geom_density(data = ndviSelect,aes(x = GRYLD.y), colour = "blue") +
  geom_vline(xintercept = mean(dat18$GRYLD), linetype = 2) +
  geom_vline(xintercept = mean(ndviSelect$GRYLD.y), linetype = 2, 
             colour = "blue") +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2, 
             colour = "#737373") +
  theme_bw() +
  labs(title = "Distribution of GRYLD in 2017/2018", 
       subtitle = "All lines - black, lines selected from NDVI_20170512 - blue, mean from 2016/2017 - grey")
mean(ndviSelect$GRYLD.x)
mean(ndviSelect$GRYLD.y)
mean(dat17$GRYLD)
mean(dat18$GRYLD)

t.test(dat18$GRYLD,ndviSelect$GRYLD.y)
t.test(dat18$GRYLD,dat17$GRYLD)

# Check RedEdge because it's opposite
dat17 %>% 
  ggplot(aes(x = RedEdge_20170512)) +
  geom_density() +
  geom_vline(xintercept = mean(dat17$RedEdge_20170512), linetype = 2) +
  theme_bw()

RE0512_top<- dat17 %>% 
  dplyr::arrange(RedEdge_20170512) %>% #use desc() if you want highest to lowest
  tidylog::select(rn,GRYLD,RedEdge_20170512) %>% 
  glimpse()
0.05*nrow(dat17)
RE0512_top<- RE0512_top[1:round(0.05*nrow(dat17)),]

dat17 %>% 
  ggplot(aes(x = RedEdge_20170512)) +
  geom_density() +
  geom_vline(xintercept = mean(dat17$RedEdge_20170512), linetype = 2) +
  geom_vline(xintercept = mean(RE0512_top$RedEdge_20170512), colour = "blue",
             linetype = 2) +
  geom_vline(xintercept = max(RE0512_top$RedEdge_20170512), colour = "blue") +
  theme_bw() +
  labs(title = "Distribution of RedEdge_20170512", 
       subtitle = "All lines - black, lowest 5% of RedEdge_20170512 - blue")

dat17 %>% 
  ggplot(aes(x = GRYLD)) +
  geom_density() +
  geom_density(data = RE0512_top,aes(x = GRYLD), colour = "blue") +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2) +
  geom_vline(xintercept = mean(RE0512_top$GRYLD), linetype = 2, 
             colour = "blue") +
  #geom_vline(xintercept = min(RE0512_top$GRYLD), colour = "blue") +
  #geom_vline(xintercept = max(RE0512_top$GRYLD), colour = "blue") +
  theme_bw() +
  labs(title = "Distribution of GRYLD in 2016/2017", 
       subtitle = "All lines - black, lines selected from RedEdge_20170512 - blue")

reSelect<- RE0512_top %>% 
  left_join(dat18, by = "rn") %>% 
  tidylog::select(rn,GRYLD.x,RedEdge_20170512,GRYLD.y) %>% 
  glimpse

dat18 %>% 
  ggplot(aes(x = GRYLD)) +
  geom_density() +
  geom_density(data = reSelect,aes(x = GRYLD.y), colour = "blue") +
  geom_vline(xintercept = mean(dat18$GRYLD), linetype = 2) +
  geom_vline(xintercept = mean(reSelect$GRYLD.y), linetype = 2, 
             colour = "blue") +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2, 
             colour = "#737373") +
  theme_bw() +
  labs(title = "Distribution of GRYLD in 2017/2018", 
       subtitle = "All lines - black, lines selected from RedEdge_20170512 - blue, mean from 2016/2017 - grey")
mean(reSelect$GRYLD.x)
mean(reSelect$GRYLD.y)
t.test(dat18$GRYLD,reSelect$GRYLD.y)
t.test(reSelect$GRYLD.y,ndviSelect$GRYLD.y)
t.test(reSelect$GRYLD.x,ndviSelect$GRYLD.x)

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

# Predict marker effects

##### Determing marker effects for each trait 2017 
gryldME<- mixed.solve(dat17$GRYLD, Z=snpMatrix17)
ndre14MayME<- mixed.solve(dat17$NDVI_20170512, Z=snpMatrix17)
tidy(cor.test(gryldME$u,ndre14MayME$u))
traitME_17<- as.data.frame(gryldME$u)

colnames(traitME_17)<- "GRYLD"

traits<- dat17 %>% 
  tidylog::select(-rn,-GRYLD) %>% 
  colnames()

for (i in traits) {
  
  print(paste("Working on trait", i))
  y=dat17[[i]]
  y
  meRes<- mixed.solve(y=y, Z=snpMatrix17)
  print(meRes$Vu)
  traitME_17[[i]] <- meRes$u
}

##### Determing marker effects for each trait 2018 

gryldME<- mixed.solve(dat18$GRYLD, Z=snpMatrix18)
re20180613ME<- mixed.solve(dat18$RE_20180613, Z=snpMatrix18)
tidy(cor.test(gryldME$u,re20180613ME$u))

traitME_18<- as.data.frame(gryldME$u)
colnames(traitME_18)<- "GRYLD"

traits<- dat18 %>% 
  tidylog::select(-rn,-GRYLD) %>% 
  colnames()

for (i in traits) {
  
  print(paste("Working on trait", i))
  y=dat18[[i]]
  y
  meRes<- mixed.solve(y=y, Z=snpMatrix18)
  print(meRes$Vu)
  traitME_18[[i]] <- meRes$u
}

beep(3)
dev.list()
graphics.off()

##### Looking at the marker effects distributions
# Trying something.... 

ggplot(data = pheno17, aes(x=GRYLD)) +
  geom_density(colour = "red") +
  geom_vline(xintercept = mean(pheno17$GRYLD),linetype = 2, colour = "red") +
  geom_density(data = pheno18, aes(x=GRYLD), colour="blue") +
  geom_vline(xintercept = mean(pheno18$GRYLD),linetype = 2, colour = "blue") +
  theme_bw() +
  labs(title = "Distribution of GRYLD over 2016/2017 and 2017/2018",
       subtitle = "2016/2017-red, 2017/2018-blue")

ggplot(data = traitME_17, aes(x=GRYLD)) +
  geom_density() +
  theme_bw() +
  labs(title = "Marker Effect distribution for GRYLD2016/2017")

ggplot(data = traitME_18, aes(x=GRYLD)) +
  geom_density() +
  theme_bw() +
  labs(title = "Marker Effect distribution for GRYLD 2017/2018")

ggplot(data = traitME_17, aes(x=GRYLD)) +
  geom_density(colour = "red") +
  geom_density(data = traitME_18, aes(x=GRYLD), colour = "blue") +
  geom_vline(xintercept = mean(traitME_17$GRYLD), colour = "red") +
  geom_vline(xintercept = mean(traitME_18$GRYLD), colour = "blue") +
  theme_bw() +
  labs(title = "Marker Effect distribution for GRYLD 2016/2017 and 2017/2018",
       subtitle = "2016/2017 - red, 2017/2018 - blue")

mean(pheno17$GRYLD)
mean(pheno18$GRYLD)
var(pheno17$GRYLD)
var(pheno18$GRYLD)

mean(traitME_17$GRYLD)
mean(traitME_18$GRYLD)
var(traitME_17$GRYLD)
var(traitME_18$GRYLD)

t.test(traitME_17$GRYLD,traitME_18$GRYLD)
t.test(pheno17$GRYLD,pheno18$GRYLD)
t.test(dat17$GRYLD,dat18$GRYLD)

ggplot() +
  geom_point(aes(x=mean(pheno17$GRYLD), y = mean(traitME_17$GRYLD)),
             colour = "red") +
  geom_point(aes(x=mean(pheno18$GRYLD), y = mean(traitME_18$GRYLD)), 
             colour = "blue") +
  theme_bw() +
  labs(x = "mean GRYLD", y = "mean marker effect",
       title = "Comparison of mean GRYLD and mean marker effect",
       subtitle = "2016/2017 - red, 2017/2018 - blue")

ggplot() +
  geom_point(aes(x=mean(pheno17$GRYLD), y = var(traitME_17$GRYLD)),
             colour = "red") +
  geom_point(aes(x=mean(pheno18$GRYLD), y = var(traitME_18$GRYLD)), 
             colour = "blue") +
  theme_bw() +
  labs(x = "mean GRYLD", y = "variance of marker effect",
       title = "Comparison of mean GRYLD and variance of marker effect",
       subtitle = "2016/2017 - red, 2017/2018 - blue")

# Correlation between ME for GRYLD in different years
tidy(rcorr(traitME_17$GRYLD,traitME_18$GRYLD))

##### Correlation Matrix examination ####

corrMatrix_17<- rcorr(as.matrix(traitME_17))
corrMatrix_18<- rcorr(as.matrix(traitME_18))

##### Distribution of Correlation Matrix 

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

correlationME17<- flattenCorrMatrix(corrMatrix_17$r,
                                    corrMatrix_17$P)
correlationME18<- flattenCorrMatrix(corrMatrix_18$r,
                                    corrMatrix_18$P)

phenoCorrMatrix_17<- dat17 %>% 
  tidylog::select(-rn)
phenoCorrMatrix_17<- rcorr(as.matrix(phenoCorrMatrix_17))

correlationsPheno17<- flattenCorrMatrix(phenoCorrMatrix_17$r,
                                        phenoCorrMatrix_17$P)

phenoCorrMatrix_18<- dat18 %>% 
  tidylog::select(-rn)
phenoCorrMatrix_18<- rcorr(as.matrix(phenoCorrMatrix_18))

correlationsPheno18<- flattenCorrMatrix(phenoCorrMatrix_18$r,
                                        phenoCorrMatrix_18$P)

## Distribution of all correlations for VI and GRYLD ME
correlationME17 %>% 
  ggplot(aes(x=cor)) +
  geom_histogram(binwidth = 0.05,
                 fill = "white",colour="black") +
  theme_bw() +
  labs(title = 
         "Distribution of the correlations between marker effects",
       subtitle = 
         "All VI and GRYLD, generated by rrBLUP 2016/2017 Season") +
  xlim(-1,1)

correlationME18 %>% 
  ggplot(aes(x=cor)) +
  geom_histogram(binwidth = 0.05,
                 fill = "white",colour="black") +
  theme_bw() +
  labs(title = 
         "Distribution of the correlations between marker effects",
       subtitle = 
         "All VI and GRYLD, generated by rrBLUP 2017/2018 Season") +
  xlim(-1,1)

## Distribution for only those related to GRYLD
correlationsGryld17<- correlationsPheno17 %>% 
  tidylog::filter(column == "GRYLD") %>% 
  tidylog::select(-column)
unique(correlationsGryld17$row)

correlationsGryld18<- correlationsPheno18 %>% 
  tidylog::filter(column == "GRYLD") %>% 
  tidylog::select(-column)
unique(correlationsGryld18$row)

correlationME17 %>% 
  tidylog::filter(row == "GRYLD") %>% 
  left_join(correlationsGryld17, by = c("column"="row")) %>% 
  dplyr::rename(CorToGeno=cor.x,Pvalue=p.x,corToPheno=cor.y) %>% 
  tidylog::select(-p.y) %>% 
  glimpse() %>% 
  ggplot(aes(x=CorToGeno)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill="white") +
  theme_bw() +
  labs(
    title = 
      "Distribution of the correlations between marker effects",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2016/2017 Season") +
  xlim(-1,1)

correlationME18 %>% 
  tidylog::filter(row == "GRYLD") %>% 
  left_join(correlationsGryld18, by = c("column"="row")) %>% 
  dplyr::rename(CorToGeno=cor.x,Pvalue=p.x,corToPheno=cor.y) %>% 
  tidylog::select(-p.y) %>% 
  glimpse() %>% 
  ggplot(aes(x=CorToGeno)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill="white") +
  theme_bw() +
  labs(
    title = 
      "Distribution of the correlations between marker effects",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2017/2018 Season") +
  xlim(-1,1)

correlationME17 %>% 
  tidylog::filter(row == "GRYLD") %>% 
  left_join(correlationsGryld17, by = c("column"="row")) %>% 
  dplyr::rename(CorToGeno=cor.x,Pvalue=p.x,corToPheno=cor.y) %>% 
  tidylog::select(-p.y) %>% 
  separate(column,c("Trait","Date"), sep = "_") %>% 
  glimpse() %>% 
  ggplot(aes(x=CorToGeno, y=corToPheno,colour=Date,shape=Trait)) +
  geom_point(alpha=0.75,size=3) +
  scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a",
                                "#66a61e","#e6ab02","#a6761d","#666666")) +
  theme_bw()  +
  xlim(-1,1) +
  ylim(-1,1) +
  labs(
    title = 
      "Correlation between Marker Effects correlations and Phenotypic measurements correlations",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2016/2017 Season")

correlationME18 %>% 
  tidylog::filter(row == "GRYLD") %>% 
  left_join(correlationsGryld18, by = c("column"="row")) %>% 
  dplyr::rename(CorToGeno=cor.x,Pvalue=p.x,corToPheno=cor.y) %>% 
  tidylog::select(-p.y) %>% 
  separate(column,c("Trait","Date"), sep = "_") %>% 
  glimpse() %>% 
  ggplot(aes(x=CorToGeno, y=corToPheno,colour=Date,shape=Trait)) +
  geom_point(alpha=0.75,size=3) +
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                                "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
                                '#cab2d6',"#6a3d9a","#ffff4d","#b15928",
                                "#9a0036","#555555","#2a2a2a")) +
  theme_bw()  +
  xlim(-1,1) +
  ylim(-1,1) +
  labs(
    title = 
      "Correlation between Marker Effects correlations and Phenotypic measurements correlations",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2017/2018 Season")

##### Converting Correlation matrix to a distance matrix ####

distmat_17<- distanceMatrix(as.matrix(traitME_17),
                            metric = "absolute pearson")
distmat_18<- distanceMatrix(as.matrix(traitME_18),
                            metric = "absolute pearson")

pcoord_17<- pcoa(distmat_17)
pcoord_18<- pcoa(distmat_18)
biplot.pcoa(pcoord_17)
biplot.pcoa(pcoord_18)

## Making better biplot
pcoOrd_17<- setDT(as.data.frame(pcoord_17$vectors),keep.rownames = T)
pcoOrd_18<- setDT(as.data.frame(pcoord_18$vectors),keep.rownames = T)

pcoOrd_17<- pcoOrd_17 %>% 
  separate(rn,c("Trait","date"),sep = "_")
pcoOrd_17[1,2]<- "20170613"
pcoOrd_17 %>% 
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0,1,19,2,5,6,7)) +
  scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3",
                                "#e7298a","#66a61e","#e6ab02","#a6761d",
                                '#666666',"#000000")) +
  theme_bw() +
  labs(
    title = 
      "Principal Coordinate analysis of genetic distance matrix",
    subtitle = "2016/2017 season") +
  theme(axis.text = element_text(size = 10)) +
  coord_fixed(xlim = c(-0.6,0.6),ylim = c(-0.6,0.6))

pcoOrd_18<- pcoOrd_18 %>% 
  separate(rn,c("Trait","date"),sep = "_")
pcoOrd_18[1,2]<- "20180615"
pcoOrd_18 %>% 
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0,1,19,2,5,6,7)) +
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c",
                                "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
                                '#cab2d6',"#6a3d9a","#9999ff","#b15928",
                                "#9a0036","#888888","#2a2a2a","#000000")) +
  theme_bw() +
  labs(title = "Principal Coordinate analysis of genetic distance matrix",
       subtitle = "2017/2018 season") +
  theme(axis.text = element_text(size = 10)) +
  coord_fixed(xlim = c(-0.6,0.6),ylim = c(-0.6,0.6))

group1_18<- pcoOrd_18 %>% 
  filter(Axis.1 > 0.04 & Axis.2 > -0.19 & Axis.2 < 0.25) %>% 
  dplyr::select(Trait,date) %>% 
  unite("trait", c("Trait","date"), sep = "_")

group2_18<- pcoOrd_18 %>% 
  filter(Axis.1 < 0.04 & Axis.2 > 0) %>% 
  dplyr::select(Trait,date) %>% 
  unite("trait", c("Trait","date"), sep = "_")

group3_18<- pcoOrd_18 %>% 
  filter(Axis.1 < 0.1 & Axis.2 < 0) %>% 
  dplyr::select(Trait,date) %>% 
  unite("trait", c("Trait","date"), sep = "_")

group1_18[group1_18 == "GRYLD_20180615"] = "GRYLD"

group1ME_18<- traitME_18 %>% 
  dplyr::select(group1_18$trait) %>% 
  as.matrix()
group2ME_18<- traitME_18 %>% 
  dplyr::select(group2_18$trait) %>% 
  as.matrix()
group3ME_18<- traitME_18 %>% 
  dplyr::select(group3_18$trait) %>% 
  as.matrix()

group1Cor_18<- cor(group1ME_18)
group2Cor_18<- cor(group2ME_18)
group3Cor_18<- cor(group3ME_18)

#mantel.test(group1Cor_18, group2Cor_18, graph = T)

## Hierarchical clustering

hClustering_17 <- hclust(distmat_17, method = 'ward.D2')
plot(hClustering_17, hang = -1, no.margin = T)

hClustering_18 <- hclust(distmat_18, method = 'ward.D2')
plot(hClustering_18, hang = -1, no.margin = T)

colors = c("#762a83",
           "#1b7837",
           "#9970ab",
           "#5aae61",
           "#c2a5cf",
           "#a6dba0",
           "#e7d4e8",
           "#d9f0d3")
clus17 = cutree(hClustering_17, 4)
plot(as.phylo(hClustering_17), tip.color = colors[clus17],
     label.offset = 0.01, cex = 0.7, no.margin = T)

clus18 = cutree(hClustering_18, 5)
plot(as.phylo(hClustering_17), tip.color = colors[clus18],
     label.offset = 0.01, cex = 0.7, no.margin = T)

## Converting hclust to dendograms

hClustDen_17<- as.dendrogram(hClustering_17)

ggdendrogram(hClustering_17, rotate = T)

hClustDen_18<- as.dendrogram(hClustering_18)

ggdendrogram(hClustering_18, rotate = T)

hclustDenData_17<- dendro_data(hClustDen_17)
hclustDenData_18<- dendro_data(hClustDen_18)

p <- ggplot(hclustDenData_17$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = hclustDenData_17$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3) +
  ylim(-0.5, 1) +
  theme_bw()
p

p <- ggplot(hclustDenData_18$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = hclustDenData_18$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 2) +
  ylim(-0.5, 1) +
  theme_bw()
p

## Neighbour-joining
nj_17<- nj(distmat_17)
plot.phylo(nj_17, type = "phylogram", use.edge.length = F,
           node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
           edge.color = "black", edge.width = 1, edge.lty = 1, font = 1,
           cex = 0.75, adj = NULL, srt = 0, no.margin = T,
           root.edge = FALSE, label.offset = 0.25, underscore = FALSE,
           x.lim = NULL, y.lim = NULL, direction = "rightwards",
           lab4ut = NULL, tip.color = "black", plot = TRUE,
           rotate.tree = 0, open.angle = 0, node.depth = 1,
           align.tip.label = T)


nj_18<- nj(distmat_18)
plot.phylo(nj_18, type = "phylogram", use.edge.length = F,
           node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
           edge.color = "black", edge.width = 1, edge.lty = 1, font = 1,
           cex = 0.5, adj = NULL, srt = 0, no.margin = T,
           root.edge = FALSE, label.offset = 0.25, underscore = FALSE,
           x.lim = NULL, y.lim = NULL, direction = "rightwards",
           lab4ut = NULL, tip.color = "black", plot = TRUE,
           rotate.tree = 0, open.angle = 0, node.depth = 1,
           align.tip.label = T)

gplots::heatmap.2(as.matrix(distmat_17),
                  margins =c(8,8),trace = "none",
                  dendrogram = "both",
                  density.info = "density",
                  col = "viridis",
                  main = "Hierarchichal Clustering of additive genetic effects 2016/2017 season")

gplots::heatmap.2(as.matrix(distmat_18),
                  margins =c(8,8),trace = "none",
                  dendrogram = "both",
                  density.info = "density",
                  col = "viridis",
                  main = "Hierarchichal Clustering of additive genetic effects 2017/2018 season")

## hierarchical clustering with significance
# Can take a long time

pvClust_17_ward<- pvclust(as.matrix(dat17[,2:ncol(dat17)]), 
                          method.hclust="ward.D2", 
                          method.dist="abscor", 
                          use.cor="pairwise.complete.obs", 
                          nboot=100000)
beep(4)

par(mar=c(0.5,2,2,0.25))

plot(pvClust_17_ward)
pvrect(pvClust_17_ward,alpha = 0.95)
par(mar=c(4,4,2,0.25))
x <- seplot(pvClust_17_ward, identify=TRUE)
print(pvClust_17_ward, which=x)

pvClust_17_ave<- pvclust(as.matrix(dat17[,2:ncol(dat17)]), 
                         method.hclust="average", 
                         method.dist="abscor", 
                         use.cor="pairwise.complete.obs", 
                         nboot=100000)
beep(5)

par(mar=c(0.5,2,2,0.25))

plot(pvClust_17_ave)
pvrect(pvClust_17_ave,alpha = 0.95)
par(mar=c(4,4,2,0.25))
x <- seplot(pvClust_17_ave, identify=TRUE)
print(pvClust_17_ave, which=x)

pvClust_17_com<- pvclust(as.matrix(dat17[,2:ncol(dat17)]), 
                         method.hclust="complete", 
                         method.dist="abscor", 
                         use.cor="pairwise.complete.obs", 
                         nboot=100000)
beep(6)

par(mar=c(0.5,2,2,0.25))

plot(pvClust_17_com)
pvrect(pvClust_17_com,alpha = 0.95)
par(mar=c(4,4,2,0.25))
x <- seplot(pvClust_17_com, identify=TRUE)
print(pvClust_17_com, which=x)

pvClust_18_ward<- pvclust(as.matrix(dat18[,2:ncol(dat18)]),
                          method.hclust="ward.D2", 
                          method.dist="abscor", 
                          use.cor="pairwise.complete.obs", 
                          nboot=100000)
beep(7)

par(mar=c(0.5,2,2,0.25))

plot(pvClust_18_ward)
pvrect(pvClust_18_ward,alpha = 0.95)
par(mar=c(4,4,2,0.25))
x <- seplot(pvClust_18_ward, identify=TRUE)
print(pvClust_18_ward, which=x)

pvClust_18_ave<- pvclust(as.matrix(dat18[,2:ncol(dat18)]), 
                         method.hclust="average", 
                         method.dist="abscor", 
                         use.cor="pairwise.complete.obs", 
                         nboot=100000)
beep(8)

par(mar=c(0.5,2,2,0.25))

plot(pvClust_18_ave)
pvrect(pvClust_18_ave,alpha = 0.95)
par(mar=c(4,4,2,0.25))
x <- seplot(pvClust_18_ave, identify=TRUE)
print(pvClust_18_ave, which=x)

pvClust_18_com<- pvclust(as.matrix(dat18[,2:ncol(dat18)]), 
                         method.hclust="complete", 
                         method.dist="abscor", 
                         use.cor="pairwise.complete.obs", 
                         nboot=100000)
beep(8)

par(mar=c(0.5,2,2,0.25))

plot(pvClust_18_com)
pvrect(pvClust_18_com,alpha = 0.95)
par(mar=c(4,4,2,0.25))
x <- seplot(pvClust_18_com, identify=TRUE)
print(pvClust_18_com, which=x)

