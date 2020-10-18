rm(list = objects()); ls()

library(MatrixModels)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(asreml)
library(rrBLUP)
library(Hmisc)
library(MVN)

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

snpChip[snpChip == snpChip$allele_a] = 0
snpChip[snpChip == snpChip$allele_b] = 2
snpChip[snpChip == "H"] = 1
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

snpChip[1:10,1:10]

myGD<- setDT(as.data.frame(snpMatrix), keep.rownames = T)
myGD[1:5,1:5]
colnames(myGD)[colnames(myGD)=="rn"] <- "Taxa"
myGM<- snpChip[,1:3]
colnames(myGM)<- c("Name", "Chromosome","Position")

##### Phenotypes ####

pheno17<- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18<- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19<- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoLong<- fread("./Phenotype_Database/Pheno_Long171819.txt")
pheno_awns<- fread("./Phenotype_Database/Pheno_LongAwns.txt")
hddt17<- fread("./Phenotype_Database/HDDT2017.txt")
hddt18<- fread("./Phenotype_Database/HDDT2018.txt")
hddt19<- fread("./Phenotype_Database/HDDT2019.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(pheno19)
glimpse(phenoLong)
glimpse(pheno_awns)

pheno_awns<- pheno_awns %>% 
  tidylog::select(Variety, phenotype_value) %>% 
  tidylog::distinct() %>% 
  arrange(Variety) %>% 
  tidylog::rename(rn = Variety) %>% 
  tidylog::mutate(Awned = if_else(phenotype_value == "Awned", 1, 0)) 
  
phenoLong<- phenoLong %>% 
  dplyr::rename(Plot_ID = entity_id) %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column) %>% 
  arrange(Variety)

pheno17$Date<- as.Date(pheno17$Date,format = "%Y-%m-%d")
pheno17$Date<- format(pheno17$Date, "%Y%m%d")
pheno18$Date<- as.Date(pheno18$Date,format = "%Y-%m-%d")
pheno18$Date<- format(pheno18$Date, "%Y%m%d")
pheno19$Date<- as.Date(pheno19$Date,format = "%Y-%m-%d")
pheno19$Date<- format(pheno19$Date, "%Y%m%d")

pheno17<- pheno17 %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20170331:RedEdge_20170609) %>% 
  tidylog::inner_join(phenoLong) %>% 
  tidylog::inner_join(hddt17, by = c("Plot_ID" = "plots17")) %>% 
  tidylog::select(-hddt17) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse() %>% 
  arrange(Variety)

pheno18<- pheno18 %>% 
  dplyr::rename(Plot_ID = entity_id)  %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20171120:RE_20180613) %>% 
  tidylog::inner_join(phenoLong) %>% 
  tidylog::inner_join(hddt18, by = c("Plot_ID" = "plots18")) %>% 
  tidylog::select(-hddt18) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse() %>% 
  arrange(Variety)

pheno19<- pheno19 %>% 
  dplyr::rename(Plot_ID = entity_id)  %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20190103:RE_20190624) %>% 
  tidylog::inner_join(phenoLong) %>% 
  tidylog::inner_join(hddt19, by = c("Plot_ID" = "plots19")) %>% 
  tidylog::select(-hddt19) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse() %>% 
  arrange(Variety)

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

pheno19$Plot_ID<- as.factor(pheno19$Plot_ID)
pheno19$Variety<- as.factor(pheno19$Variety)
pheno19$block<- as.factor(pheno19$block)
pheno19$rep<- as.factor(pheno19$rep)
pheno19$range<- as.factor(pheno19$range)
pheno19$column<- as.factor(pheno19$column)

##### Normality tests ####
pheno17NormT<- pheno17 %>% 
  tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
norm17<- mvn(pheno17NormT)
norm17$univariateNormality
boxplot.stats(pheno17NormT$GRYLD)$out
boxplot.stats(pheno17NormT$GNDVI_20170331)$out

pheno17clean<- pheno17 %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column) %>% 
  bind_cols(pheno17NormT)

remove_outliers <- function(x, na.rm = T, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  caps <- quantile(x, probs = c(.05, .95), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  #abs(scale(x)) >= 3
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

pheno17clean[,7:ncol(pheno17clean)]<- as.data.frame(
  lapply(pheno17clean[,7:ncol(pheno17clean)],remove_outliers))
pheno17NormTclean<- pheno17clean %>% 
  tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
norm17clean<- mvn(pheno17NormTclean)
norm17clean$univariateNormality
boxplot.stats(pheno17NormTclean$GRYLD)$out
boxplot.stats(pheno17NormTclean$GNDVI_20170331)$out

pheno18NormT<- pheno18 %>% 
  tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
norm18<- mvn(pheno18NormT)
norm18$univariateNormality
boxplot.stats(pheno18NormT$GRYLD)$out
boxplot.stats(pheno18NormT$GNDVI_20171120)$out

pheno18clean<- pheno18 %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column) %>% 
  bind_cols(pheno18NormT)

pheno18clean[,7:ncol(pheno18clean)]<- as.data.frame(
  lapply(pheno18clean[,7:ncol(pheno18clean)],remove_outliers))
pheno18NormTclean<- pheno18clean %>% 
  tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
norm18clean<- mvn(pheno18NormTclean)
norm18clean$univariateNormality
boxplot.stats(pheno18NormTclean$GRYLD)$out
boxplot.stats(pheno18NormTclean$GNDVI_20171120)$out

pheno19NormT<- pheno19 %>% 
  tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
norm19<- mvn(pheno19NormT)
norm19$univariateNormality
boxplot.stats(pheno19NormT$GRYLD)$out
boxplot.stats(pheno19NormT$GNDVI_20190103)$out

pheno19clean<- pheno19 %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column) %>% 
  bind_cols(pheno19NormT)

pheno19clean[,7:ncol(pheno19clean)]<- as.data.frame(
  lapply(pheno19clean[,7:ncol(pheno19clean)],remove_outliers))
pheno19NormTclean<- pheno19clean %>% 
  tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
norm19clean<- mvn(pheno19NormTclean)
norm19clean$univariateNormality
boxplot.stats(pheno19NormTclean$GRYLD)$out
boxplot.stats(pheno19NormTclean$GNDVI_20190103)$out

##### ASREML BLUEs ####

asreml.license.status()

##### Trialing something ####

set.seed(1962)

dat<- pheno17 %>% 
  dplyr::select(GRYLD, Variety, rep, block)

t17<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,
             data = dat)
plot(t17)
blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat17<- blues %>% 
  tidylog::rename(GRYLD = effect) %>% 
  glimpse() 

##### Generating all 2017 VI BLUEs

effectvars <-  c("block", "rep", "Variety", "year", 
                                    "column","range", "Plot_ID","GRYLD")
traits <- pheno17 %>% 
  dplyr::select(!any_of(effectvars)) %>% 
  colnames()
traits
fieldInfo<- pheno17 %>% 
  tidylog::select(Variety, rep, block, column, range)

for (i in traits) {
  print(paste("Working on trait", i))
  j<- pheno17 %>% 
    pull(i)
  
  data<- cbind(fieldInfo, j)
  names(data)<- c("Variety","rep","block","column","range","Trait")
  print(colnames(data))
  
  t17<- asreml(fixed = Trait ~ 0 + Variety,
               random = ~ rep + rep:block,
               data = data)
  
  blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat17<- blues %>% 
    inner_join(dat17)
  
}

dat17[1:5,1:15]

#### All of 2018
dat<- pheno18 %>% 
  dplyr::select(GRYLD,Variety,rep,block) %>% 
  arrange(Variety)

t18<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,
             data = dat)
plot(t18)
blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat18<- blues %>% 
  tidylog::rename(GRYLD = effect) %>% 
  glimpse() 

##### Generating all 2018 VI BLUEs

traits <- pheno18 %>% 
  dplyr::select(!any_of(effectvars)) %>% 
  colnames()
traits
fieldInfo<- pheno18 %>% 
  tidylog::select(Variety, rep, block, column, range)

for (i in traits) {
  print(paste("Working on trait", i))
  j<- pheno18 %>% 
    pull(i)
  data<- cbind(fieldInfo, j)
  names(data)<- c("Variety","rep","block","column","range","Trait")
  print(colnames(data))
  
  t18<- asreml(fixed = Trait ~ 0 + Variety,
               random = ~ rep + rep:block,
               data = data)
  
  blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat18<- blues %>% 
    inner_join(dat18)
  
}

dat18[1:5,1:15]

#### All of 2019
dat<- pheno19 %>% 
  dplyr::select(GRYLD,Variety,rep,block) %>% 
  arrange(Variety)

t19<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,
             data = dat)
plot(t19)
blues<- setDT(as.data.frame(coef(t19)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat19<- blues %>% 
  tidylog::rename(GRYLD = effect) %>% 
  glimpse() 

##### Generating all 2019 VI BLUEs

traits <- pheno19 %>% 
  dplyr::select(!any_of(effectvars)) %>% 
  colnames()
traits

fieldInfo<- pheno19 %>% 
  tidylog::select(Variety, rep, block, column, range)

for (i in traits) {
  print(paste("Working on trait", i))
  j<- pheno19 %>% 
    pull(i)
  
  data<- cbind(fieldInfo, j)

  names(data)<- c("Variety","rep","block","column","range","Trait")
  print(colnames(data))
  
  t19<- asreml(fixed = Trait ~ 0 + Variety,
               random = ~ rep + rep:block,
               data = data)
  
  blues<- setDT(as.data.frame(coef(t19)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat19<- blues %>% 
    inner_join(dat19)
  
}

dat19[1:5,1:15]

snpChip[1:10,1:10]
snpMatrix[1:10,1:10]

snpLines<- setDT(
  as.data.frame(colnames(snpChip[,3:ncol(snpChip)])))
colnames(snpLines)<- "rn"
dat17<- semi_join(dat17,snpLines, by = "rn") %>% 
  dplyr::arrange(var = "rn")
colnames(dat17)[colnames(dat17)=="rn"] <- "Taxa"

dat18<- semi_join(dat18,snpLines, by = "rn") %>% 
 dplyr::arrange(var = "rn")
colnames(dat18)[colnames(dat18)=="rn"] <- "Taxa"

dat19<- semi_join(dat19,snpLines, by = "rn") %>% 
  dplyr::arrange(var = "rn")
colnames(dat19)[colnames(dat19)=="rn"] <- "Taxa"

pheno_awns<- pheno_awns %>% 
  dplyr::select(phenotype_value,rn) %>% 
  distinct() 

#Function to add a column based on a portion of text in another column
ff = function(x, patterns, replacements = patterns, fill = NA, ...)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  
  return(ans)
}

#Adding a year column based on the entity_id
pheno_awns$phenotype_value<- ff(pheno_awns$phenotype_value, 
                                c("Awned","Awnless"), 
                                c("1", "0"),
                     "NA", ignore.case = TRUE)

dat_awns<- semi_join(pheno_awns,snpLines) %>% 
  mutate(phenotype_value = as.numeric(phenotype_value)) %>% 
  tidylog::select(rn, phenotype_value) %>% 
  dplyr::arrange(var = "rn")
colnames(dat_awns)[colnames(dat_awns)=="rn"] <- "Taxa"

write.table(dat17,file = "./Phenotype_Database/ASREMLBlup_2017.txt",quote = F,
            sep = "\t",row.names = F,col.names = T)
write.table(dat18,file = "./Phenotype_Database/ASREMLBlup_2018.txt",quote = F,
            sep = "\t",row.names = F,col.names = T)
write.table(dat19,file = "./Phenotype_Database/ASREMLBlup_2019.txt",quote = F,
            sep = "\t",row.names = F,col.names = T)
write.table(dat_awns,file = "./Phenotype_Database/ASREMLBlup_awns.txt",quote = F,
            sep = "\t",row.names = F,col.names = T)

numberedCols<- paste(1:(ncol(myGD)-1), sep = ",")

myGM<- myGM %>% 
  tidylog::mutate(V = "V") %>% 
  unite(Name, c("V","Name"),sep ="") %>% 
  glimpse

write.table(myGD,"./R/Gapit/HypothesisEleven/myGD.txt",sep = "\t", quote = F,
            col.names = T, row.names = F)
write.table(myGM,"./R/Gapit/HypothesisEleven/myGM.txt",sep = "\t", quote = F,
            col.names = T, row.names = F)

##### GAPIT Trial ####

library(multtest)
library(readr)
library(data.table)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")

source("http://zzlab.net/GAPIT/emma.txt")

setwd("~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/HypothesisEleven/")

#Step 1: Set working directory and import data
myY <- read.table(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/ASREMLBlup_2017.txt", head = TRUE)
myY[1:10,]

myGD <- read.table("./myGD.txt", head = TRUE)
myGD<- myGD %>% 
  dplyr::arrange(var = "Taxa")
myGD[1:5,1:5]

myGM <- read.table("./myGM.txt", head = TRUE)
myGM[1:5,]

getwd()
setwd(
  "~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/HypothesisEleven/PC4_01_2017/")

# #Step 2: Run GAPIT 
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM ,
  PCA.total = 4,
  cutOff = 0.05
)

## 2018
myY <- read.table(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/ASREMLBlup_2018.txt", head = TRUE)
myY[1:10,]

setwd(
  "~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/HypothesisEleven/PC4_01_2018/")

#Step 2: Run GAPIT 
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM ,
  PCA.total = 4,
  cutOff = 0.05
)

## 2019
myY <- read.table(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/ASREMLBlup_2019.txt", head = TRUE)
myY[1:10,]

setwd(
  "~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/HypothesisEleven/PC4_2019")

#Step 2: Run GAPIT 
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM ,
  PCA.total = 4,
  cutOff = 0.05
)
