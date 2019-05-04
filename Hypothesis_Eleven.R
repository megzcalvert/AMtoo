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
library(MVN)
library(qvalue)

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
phenoLong<- fread("./Phenotype_Database/Pheno_Long1718.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(phenoLong)

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

##### Normality tests ####
# pheno17NormT<- pheno17 %>% 
#   tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
# norm17<- mvn(pheno17NormT)
# norm17$univariateNormality
# boxplot.stats(pheno17NormT$GRYLD)$out
# boxplot.stats(pheno17NormT$GNDVI_20170331)$out
# 
# pheno17clean<- pheno17 %>% 
#   tidylog::select(Plot_ID,Variety,block,rep,range,column) %>% 
#   bind_cols(pheno17NormT)
# 
# remove_outliers <- function(x, na.rm = T, ...) {
#   qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
#   caps <- quantile(x, probs = c(.05, .95), na.rm = na.rm)
#   H <- 1.5 * IQR(x, na.rm = na.rm)
#   y <- x
#   #abs(scale(x)) >= 3
#   y[x < (qnt[1] - H)] <- NA
#   y[x > (qnt[2] + H)] <- NA
#   y
# }
# 
# pheno17clean[,7:ncol(pheno17clean)]<- as.data.frame(
#   lapply(pheno17clean[,7:ncol(pheno17clean)],remove_outliers))
# pheno17NormTclean<- pheno17clean %>% 
#   tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
# norm17clean<- mvn(pheno17NormTclean)
# norm17clean$univariateNormality
# boxplot.stats(pheno17NormTclean$GRYLD)$out
# boxplot.stats(pheno17NormTclean$GNDVI_20170331)$out
# 
# 
# pheno18NormT<- pheno18 %>% 
#   tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
# norm18<- mvn(pheno18NormT)
# norm18$univariateNormality
# boxplot.stats(pheno18NormT$GRYLD)$out
# boxplot.stats(pheno18NormT$GNDVI_20171120)$out
# 
# pheno18clean<- pheno18 %>% 
#   tidylog::select(Plot_ID,Variety,block,rep,range,column) %>% 
#   bind_cols(pheno18NormT)
# 
# pheno18clean[,7:ncol(pheno18clean)]<- as.data.frame(
#   lapply(pheno18clean[,7:ncol(pheno18clean)],remove_outliers))
# pheno18NormTclean<- pheno18clean %>% 
#   tidylog::select(-Plot_ID,-Variety,-block,-rep,-range,-column)
# norm18clean<- mvn(pheno18NormTclean)
# norm18clean$univariateNormality
# boxplot.stats(pheno18NormTclean$GRYLD)$out
# boxplot.stats(pheno18NormTclean$GNDVI_20171120)$out

##### ASREML BLUEs ####

asreml.license.status()

##### Trialing something ####

set.seed(1962)

t17<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,
             data = pheno17)
plot(t17)
blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat17<- blues %>% 
  rename(GRYLD = effect) %>% 
  glimpse() 

##### Generating all 2017 VI BLUEs

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
  
  blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat17<- blues %>% 
    inner_join(dat17)
  
}

beep(2)

dat17[1:5,1:15]

#### All of 2018
t18<- asreml(fixed = GRYLD ~ 0 + Variety,
             random = ~ rep + rep:block,
             data = pheno18)
plot(t18)
blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
blues$rn<- str_remove(blues$rn,"Variety_")
dat18<- blues %>% 
  rename(GRYLD = effect) %>% 
  glimpse() 

##### Generating all 2018 VI BLUEs

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
  
  blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat18<- blues %>% 
    inner_join(dat18)
  
}

beep(2)

dat18[1:5,1:15]

snpChip[1:10,1:10]
snpMatrix[1:10,1:10]

snpLines<- setDT(
  as.data.frame(colnames(snpChip[,3:ncol(snpChip)])))
colnames(snpLines)<- "rn"
dat17<- semi_join(dat17,snpLines, by = "rn")
colnames(dat17)[colnames(dat17)=="rn"] <- "Taxa"

dat18<- semi_join(dat18,snpLines, by = "rn")
colnames(dat18)[colnames(dat18)=="rn"] <- "Taxa"

write.table(dat17,file = "./Phenotype_Database/ASREMLBlup_2017.txt",quote = F,
            sep = "\t",row.names = F,col.names = T)
write.table(dat18,file = "./Phenotype_Database/ASREMLBlup_2018.txt",quote = F,
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

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("multtest")
# biocLite("snpStats")
# 
# install.packages("gplots")
# install.packages("LDheatmap")
# install.packages("genetics")
# install.packages("ape")
# install.packages("EMMREML")
# install.packages("scatterplot3d")

# library(multtest)
# library(readr)
# library(data.table)
# library(multtest)
# library(gplots)
# library(LDheatmap)
# library(genetics)
# library(ape)
# library(EMMREML)
# library(compiler) #this library is already installed in R
# library("scatterplot3d")
# 
# source("http://zzlab.net/GAPIT/gapit_functions.txt")
# 
# source("http://zzlab.net/GAPIT/emma.txt")
# 
# setwd("~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/HypothesisEleven/")
# 
# #Step 1: Set working directory and import data
# myY <- read.table(
#   "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/ASREMLBlup_2017.txt", head = TRUE)
# myY[1:10,1:10]
# 
# myGD <- read.table("./myGD.txt", head = TRUE)
# myGD[1:5,1:5]
# myGM <- read.table("./myGM.txt", head = TRUE)
# myGM[1:5,]
# 
# setwd(
#   "~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/HypothesisEleven/PC3_05_2017/")
# 
# #Step 2: Run GAPIT 
# myGAPIT <- GAPIT(
#   Y = myY,
#   GD = myGD,
#   GM = myGM ,
#   PCA.total = 3,
#   cutOff = 0.05
# )
# 
# beep(1)

##### rrBLUP ####
library(rrBLUP)
library(tidyverse)

setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

snpPos<- snpChip %>% 
  tidylog::select(rs_number,chrom,pos)
snpMarkers<- snpChip %>% 
  tidylog::select(-rs_number,-chrom,-pos) %>% 
  as.matrix()

snpMarkers[snpMarkers == 0] = -1
snpMarkers[snpMarkers == 1] = 0
snpMarkers[snpMarkers == 2] = 1
colnames(snpMarkers)<- colnames(snpChip[,4:ncol(snpChip)])

snpRR<- cbind(snpPos,snpMarkers)

effectvars <- names(dat17) %in% c("Taxa")
traits <- colnames(dat17[ , !effectvars])
traits

#png(filename = "./R/rrBlup/HypothesisEleven/Plot_2017.png")

for (i in traits) {
  print(paste("working on ", i))
  
  dat<- dat17 %>% 
    tidylog::select("Taxa", paste(i))
  
  png(filename = 
        paste0("./R/rrBlup/HypothesisEleven/PC4_2017/",i,"_2017_4PC.png"))
  
  res<- GWAS(pheno = dat, 
             geno = snpRR, n.PC = 3)
  dev.off()
  
  write.table(res,file = paste0("./R/rrBlup/HypothesisEleven/PC4_2017/",
                                i,"_2017_4PC.txt"), quote = F, sep = "\t", 
              row.names = F,col.names = T)
}

graphics.off()
dev.off()

beep(6)

effectvars <- names(dat18) %in% c("Taxa")
traits <- colnames(dat18[ , !effectvars])
traits


for (i in traits) {
  print(paste("working on ", i))
  
  dat<- dat18 %>% 
    tidylog::select("Taxa", paste(i))
  
  png(filename = 
        paste0("./R/rrBlup/HypothesisEleven/PC4_2018/",i,"_2018_4PC.png"))
  res<- GWAS(pheno = dat, 
             geno = snpRR, n.PC = 3)
  dev.off()
  write.table(res,file = paste0("./R/rrBlup/HypothesisEleven/PC4_2018/",
                                i,"_2018_4PC.txt"), quote = F, sep = "\t", 
              row.names = F,col.names = T)
}

beep(9)

###############################################################################
#####       Comparing the GAPIT and rrBLUP results on the same scale       ####

######### Reading in files as a list of data frames

load.file<- function (filename) {
  d<- fread(file = filename,header = TRUE,check.names = F,data.table = F)
  d
}

## rrBLUP results
fileNames<- list.files(path = "./R/rrBlup/HypothesisEleven/PC3_2017",
                       full.names = T,
                       pattern = "_2017_3PC.txt$")
traitNames<- basename(fileNames) %>%
  str_remove_all("_2017_3PC.txt")

rrB3_17<- lapply(fileNames, load.file)

names(rrB3_17)<- traitNames

fileNames<- list.files(path = "./R/rrBlup/HypothesisEleven/PC4_2017",
                       full.names = T,
                       pattern = "_2017_4PC.txt$")
traitNames<- basename(fileNames) %>%
  str_remove_all("_2017_4PC.txt")

rrB4_17<- lapply(fileNames, load.file)

names(rrB4_17)<- traitNames

fileNames<- list.files(path = "./R/rrBlup/HypothesisEleven/PC3_2018",
                       full.names = T,
                       pattern = "_2018_3PC.txt$")
traitNames<- basename(fileNames) %>%
  str_remove_all("_2018_3PC.txt")

rrB3_18<- lapply(fileNames, load.file)

names(rrB3_18)<- traitNames

fileNames<- list.files(path = "./R/rrBlup/HypothesisEleven/PC4_2018",
                       full.names = T,
                       pattern = "_2018_4PC.txt$")
traitNames<- basename(fileNames) %>%
  str_remove_all("_2018_4PC.txt")

rrB4_18<- lapply(fileNames, load.file)

names(rrB4_18)<- traitNames

snpPos<- rrB3_17$GNDVI_20170331[,1:3]

rrB3_17F<- snpPos %>% 
  bind_cols(map_dfr(rrB3_17,4)) %>% 
  glimpse

rrB3_18F<- snpPos %>% 
  bind_cols(map_dfr(rrB3_18,4)) %>% 
  glimpse

rrB4_17F<- snpPos %>% 
  bind_cols(map_dfr(rrB4_17,4)) %>% 
  glimpse

rrB4_18F<- snpPos %>% 
  bind_cols(map_dfr(rrB4_18,4)) %>% 
  glimpse

rrB3_17F[,4:ncol(rrB3_17F)]<- 1/(10^rrB3_17F[,4:ncol(rrB3_17F)])
rrB3_18F[,4:ncol(rrB3_18F)]<- 1/(10^rrB3_18F[,4:ncol(rrB3_18F)])
rrB4_17F[,4:ncol(rrB4_17F)]<- 1/(10^rrB4_17F[,4:ncol(rrB4_17F)])
rrB4_18F[,4:ncol(rrB4_18F)]<- 1/(10^rrB4_18F[,4:ncol(rrB4_18F)])

## GAPIT results

fileNames<- list.files(path = "./R/Gapit/HypothesisEleven/PC3_01_2017",
                       full.names = T,
                       pattern = "GWAS.Results.csv$")
traitNames<- basename(fileNames) %>%
  str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")

gap3_17<- lapply(fileNames, load.file)

names(gap3_17)<- traitNames

gap3_17F<- gap3_17$GNDVI_20170331[,1:3] %>% 
  glimpse %>% 
  bind_cols(map_dfr(gap3_17,4)) %>% 
  dplyr::rename(chrom = Chromosome, pos = Position) %>% 
  glimpse()

fileNames<- list.files(path = "./R/Gapit/HypothesisEleven/PC4_01_2017",
                       full.names = T,
                       pattern = "GWAS.Results.csv$")
traitNames<- basename(fileNames) %>%
  str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")

gap4_17<- lapply(fileNames, load.file)

names(gap4_17)<- traitNames

gap4_17F<- gap4_17$GNDVI_20170331[,1:3] %>% 
  glimpse %>% 
  bind_cols(map_dfr(gap4_17,4)) %>% 
  dplyr::rename(chrom = Chromosome, pos = Position) %>% 
  glimpse()

fileNames<- list.files(path = "./R/Gapit/HypothesisEleven/PC3_01_2018",
                       full.names = T,
                       pattern = "GWAS.Results.csv$")
traitNames<- basename(fileNames) %>%
  str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")

gap3_18<- lapply(fileNames, load.file)

names(gap3_18)<- traitNames

gap3_18F<- gap3_18$GNDVI_20171120[,1:3] %>% 
  glimpse %>% 
  bind_cols(map_dfr(gap3_18,4)) %>% 
  dplyr::rename(chrom = Chromosome, pos = Position) %>% 
  glimpse()

fileNames<- list.files(path = "./R/Gapit/HypothesisEleven/PC4_01_2018",
                       full.names = T,
                       pattern = "GWAS.Results.csv$")
traitNames<- basename(fileNames) %>%
  str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")

gap4_18<- lapply(fileNames, load.file)

names(gap4_18)<- traitNames

gap4_18F<- gap4_18$GNDVI_20171120[,1:3] %>% 
  glimpse %>% 
  bind_cols(map_dfr(gap4_18,4)) %>% 
  dplyr::rename(chrom = Chromosome, pos = Position) %>% 
  glimpse()

##### Comparison plots between GWAS methods and PCs ####

don_rrB3_17 <- rrB3_17F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB3_17F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

axisdf = don_rrB3_17 %>% 
  group_by(chrom) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don_rrB4_17 <- rrB4_17F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_17F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

don_rrB3_18 <- rrB3_18F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB3_18F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

don_rrB4_18 <- rrB4_18F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_18F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

don_gap3_17 <- gap3_17F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gap3_17F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

don_gap4_17 <- gap4_17F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gap4_17F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

don_gap3_18 <- gap3_18F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gap3_18F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

don_gap4_18 <- gap4_18F %>% 
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(pos)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gap4_18F, ., by=c("chrom"="chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate( BPcum=pos+tot)

ggplot(don_rrB3_17, aes(x=BPcum)) +
  # Show all points
  geom_point(aes(y=-log10(don_rrB3_18$RE_20180606)), 
             colour = "#4d9221",
             alpha=0.5, size=1) +
  geom_point(aes(y = -log10(don_gap3_18$RE_20180606)),
             colour = "#c51b7d",
             alpha = 0.5, size = 1) +
  # Significance Threshold
  geom_hline(yintercept = -log10(0.05/nrow(don_rrB3_17)), linetype = 2) +
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0.05) ) +# remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 16),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16)
  ) +
  labs(title = "GWAS results RE_20180606",
       subtitle = "rrBLUP results - green, GAPIT results - pink, Bonferroni corection alpha = 0.05",
       x = "Chromosome",
       y = "-log10(P)")

ggplot(don_rrB3_18, aes(x=BPcum, colour = as.factor(don_rrB3_18$chrom))) +
  # Show all points
  geom_point(aes(y=-log10(GNDVI_20180613)),
             alpha=0.5, size=1) +
  scale_color_manual(values = rep(c("#2ca25f", "#8856a7","#43a2ca"), 22 )) +
  # Significance Threshold
  geom_hline(yintercept = -log10(0.05/nrow(don_rrB3_18)), linetype = 2) +
  # custom X axis:
  scale_x_continuous( label = c("1A","1B","1D",
                                "2A","2B","2D",
                                "3A","3B","3D",
                                "4A","4B","4D",
                                "5A","5B","5D",
                                "6A","6B","6D",
                                "7A","7B","7D"), 
                      breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0.05) ) +# remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 16),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16)
  ) +
  labs(title = "GWAS results GNDVI_20180613",
       subtitle = "rrBLUP Bonferroni Threshold alpha = 0.05",
       x = "Chromosome",
       y = "-log10(P)") 
