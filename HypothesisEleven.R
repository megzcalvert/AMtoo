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
pheno17$Date<- format(pheno17$Date, "%d%b")
pheno18$Date<- as.Date(pheno18$Date,format = "%Y-%m-%d")
pheno18$Date<- format(pheno18$Date, "%d%b")

pheno17<- pheno17 %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_02Jun:RedEdge_31Mar) %>% 
  tidylog::inner_join(phenoLong) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse()

pheno18<- pheno18 %>% 
  dplyr::rename(Plot_ID = entity_id)  %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_04Apr:RE_29May) %>% 
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

##### Making it a function for all of the VI's

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

beep()
dat17[1:5,1:15]
snpChip[1:10,1:10]
snpMatrix[1:10,1:10]
snpLines<- setDT(
  as.data.frame(colnames(snpChip[,3:ncol(snpChip)])))
colnames(snpLines)<- "rn"
dat17<- semi_join(dat17,snpLines, by = "rn")
colnames(dat17)[colnames(dat17)=="rn"] <- "Taxa"

write.table(dat17,file = "./Phenotype_Database/ASREMLBlup_2017.txt",quote = F,
            sep = "\t",row.names = F,col.names = T)

numberedCols<- paste(1:(ncol(myGD)-1), sep = ",")

colnames(myGD[,2:ncol(myGD)])<- numberedCols
myGD[1:5,1:5]
#write.table(myGD,"./R/Gapit/HypothesisEleven/myGD.txt",sep = "\t", quote = F,
#            col.names = T, row.names = F)
#write.table(myGM,"./R/Gapit/HypothesisEleven/myGM.txt",sep = "\t", quote = F,
#            col.names = T, row.names = F)

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
myY[1:10,1:10]

myGD <- read.table("./myGD.txt", head = TRUE)
myGD[1:5,1:5]
myGM <- read.table("./myGM.txt", head = TRUE)
myGM[1:5,]

#Step 2: Run GAPIT 
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM ,
  PCA.total = 10,
  Model.selection = TRUE,
  group.from = 252,
  group.to = 251
)

### How many PCs based on model selection
fileNames<- list.files(path = "./ModelSelection2017",
                       full.names = T,
                       pattern = ".BIC.Model.Selection.Results.csv$")
traitNames<- basename(fileNames) %>%
  str_remove_all("GAPIT.MLM.|.BIC.Model.Selection.Results.csv")

load.file<- function (filename) {
  d<- fread(file = filename,header = TRUE,check.names = F,data.table = F)
  d
}

data17<- fileNames %>% 
  map(load.file) 

names(data17)<- traitNames

bic17<- map2_df(data17, names(data17), ~ mutate(.x, ID = .y)) %>% 
  rename(PC = `Number of PCs/Covariates`, 
         BIC = `BIC (larger is better) - Schwarz 1978`,
         LogLikelihood = `log Likelihood Function Value`) %>% 
  group_by(ID) %>% 
  glimpse()

bic17max<- bic17 %>% 
  tidylog::filter(BIC == max(BIC)) %>% 
  glimpse()

bic17max %>%
  ggplot(aes(x=PC)) +
  geom_histogram() +
  theme_bw()

colnames(bic17)


par(mfcol = c(1,2))


tnoKnoPC<- GWAS(dat17[,c("rn","GRYLD")],geno = snpChip, fixed = NULL, K = NULL, 
                n.PC = 0, min.MAF = 0.05, n.core = 1, P3D = TRUE, plot = TRUE)

K = A.mat(snpMatrix)

tKnoPC<- GWAS(dat17[,c("rn","GRYLD")],geno = snpChip, fixed = NULL, K = K, 
              n.PC = 0, min.MAF = 0.05, n.core = 1, P3D = TRUE, plot = TRUE)
