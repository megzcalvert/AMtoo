rm(list = objects()); ls()

library(rrBLUP)
library(data.table)
library(tidyverse)
library(readr)
library(janitor)
library(readr)
library(ggplot2)
library(stringr)
library(tidylog)

#### Initial GWAS analysis
getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

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

snpChip[snpChip == snpChip$allele_a] = -1
snpChip[snpChip == snpChip$allele_b] = 1
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

pcaMethods::checkData(snpMatrix)  #Check PCA assumptions

pcaAM<- pcaMethods::pca(snpMatrix, nPcs = 10) #SVD PCA

sumPCA<- as.data.frame(summary(pcaAM))
knitr::kable(sumPCA)

lineInfo <- fread(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/LineDetailsAMPanel.txt", 
  header = T, check.names = F, sep = '\t', data.table = F)
lineInfo$Name <- tolower(lineInfo$Name)
lineInfo$Name<- str_replace_all(lineInfo$Name, " ", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "-", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "'", "")

program<- lineInfo[,c("Name","Program")]

Scores<- as.data.frame(pcaMethods::scores(pcaAM)) 
Scores<- setDT(Scores, keep.rownames = TRUE)
Scores<- left_join(Scores, program, by = c("rn" = "Name"))

pca.plot <- function(x, p, q, results, ...) {
  
  plots<-ggplot(data = x, aes_string(x = p, y = q)) +
    geom_point(position = "jitter",aes(colour = factor(Program))) +
    theme_bw() +
    labs(title = paste0("PCA Plot ", p, " and ", q), 
         x = paste(p, "R2 = ",(round(sumPCA[1,paste0(p)],3))*100,"%"), 
         y = paste(q, "R2 = ",(round(sumPCA[1,paste0(q)],3))*100,"%")) +
    theme(legend.title=element_blank())
  ggsave(paste0("Biplot", p, "and", q,".pdf"), 
         path=paste(results, sep=''))
  print(plots)
  
}

pca.plot(Scores, "PC1", "PC2", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA/')
pca.plot(Scores, "PC1", "PC3", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC2", "PC3", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC1", "PC4", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC2", "PC4", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC3", "PC4", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')


rm(all.opts,blues,chipLines,chrSum,cmg_dot,cor,cordist,cpg_dot,cpm_dot,geno,
   geno.opts,geno.optsDesc,genoOb,her2017,Her2017,lineInfo,misc.opts,
   misc.optsDesc,missingGenotype,missingPhenotype,opts,opts.desc,output,
   pcaAM,pd,pheno,pheno.opts,pheno.optsDesc,phenoDat,phenoLines,phenoOb,
   plotopts,plotopts.optsDesc,program,res2,resultsSingleUnP18,resultsVarSel18,
   Scores,snpMatrix,sumPCA,toplot,towrite,effectvars,hetCodes,missAmbiguous,
   myvars,numPhenos,resDir,traits)


blues<- read.table("./Phenotype_Database/CollaboratorPAconc.txt", 
                   sep = "\t", 
                   header = TRUE, 
                   stringsAsFactors = TRUE)

blues$Taxa<- tolower(blues$Taxa)
blues$Taxa<- str_replace_all(blues$Taxa, " ", "_")
blues$Taxa<- str_replace_all(blues$Taxa, "-", "_")
blues$Taxa<- str_replace_all(blues$Taxa, "'", "")

genoLines<- as.data.frame(colnames(snpChip[,4:ncol(snpChip)]))
names(genoLines)<- "Taxa"

missingGeno<- anti_join(blues,genoLines)

write.table(missingGeno,"~/Desktop/missingGeno_PinAcid.txt", quote = F, 
            sep = "\t", col.names = T,row.names = F)
blues<- blues %>%
  semi_join(genoLines) %>%
  dplyr::select(-std_dev)

gwaBlue<- rrBLUP::GWAS(pheno = blues,
                       geno = snpChip, 
                       fixed = NULL,
                       K = NULL,
                       n.PC = 0,
                       min.MAF = 0.05,
                       P3D = F, 
                       plot = T)

write.table(gwaBlue,"~/Dropbox/Research_Poland_Lab/AM Panel/R/rrBlup/PinAcidgwas.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)


bluesNA<- na.omit(blues)
bluesNA$Taxa<- toupper(bluesNA$Taxa)

write.table(bluesNA, 
            "./Phenotype_Database/cleanPinAcid_Gapit.txt",
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)

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



setwd("~/Dropbox/Research_Poland_Lab/AM Panel/R/Gapit/Collab/")

#Step 1: Set working directory and import data
myY <- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPinAcid_Gapit.txt", 
                  head = TRUE)
myGD <- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/myGD.txt", head = TRUE)
myGM <- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/myGM.txt", head = TRUE)

#Step 2: Run GAPIT 
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM ,
  #KI = FALSE,
  PCA.total = 10,
  Model.selection = TRUE
)


## GS with rrBLUP

# Define training and validation population, use a 60:40 split? 
# 161 lines, 96:65

setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

blues<- read.table("./Phenotype_Database/CollaboratorPAconc.txt", 
                   sep = "\t", 
                   header = TRUE, 
                   stringsAsFactors = TRUE)

blues$Taxa<- tolower(blues$Taxa)
blues$Taxa<- str_replace_all(blues$Taxa, " ", "_")
blues$Taxa<- str_replace_all(blues$Taxa, "-", "_")
blues$Taxa<- str_replace_all(blues$Taxa, "'", "")

genoLines<- as.data.frame(colnames(snpChip[,4:ncol(snpChip)]))
names(genoLines)<- "Taxa"

missingGeno<- anti_join(blues,genoLines)

write.table(missingGeno,"~/Desktop/missingGeno_PinAcid.txt", quote = F, 
            sep = "\t", col.names = T,row.names = F)
blues<- blues %>%
  semi_join(genoLines) %>%
  dplyr::select(-std_dev)

blues<- blues %>% 
  arrange(Taxa)

markerMatrix<- as.data.frame(t(snpChip[,c(-1,-2,-3)]))
markerMatrix<- setDT(markerMatrix, keep.rownames = T)
markerMatrix[1:5,1:5]
markerMatrix<- markerMatrix %>% 
  arrange(rn) %>% 
  semi_join(blues, by = c("rn" = "Taxa"))
markerMatrix[1:5,1:5]

blues<- as.matrix(blues[,2])
markerMatrix<- markerMatrix %>% 
  tidylog::select(-rn)
markerMatrix<- as.matrix(markerMatrix)
markerMatrix[1:5,1:5]

train<- as.matrix(sample(1:161,96))
head(train)
test<- setdiff(1:161,train)
test

pheno_train<- blues[train,]
m_train<- markerMatrix[train,]
str(m_train)
pheno_valid<- blues[test,]
m_valid<- markerMatrix[test,]
str(m_valid)

pinA_answer<- mixed.solve(pheno_train,Z=m_train,K=NULL,SE=F,return.Hinv = F)

PinA<- pinA_answer$u
e<- as.matrix(PinA)
head(e)
pred_pinA_valid<- m_valid %*% e
str(pred_pinA_valid)
pred_pinA<- as.vector(pred_pinA_valid) + pinA_answer$beta
pred_pinA
PinA_accuracy<- cor(pred_pinA_valid,pheno_valid)
PinA_accuracy

#### cross validation for many cycles for yield only
traits=1
cycles=1000
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles) {
  
  train= as.matrix(sample(1:161, 96))
  test<-setdiff(1:161,train)
  Pheno_train=blues[train,]
  m_train=markerMatrix[train,]
  Pheno_valid=blues[test,]
  m_valid=markerMatrix[test,]
  
  pinA=Pheno_train
  pinA_answer<-mixed.solve(pinA, Z=m_train, K=NULL, 
                            SE = FALSE, return.Hinv=FALSE)
  PinA = pinA_answer$u
  e = as.matrix(PinA)
  pred_PinA_valid =  m_valid %*% e
  
  pinA_valid = Pheno_valid
  print(r)
  accuracy[r,1] <-cor(pred_PinA_valid, pinA_valid, use="complete" )
}
mean(accuracy)
sd(accuracy)
sd(accuracy)/sqrt(length(accuracy))
