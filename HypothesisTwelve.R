rm(list = objects()); ls()

library(MatrixModels)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(beepr)
library(rrBLUP)
library(broom)


#### Theme set
custom_theme<- theme_minimal() %+replace%
  theme(axis.title = element_text(colour = "black",
                                  size = rel(2)),
        axis.title.x = element_text(vjust = 0,
                                    margin = margin(t = 0, r = 0.25,
                                                    b = 0, l = 0,
                                                    unit = "cm")),
        axis.title.y = element_text(vjust = 1,
                                    angle = 90,
                                    margin = margin(t = 0, r = 0.25,
                                                    b = 0, l = 0.1,
                                                    unit = "cm")),
        axis.text = element_text(colour = "black",
                                 size = rel(1.5)),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(3,"pt"),
        axis.line = element_line(color = "black",
                                 size = 0.5),
        legend.key.size = unit(4,"lines"),
        # legend.background = element_rect(fill = NULL, colour = NULL),
        # legend.box = NULL,
        legend.margin = margin(t = 0, r = 0.75,
                               b = 0, l = 0.75,
                               unit = "cm"),
        legend.text = element_text(size = rel(2)),
        legend.title = element_text(size = rel(1.5)),
        panel.grid.major = element_line(colour = "#969696",
                                        linetype = 3),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = rel(2),
                                margin = margin(t = 0.1, r = 0.1, 
                                                b = 0.1, l = 0.1,
                                                unit = "cm")),
        plot.margin = margin(t = 0.5, r = 0.5, 
                             b = 0.5, l = 0,
                             unit = "cm"),
        plot.title = element_text(colour = "black",
                                  size = rel(3),
                                  vjust = 0,
                                  hjust = 0,
                                  margin = margin(t = 0.25, r = 0.25, 
                                                  b = 0.5, l = 0.25,
                                                  unit = "cm")),
        strip.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1),
        strip.text = element_text(colour = "black",
                                  size = rel(1)),
        complete = F)

theme_set(custom_theme)

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
pheno19<- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoLong<- fread("./Phenotype_Database/Pheno_Long171819.txt")

dat17<- fread("./Phenotype_Database/Hyp10BLUEs_17.txt")
dat18<- fread("./Phenotype_Database/Hyp10BLUEs_18.txt")
dat19<- fread("./Phenotype_Database/Hyp10BLUEs_19.txt")
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

dat19<- dat19 %>% 
  semi_join(snpMatrix, by = "rn")

snpMatrix19<- snpMatrix %>% 
  semi_join(dat19,by = "rn") 
rownames(snpMatrix19) <- snpMatrix19[,1]
snpMatrix19[,1] <- NULL

snpMatrix19<- as.matrix(snpMatrix19)

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
yield_answer<-mixed.solve(yield, Z = m_train17, #X = covar,
                          SE = TRUE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid17 %*% e
pred_yield = (pred_yield_valid[,1]) + yield_answer$beta
pred_yield
yield_valid = Pheno_valid17[,"GRYLD"]
YLD_accuracy <-cor(pred_yield_valid, yield_valid )
YLD_accuracy

## VI
ndvi0512 = (Pheno_train17[,"NDVI_20170512"])
ndvi0512_answer<- mixed.solve(ndvi0512, Z = m_train17, SE = TRUE)
NDVI0512 = ndvi0512_answer$u
e = as.matrix(NDVI0512)
pred_ndvi0512_valid =  m_valid17 %*% e
pred_ndvi0512 = (pred_ndvi0512_valid[,1]) + ndvi0512_answer$beta
pred_ndvi0512
ndvi0512_valid = Pheno_valid17[,"NDVI_20170512"]
NDVI_20170512_accuracy <-cor(pred_ndvi0512_valid, ndvi0512_valid )
NDVI_20170512_accuracy

re0512 = (Pheno_train17[,"RedEdge_20170512"])
re0512_answer<- mixed.solve(re0512, Z = m_train17, SE = TRUE)
RE0512 = re0512_answer$u
e = as.matrix(RE0512)
pred_re0512_valid =  m_valid17 %*% e
pred_re0512 = (pred_re0512_valid[,1]) + re0512_answer$beta
pred_re0512
re0512_valid = Pheno_valid17[,"RedEdge_20170512"]
RedEdge_20170512_accuracy <-cor(pred_re0512_valid, re0512_valid )
RedEdge_20170512_accuracy

##### Cross-Validation ####

cvWithHistogram<- function(dat,columns_to_remove,cycles,saveFileFigures,
                           selectionTrait, snpMatrix, identifier, ...) {
  effectvars <- names(dat) %in% columns_to_remove
  traits <- colnames(dat[ , !effectvars])
  (traits)
  cycles = cycles
  accuracy = matrix(nrow=cycles, ncol=length(traits))
  colnames(accuracy)<- traits
  
  for(r in 1:cycles) {
    print(paste("Rep cycle: ",r))
    
    Pheno_train<- dat %>%
      dplyr::sample_frac(0.8)
    Pheno_valid<- dat %>%
      anti_join(Pheno_train, by = "rn")
    
    m_train=snpMatrix[Pheno_train$rn,]
    m_valid=snpMatrix[Pheno_valid$rn,]
    
    for (i in traits) {
      print(paste(i))
      trait=(Pheno_train[,paste(i)])
      trait_answer<-mixed.solve(trait, Z=m_train, SE = F)
      TRT = trait_answer$u
      e = as.matrix(TRT)
      pred_trait_valid =  m_valid %*% e
      pred_trait = as.data.frame((pred_trait_valid[,1]) + trait_answer$beta)
      colnames(pred_trait)<- paste(i)
      pred_trait
      trait_valid = Pheno_valid[,paste(i)]
      accuracy[r,paste(i)] <-(cor(pred_trait_valid, trait_valid, 
                                  use="complete"))
      if (i != selectionTrait) {
        
        selectTrait<- dat %>% 
          dplyr::select(rn, paste(selectionTrait))
        
        topSelec<- setDT(pred_trait, keep.rownames = T) 
        topSelec<- topSelec %>% 
          dplyr::top_n(nrow(topSelec) * 0.25) %>% 
          inner_join(selectTrait)
        
        mT<- t.test(topSelec[,paste(selectionTrait)],
                    dat[,paste(selectionTrait)])
        
        plot1<- ggplot(data = pred_trait,
                       mapping = aes_string(i)) +
          geom_histogram(colour = "black",
                         fill = "white",
                         bins = 100) +
          geom_vline(xintercept = mean(pred_trait[,paste(i)]), linetype = 2) +
          geom_histogram(data = topSelec, 
                         mapping = aes_string(paste(i)),
                         colour = "red",
                         bins = 100) +
          geom_vline(xintercept = mean(topSelec[,paste(i)]),
                     colour = "red", linetype = 2) 
        
        plot2<- ggplot(data = dat,
                       mapping = aes_string(selectionTrait)) +
          geom_histogram(colour = "black",
                         fill = "white",
                         bins = 100) + 
          geom_vline(xintercept = mean(dat[,paste(selectionTrait)]), 
                     linetype = 2) +
          geom_histogram(data = topSelec, 
                         mapping = aes_string(selectionTrait),
                         colour = "red",
                         bins = 100) +
          geom_vline(xintercept = mean(topSelec[,paste(selectionTrait)]),
                     colour = "red", linetype = 2) +
          labs(subtitle = paste("T-test p-value = ",mT$p.value))
        
        plots<- ggpubr::ggarrange(plot1,plot2)
        print(plots)
        ggpubr::ggexport(plots, filename = paste0(saveFileFigures,
                                                  identifier,"_",i,"_",r,".png"),
                         width = 1000, height = 550)
      }
      else {
        print("GRYLD")
      }
    }
  }
  return(accuracy)
}

accuracy17<- cvWithHistogram(dat = dat17, columns_to_remove = c("rn"),
                             cycles = 100, 
                             saveFileFigures = "./Figures/Selection/",
                             selectionTrait = "GRYLD",
                             snpMatrix = snpMatrix17,
                             identifier = "Selection17")

write.csv(accuracy17, 
          "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy17.txt", 
          quote = F,row.names = F)
colMeans(accuracy17)

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

accuracy18<- cvWithHistogram(dat = dat18, columns_to_remove = c("rn"),
                             cycles = 100, 
                             saveFileFigures = "./Figures/Selection/",
                             selectionTrait = "GRYLD",
                             snpMatrix = snpMatrix18,
                             identifier = "Selection18")

write.csv(accuracy18, 
          "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy18.txt", 
          quote = F,row.names = F)
colMeans(accuracy18)

# 2019
Pheno_train19<- dat19 %>% 
  dplyr::sample_frac(0.8)
Pheno_valid19<- dat19 %>% 
  anti_join(Pheno_train19, by = "rn")

m_train19=snpMatrix19[Pheno_train19$rn,]
m_valid19=snpMatrix19[Pheno_valid19$rn,]

##### Predicting Phenotypes ####
## GRYLD
yield=(Pheno_train19[,"GRYLD"])
yield_answer<-mixed.solve(yield, Z = m_train19, SE = TRUE)
YLD = yield_answer$u
e = as.matrix(YLD)
pred_yield_valid =  m_valid19 %*% e
pred_yield=(pred_yield_valid[,1]) + yield_answer$beta
pred_yield
yield_valid = Pheno_valid19[,"GRYLD"]
YLD_accuracy <-cor(pred_yield_valid, yield_valid )
YLD_accuracy

##### Cross-Validation ####

accuracy19<- cvWithHistogram(dat = dat19, columns_to_remove = c("rn"),
                             cycles = 100, 
                             saveFileFigures = "./Figures/Selection/",
                             selectionTrait = "GRYLD",
                             snpMatrix = snpMatrix19,
                             identifier = "Selection19")
write.csv(accuracy19, 
          "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy19.txt", 
          quote = F,row.names = F)
colMeans(accuracy19)

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

gs19_100<- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy19.txt",
  header = T)
gs19_100<- gs19_100 %>% 
  mutate(Rep = row.names(gs19_100)) %>% 
  gather(key = "Trait", value = "Accuracy",RE_20190624:GRYLD) %>% 
  # group_by(Trait) %>% 
  # summarise(average = mean(Accuracy), standardDeviation = sd(Accuracy)) %>% 
  separate(Trait, c("Trait","Date"), sep = "_") %>% 
  filter(Trait != "GRYLD")

gs17_100$Date<- as.Date(gs17_100$Date, format = "%Y%m%d")
gs18_100$Date<- as.Date(gs18_100$Date, format = "%Y%m%d")
gs19_100$Date<- as.Date(gs19_100$Date, format = "%Y%m%d")

ggplot(gs19_100, aes(x = Date, y = Accuracy)) + 
  geom_boxplot(aes(group = Date)) +
  facet_wrap(~Trait, ncol = 2, scales = "free") + 
  theme_bw() +
  theme(axis.text = element_text(colour = "black")) +
  coord_cartesian(ylim = c(-1,1)) +
  #scale_x_date(labels = "%m%d") +
  labs(title = "Genomic Prediction Accuracies",
       subtitle = "2018/2019 season",
       y = "Average prediction accuracy")

