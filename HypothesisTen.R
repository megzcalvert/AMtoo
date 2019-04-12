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

# normalisedPheno17<-mvn(pheno17[,3:51], univariateTest = "SW", desc = FALSE,
#                        bc = TRUE, bcType = "rounded")
# p1<-powerTransform(pheno17$GNDVI_20170331 ~ 1)
# summary(p1)
# coef(p1)
# 
# transformed<- as.data.frame((pheno17$GNDVI_20170331 ^ (p1$lambda) -1)/p1$lambda)
# mvn(transformed)

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
  
  blues<- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat17<- blues %>% 
    inner_join(dat17)
  
}

dat17[1:5,1:15]

beep()

##### Making it a function for all of the VI's 2017

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
  
  blues<- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
  blues$rn<- str_remove(blues$rn,"Variety_")
  colnames(blues)[colnames(blues)=="effect"] <- paste(i)
  dat18<- blues %>% 
    inner_join(dat18)
  
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

# Predict marker effects
gryldME<- mixed.solve(dat17$GRYLD, Z=snpMatrix17)
ndre14MayME<- mixed.solve(dat17$NDVI_20170512, Z=snpMatrix17)
tidy(cor.test(gryldME$u,ndre14MayME$u))

gryldME<- mixed.solve(dat18$GRYLD, Z=snpMatrix18)
re20180613ME<- mixed.solve(dat18$RE_20180613, Z=snpMatrix18)
tidy(cor.test(gryldME$u,re20180613ME$u))

##### Determing marker effects for each trait 2017 

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

##### Determing marker effects for each trait 2017 

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





##### G Matrix examination ####

gMatrix_17<- rcorr(as.matrix(traitME_17))
gMatrix_18<- rcorr(as.matrix(traitME_18))

##### Distribution of G Matrix 

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

correlationME17<- flattenCorrMatrix(gMatrix_17$r,gMatrix_17$P)
correlationME18<- flattenCorrMatrix(gMatrix_18$r,gMatrix_18$P)

phenogMatrix_17<- dat17 %>% 
  tidylog::select(-rn)
phenogMatrix_17<- rcorr(as.matrix(phenogMatrix_17))

correlationsPheno17<- flattenCorrMatrix(phenoCor17$r,phenoCor17$P)

phenoCor18<- dat18 %>% 
  tidylog::select(-rn)
phenoCor18<- rcorr(as.matrix(phenoCor18))

correlationsPheno18<- flattenCorrMatrix(phenoCor18$r,phenoCor18$P)

## Distribution of all correlations for VI and GRYLD ME
correlationME17 %>% 
  ggplot(aes(x=cor)) +
  geom_histogram(binwidth = 0.05,
                 fill = "white",colour="black") +
  theme_bw() +
  labs(title = "Distribution of the correlations between marker effects",
       subtitle = "All VI and GRYLD, generated by rrBLUP 2016/2017 Season")

correlationME18 %>% 
  ggplot(aes(x=cor)) +
  geom_histogram(binwidth = 0.05,
                 fill = "white",colour="black") +
  theme_bw() +
  labs(title = "Distribution of the correlations between marker effects",
       subtitle = "All VI and GRYLD, generated by rrBLUP 2017/2018 Season")

## Distribution for only those related to GRYLD
correlationsGryld17<- correlationsPheno17 %>% 
  tidylog::filter(column == "GRYLD") %>% 
  tidylog::select(-column)
unique(correlationsGryld17$row)
unique(correlationsGryld17$column)

correlationsGryld18<- correlationsPheno18 %>% 
  tidylog::filter(column == "GRYLD") %>% 
  tidylog::select(-column)
unique(correlationsGryld18$row)
unique(correlationsGryld18$column)

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
    title = "Distribution of the correlations between marker effects",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2016/2017 Season")

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
    title = "Distribution of the correlations between marker effects",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2017/2018 Season")

correlationME17 %>% 
  tidylog::filter(row == "GRYLD") %>% 
  left_join(correlationsGryld17, by = c("column"="row")) %>% 
  dplyr::rename(CorToGeno=cor.x,Pvalue=p.x,corToPheno=cor.y) %>% 
  tidylog::select(-p.y) %>% 
  glimpse() %>% 
  ggplot(aes(x=CorToGeno, y=corToPheno,colour=column)) +
  geom_point() +
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
  glimpse() %>% 
  ggplot(aes(x=CorToGeno, y=corToPheno,colour=column)) +
  geom_point() +
  theme_bw()  +
  xlim(-1,1) +
  ylim(-1,1) +
  labs(
    title = 
      "Correlation between Marker Effects correlations and Phenotypic measurements correlations",
    subtitle = 
      "GRYLD correlations with VI, generated by rrBLUP 2017/2018 Season")

##### Converting G matrix to a distance matrix ####

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
pcoOrd_17 %>% 
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0,1,19,2,5,6,7)) +
  theme_bw() +
  labs(title = "Principal Coordinate analysis of genetic distance matrix",
       subtitle = "2016/2017 season") +
  theme(axis.text = element_text(size = 10))

pcoOrd_18<- pcoOrd_18 %>% 
  separate(rn,c("Trait","date"),sep = "_")
pcoOrd_18 %>% 
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0,1,19,2,5,6,7)) +
  theme_bw() +
  labs(title = "Principal Coordinate analysis of genetic distance matrix",
       subtitle = "2017/2018 season") +
  theme(axis.text = element_text(size = 10))

## Hierarchical clustering
hClustering_17 <- hclust(distmat_17, method = 'complete')
plot(hClustering_17, hang = -1)

hClustering_18 <- hclust(distmat_18, method = 'complete')
plot(hClustering_18, hang = -1)

## Converting hclust to dendograms
hClustDen_17<- as.dendrogram(hClustering_17)

ape::plot.phylo(as.phylo(hClustering_17),
                label.offset = 0.01,
                #type = "unrooted",
                cex = 0.7)

ggdendrogram(hClustering_17, rotate = T)

hClustDen_18<- as.dendrogram(hClustering_18)

ape::plot.phylo(as.phylo(hClustering_18),
                label.offset = 0.01,
                #type = "unrooted",
                cex = 0.5)

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
plot.phylo(as.phylo(nj_17),label.offset = 0.01 , cex = 0.5,
           type = "fan")


nj_18<- nj(distmat_18)
plot.phylo(as.phylo(nj_18),label.offset = 0.01 , cex = 0.5,
           type = "fan")


