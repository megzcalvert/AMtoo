rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(asreml)

set.seed(1966)

pheno<- fread("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/heritabilityTest_20190530_PTHT.csv")

pheno<- pheno %>% 
  mutate(range = as.factor(range),
         column = as.factor(column),
         plot_name = as.factor(plot_name),
         rep = as.factor(rep),
         block = as.factor(block),
         ptht_1 = as.numeric(ptht_1))

t19_ptht1<- asreml(fixed = ptht_1 ~ 1,
             random = ~plot_name + rep + rep:block,
             data = pheno)

plot(t19_ptht1)

h2_ptht1<- as.data.frame(summary(t19_ptht1)$varcomp)
print(h2_ptht1)

h2_ptht1<- (h2_ptht1[3,1] / (h2_ptht1[3,1] + (h2_ptht1[4,1]/2)))
h2_ptht1

t19_ptht2<- asreml(fixed = ptht_2 ~ 1,
             random = ~plot_name + rep + rep:block,
             data = pheno)

plot(t19_ptht2)

h_ptht2<- as.data.frame(summary(t19_ptht2)$varcomp)
print(h_ptht2)

h2_ptht2<- (h_ptht2[3,1] / (h_ptht2[3,1] + (h_ptht2[4,1]/2)))
h2_ptht2

#### New data with manual measurements ####
pheno<- fread("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/20190530_RF_AM_PTHT_summary.csv")
glimpse(pheno)
pheno<- pheno %>% 
  rename(Variety = plot_name) %>% 
  mutate(range = as.factor(range),
         column = as.factor(column),
         Variety = as.factor(Variety),
         rep = as.factor(rep),
         block = as.factor(block))

effectvars <- names(pheno) %in% c("block", "rep", "Variety", "column", 
                                    "range", "plot_id")
traits <- colnames(pheno[ , !effectvars])
H2<- data.frame(traits)
H2$Heritability<- NA
fieldInfo<- pheno %>% 
  tidylog::select(Variety, rep, block, column, range)
ntraits<- 1:nrow(H2)

for (i in ntraits) {
  print(paste("Working on trait", H2[i,1]))
  j<- H2[i,1]
  
  data<- cbind(fieldInfo, pheno[,paste(j)])
  names(data)<- c("Variety","rep","block","column","range","Trait")
  
  t17<- asreml(fixed = Trait ~ 1,
               random = ~Variety + rep + rep:block,
               data = data)
  
  (plot(t17))
  dev.off()
  h<- as.data.frame(summary(t17)$varcomp)
  print(paste("Creating Data Frame", j))
  print(h)
  
  h2<- (h[3,1] / (h[3,1] + (h[4,1]/2)))
  h2
  H2[i,2]<- h2
}

H2






