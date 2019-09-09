rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(asreml)

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




