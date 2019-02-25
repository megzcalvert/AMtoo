rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(asreml)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

pheno17<- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18<- fread("./Phenotype_Database/pheno18_htpLong.txt")
phenoLong<- fread("./Phenotype_Database/Pheno_Long1718.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(phenoLong)

phenoLong<- phenoLong %>% 
  rename(Plot_ID = entity_id) %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column)

pheno17<- pheno17 %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Date,Variety,GRYLD,GNDVI,
                  GRVI,NDRE,NDVI,NIR,RedEdge) %>% 
  tidylog::inner_join(phenoLong) %>% 
  glimpse()

pheno18<- pheno18 %>% 
  spread(key = ID, value = value) %>% 
  rename(Plot_ID = entity_id, NIR = Nir, RedEdge = RE)  %>% 
  tidylog::select(Plot_ID,Date,Variety,GRYLD,GNDVI,
                  GRVI,NDRE,NDVI,NIR,RedEdge) %>% 
  tidylog::inner_join(phenoLong) %>% 
  glimpse()

###############################################################################
####                  ASREML to calculate heritability                  ####


