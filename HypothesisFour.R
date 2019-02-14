rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(Hmisc)
library(psych)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

pheno17<- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18<- fread("./Phenotype_Database/pheno18_htpLong.txt")

colnames(pheno17)

# Check
trial<- pheno17 %>% 
  filter(ID == "NDVI") %>% 
  filter(Date == "2017-03-31") 
cor.test(trial$value,trial$GRYLD)

nested17<- pheno17 %>% 
  tidylog::select(-Plot_ID, -Variety) %>% 
  group_by(Date, ID) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ cor.test(.x$GRYLD,.x$value)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation)

nested18<- pheno18 %>% 
  tidylog::select(-entity_id, -Variety, -year) %>% 
  group_by(Date, ID) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ cor.test(.x$GRYLD,.x$value)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation)



