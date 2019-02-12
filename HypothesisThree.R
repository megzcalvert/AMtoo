rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
require(lubridate)
library(car)
library(tidylog)
library(broom)
library(readxl)
library(lme4)
library(Hmisc)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

pheno<- fread("./Phenotype_Database/Pheno_Long1718.txt")

phenoGryld<- pheno %>% 
  filter(trait_id == "GRYLD") %>% 
  tidylog::select(entity_id,phenotype_value,Variety,year) %>% 
  rename(GRYLD = phenotype_value)

phenoGryld17<- phenoGryld %>% 
  filter(year == "17") 

phenoGryld18<- phenoGryld %>% 
  filter(year == "18")

##### 2017 HTP VI data load ####
path <- "./Phenotype_Database/2017_Ashland_AM3_traits_UAS/2017_ASH_AM_vis.xlsx"
htp17<- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path) 

htp17long<- map2_df(htp17, names(htp17), ~ mutate(.x, ID = .y)) %>%
  group_by(ID) %>% 
  gather(key = Date,value = value,`20170331`:`20170609`) %>% 
  left_join(phenoGryld17, by = c("Plot_ID" = "entity_id")) %>% 
  select(-year)

htp17long$Date<- as.Date(htp17long$Date,"%Y%m%d")

str(htp17long)

htp17Wide<- htp17long %>% 
  unite("trait_date",c("ID","Date")) %>% 
  spread(trait_date,value) %>% 



##### 2018 HTP VI data load ####

htpPheno<- c("GNDVI","GRVI","height","NDRE","NDVI","Nir","RE")

htpFileLoad<- function(htp, f, ...) {
  for (i in htp) {
    fileNames<- list.files(path = "./Phenotype_Database/2018_Ashland_AM3_traits_UAS",
                           full.names = T,
                           pattern = paste0("_",i))
    
    traitNames<- basename(fileNames) %>%
      str_remove_all(c(".csv"))
    load.file<- function (filename) {
      d<- fread(file = filename,header = TRUE,check.names = F,data.table = F)
      d
    }
    
    data<- lapply(fileNames, load.file)
    names(data)<- traitNames
    data<- plyr::ldply(data, data.frame, .id = "Phenotype")
    print(colnames(data))
    t <- mutate(dcast(data,  
                      Plot_ID ~ Phenotype, 
                      value.var = paste0(colnames(data)[3]), 
                      fun.aggregate = NULL, 
                      na.rm = TRUE))
    head(t)
    f<- left_join(f, t, by = c("entity_id" = "Plot_ID"))
    
  }
  return(f)
}

phenoGryld18<- htpFileLoad(htpPheno,phenoGryld18)

htp18Long<- phenoGryld18 %>% 
  gather(key = Trait, value = value, `20171120_GNDVI`:`20180613_RE`) %>% 
  separate(Trait, c("Date","ID"), sep = "_")
htp18Long$Date<- as.Date(htp18Long$Date,"%Y%m%d")

#### Scatterplot if VI vs GRYLD ####

#### Data correlation with significance ####

phenoMatrix17<- phenoGryld17 %>% 
  tidylog::select(-entity_id,-Variety,-year)

corMat17<- rcorr(as.matrix(phenoMatrix17))

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flatCor17<- flattenCorrMatrix(corMat17$r,corMat17$P)
