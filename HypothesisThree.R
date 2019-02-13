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
library(psych)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

pheno<- fread("./Phenotype_Database/Pheno_Long1718.txt")

phenoGryld<- pheno %>% 
  filter(trait_id == "GRYLD") %>% 
  filter(phenotype_value > 0) %>% 
  filter(phenotype_value < 10) %>% 
  tidylog::select(entity_id,phenotype_value,Variety,year) %>% 
  dplyr::rename(GRYLD = phenotype_value)

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
  tidylog::select(-year)

htp17long$Date<- as.Date(htp17long$Date,"%Y%m%d")

str(htp17long)

htp17Wide<- htp17long %>% 
  unite("trait_date",c("ID","Date")) %>% 
  spread(trait_date,value)

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

htp18Wide<- htpFileLoad(htpPheno,phenoGryld18)

htp18Long<- htp18Wide %>% 
  gather(key = Trait, value = value, `20171120_GNDVI`:`20180613_RE`) %>% 
  separate(Trait, c("Date","ID"), sep = "_")
htp18Long$Date<- as.Date(htp18Long$Date,"%Y%m%d")

str(htp18Long)

#### Scatterplot of VI vs GRYLD ####

#GNDVI
htp17long %>% 
  filter(ID == "GNDVI") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "GNDVI", title = "GNDVI vs GRYLD 2017") +
  theme_bw()

htp18Long %>% 
  filter(ID == "GNDVI") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "GNDVI", title = "GNDVI vs GRYLD 2018") +
  theme_bw()

#GRVI
htp17long %>% 
  filter(ID == "GRVI") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "GRVI", title = "GRVI vs GRYLD 2017") +
  theme_bw()

htp18Long %>% 
  filter(ID == "GRVI") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "GRVI", title = "GRVI vs GRYLD 2018") +
  theme_bw()

#NDVI
htp17long %>% 
  filter(ID == "NDVI") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "NDVI", title = "NDVI vs GRYLD 2017") +
  theme_bw()

htp18Long %>% 
  filter(ID == "NDVI") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "NDVI", title = "NDVI vs GRYLD 2018") +
  theme_bw()

#NDRE
htp17long %>% 
  filter(ID == "NDRE") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "NDRE", title = "NDRE vs GRYLD 2017") +
  theme_bw()

htp18Long %>% 
  filter(ID == "NDRE") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "NDRE", title = "NDRE vs GRYLD 2018") +
  theme_bw()

#NIR
htp17long %>% 
  filter(ID == "NIR") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "NIR", title = "NIR vs GRYLD 2017") +
  theme_bw()

htp18Long %>% 
  filter(ID == "Nir") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "NIR", title = "NIR vs GRYLD 2018") +
  theme_bw()

#Rededge
htp17long %>% 
  filter(ID == "RedEdge") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "RedEdge", title = "RedEdge vs GRYLD 2017") +
  theme_bw()

htp18Long %>% 
  filter(ID == "RE") %>% 
  ggplot(aes(x = GRYLD, y = value)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~Date, scales = "free") +
  labs(ylab = "RedEdge", title = "RedEdge vs GRYLD 2018") +
  theme_bw()

#### Data correlation with significance ####

phenoMatrix17<- htp17Wide %>% 
  tidylog::select(-"Plot_ID",-"Variety")

corMat17<- corr.test(as.matrix(phenoMatrix17), method = "pearson",
                     adjust = "holm")

phenoMatrix18<- htp18Wide %>% 
  tidylog::select(-entity_id,-Variety,-year)

corMat18<- corr.test(as.matrix(phenoMatrix18), method = "pearson",
                     adjust = "holm")

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

flatCor17<- flattenCorrMatrix(corMat17$r,corMat17$p)
viGryld17cor<- flatCor17 %>% 
  filter(row == "GRYLD")
Gryld17corCI<- corMat17$ci

flatCor18<- flattenCorrMatrix(corMat18$r,corMat18$p)
viGryld18cor<- flatCor18 %>% 
  filter(row == "GRYLD")
Gryld18corCI<- corMat18$ci

#### Proportion of variance in GRYLD explained by VI by year ####

htp17long<- as_tibble(htp17long)
#check analysis
fit<- lm(htp17Wide$GRYLD ~ htp17Wide$`GNDVI_2017-03-31`, data = htp17Wide)
summary(fit)

reg17GNDVI<- htp17long %>% 
  tidylog::filter(ID == "GNDVI") %>% 
  tidylog::select(-Variety,-Plot_ID,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "GNDVI") %>% 
  unite("VI_date", c("VI","Date"))
  
reg17GRVI<- htp17long %>% 
  tidylog::filter(ID == "GRVI") %>% 
  tidylog::select(-Variety,-Plot_ID,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "GRVI") %>% 
  unite("VI_date", c("VI","Date"))

reg17NDVI<- htp17long %>% 
  tidylog::filter(ID == "NDVI") %>% 
  tidylog::select(-Variety,-Plot_ID,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "NDVI") %>% 
  unite("VI_date", c("VI","Date"))

reg17NDRE<- htp17long %>% 
  tidylog::filter(ID == "NDRE") %>% 
  tidylog::select(-Variety,-Plot_ID,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "NDRE") %>% 
  unite("VI_date", c("VI","Date"))

reg17NIR<- htp17long %>% 
  tidylog::filter(ID == "NIR") %>% 
  tidylog::select(-Variety,-Plot_ID,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "NIR") %>% 
  unite("VI_date", c("VI","Date"))

reg17RE<- htp17long %>% 
  tidylog::filter(ID == "RedEdge") %>% 
  tidylog::select(-Variety,-Plot_ID,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "RE") %>% 
  unite("VI_date", c("VI","Date"))

htp18Long<- as_tibble(htp18Long)
#check analysis
fit<- lm(htp18Wide$GRYLD ~ htp18Wide$`20171120_GNDVI`)
summary(fit)

reg18GNDVI<- htp18Long %>% 
  tidylog::filter(ID == "GNDVI") %>% 
  tidylog::select(-Variety,-entity_id,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "GNDVI") %>% 
  unite("VI_date", c("VI","Date"))

reg18GRVI<- htp18Long %>% 
  tidylog::filter(ID == "GRVI") %>% 
  tidylog::select(-Variety,-entity_id,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "GRVI") %>% 
  unite("VI_date", c("VI","Date"))

reg18NDVI<- htp18Long %>% 
  tidylog::filter(ID == "NDVI") %>% 
  tidylog::select(-Variety,-entity_id,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "NDVI") %>% 
  unite("VI_date", c("VI","Date"))

reg18NDRE<- htp18Long %>% 
  tidylog::filter(ID == "NDRE") %>% 
  tidylog::select(-Variety,-entity_id,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "NDRE") %>% 
  unite("VI_date", c("VI","Date"))

reg18NIR<- htp18Long %>% 
  tidylog::filter(ID == "Nir") %>% 
  tidylog::select(-Variety,-entity_id,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "NIR") %>% 
  unite("VI_date", c("VI","Date"))

reg18RE<- htp18Long %>% 
  tidylog::filter(ID == "RE") %>% 
  tidylog::select(-Variety,-entity_id,-ID) %>% 
  nest(-Date) %>% 
  mutate(test = map(data, ~lm(.x$GRYLD ~ .x$value)),
         tidied = map(test, glance)) %>% 
  unnest(tidied) %>% 
  tidylog::select(-data,-test) %>% 
  add_column(VI = "RE") %>% 
  unite("VI_date", c("VI","Date"))

VarGryldVI<- bind_rows(reg17GNDVI,reg17GRVI,reg17NDRE,reg17NDVI,reg17NIR,
                       reg17RE,reg18GNDVI,reg18GRVI,reg18NDRE,reg18NDVI,
                       reg18NIR,reg17RE) %>% 
  separate(VI_date, c("VI","Date"), sep = "_")
