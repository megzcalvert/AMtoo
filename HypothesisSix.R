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
pheno17<- pheno17 %>% 
  tidylog::select(-GRYLD)
pheno17$Date<- as.Date(pheno17$Date)

colnames(pheno18)
pheno18<- pheno18 %>% 
  tidylog::select(-GRYLD,-year)
pheno18$Date<- as.Date(pheno18$Date)

pheno17 %>% 
  ggplot(aes(x = Date, y = value, colour = Plot_ID)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  theme(legend.position = "none")

pheno18 %>% 
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  theme(legend.position = "none")

pheno18 %>% 
  filter(Date > "2018-04-1") %>% 
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  theme(legend.position = "none")

#check
trial<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  filter(ID == "NDVI") %>% 
  tidylog::select(-ID) 
trial$Date<- as.factor(trial$Date)
fit.lm<- lm(trial$value ~ trial$Date + trial$Variety + trial$Date:trial$Variety)
summary(fit.lm)
tidyFit.lm<- tidy(fit.lm)
glance(fit.lm)
fit.anova<- anova(fit.lm)
tidy(fit.anova)

pheno17$Date<- as.Date(pheno17$Date)
nested17<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  group_by(ID) %>% 
  nest() %>% 
  mutate(regression = 
           map(data, 
               ~lm(.x$value ~ .x$Date + .x$Variety + .x$Date:.x$Variety)),
         Res = map(regression,tidy)) %>% 
  unnest(Res)

pheno18$Date<- as.Date(pheno18$Date)
nested18<- pheno18 %>% 
  tidylog::select(-Plot_ID) %>% 
  group_by(ID) %>% 
  nest() %>% 
  mutate(regression = 
           map(data, 
               ~lm(.x$value ~ .x$Date + .x$Variety + .x$Date:.x$Variety)),
         Res = map(regression,tidy)) %>% 
  unnest(Res)

