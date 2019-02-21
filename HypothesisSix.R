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

#Anova of linear reg
pheno17$Date<- as.factor(pheno17$Date)
nested17<- pheno17 %>%
  tidylog::select(-Plot_ID) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

pheno18$Date<- as.factor(pheno18$Date)
nested18All<- pheno18 %>%
  tidylog::select(-entity_id) %>%
  filter(ID == "height") %>% 
  group_by(ID) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

nested18AfterV<- pheno18 %>%
  tidylog::select(-entity_id) %>%
  filter(Date != "2017-11-20") %>% 
  filter(Date != "2017-11-27") %>%
  filter(Date != "2017-12-05") %>% 
  filter(Date != "2017-12-15") %>% 
  filter(Date != "2017-12-18") %>% 
  filter(ID != "height") %>% 
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

