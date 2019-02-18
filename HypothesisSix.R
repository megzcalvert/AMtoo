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
fit.NDVI<- lm(trial$value ~ trial$Date * trial$Variety)
summary(fit.NDVI)
tidyFit.NDVI<- tidy(fit.NDVI)
glance(fit.NDVI)
anovaNDVI<- anova(fit.NDVI)
anovaNDVI<- anovaNDVI %>% 
  tidy() %>% 
  add_column(ID = "NDVI")

trial<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  filter(ID == "NDRE") %>% 
  tidylog::select(-ID) 
trial$Date<- as.factor(trial$Date)
fit.NDRE<- lm(trial$value ~ trial$Date * trial$Variety)
summary(fit.NDRE)
tidyFit.NDRE<- tidy(fit.NDRE)
glance(fit.NDRE)
anovaNDRE<- anova(fit.NDRE)
anovaNDRE<- anovaNDRE %>% 
  tidy() %>% 
  add_column(ID = "NDRE")

trial<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  filter(ID == "GNDVI") %>% 
  tidylog::select(-ID) 
trial$Date<- as.factor(trial$Date)
fit.GNDVI<- lm(trial$value ~ trial$Date * trial$Variety)
summary(fit.GNDVI)
tidyFit.GNDVI<- tidy(fit.GNDVI)
glance(fit.GNDVI)
anovaGNDVI<- anova(fit.GNDVI)
anovaGNDVI<- anovaGNDVI %>% 
  tidy() %>% 
  add_column(ID = "GNDVI")

trial<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  filter(ID == "GRVI") %>% 
  tidylog::select(-ID) 
trial$Date<- as.factor(trial$Date)
fit.GRVI<- lm(trial$value ~ trial$Date * trial$Variety)
summary(fit.GRVI)
tidyFit.GRVI<- tidy(fit.GRVI)
glance(fit.GRVI)
anovaGRVI<- anova(fit.GRVI)
anovaGRVI<- anovaGRVI %>% 
  tidy() %>% 
  add_column(ID = "GRVI")

trial<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  filter(ID == "NIR") %>% 
  tidylog::select(-ID) 
trial$Date<- as.factor(trial$Date)
fit.NIR<- lm(trial$value ~ trial$Date * trial$Variety)
summary(fit.NIR)
tidyFit.NIR<- tidy(fit.NIR)
glance(fit.NIR)
anovaNIR<- anova(fit.NIR)
anovaNIR<- anovaNIR %>% 
  tidy() %>% 
  add_column(ID = "NIR")

trial<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  filter(ID == "RedEdge") %>% 
  tidylog::select(-ID) 
trial$Date<- as.factor(trial$Date)
fit.RE<- lm(trial$value ~ trial$Date * trial$Variety)
summary(fit.RE)
tidyFit.RE<- tidy(fit.RE)
glance(fit.RE)
anovaRE<- anova(fit.RE)
anovaRE<- anovaRE %>% 
  tidy() %>% 
  add_column(ID = "RE")



pheno17$Date<- as.Date(pheno17$Date)
nested17<- pheno17 %>% 
  tidylog::select(-Plot_ID) %>% 
  group_by(ID) %>% 
  nest() %>% 
  mutate(regression = 
           map(data, 
               ~lm(.x$value ~ .x$Date * .x$Variety)),
         Res = map(regression,anova)) %>% 
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

