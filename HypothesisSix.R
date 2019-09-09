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
pheno19<- fread("./Phenotype_Database/pheno19_htpLong.txt")

colnames(pheno17)
pheno17<- pheno17 %>% 
  tidylog::select(-GRYLD)
pheno17$Date<- as.Date(pheno17$Date)

colnames(pheno18)
pheno18<- pheno18 %>% 
  tidylog::select(-GRYLD,-year) %>% 
  filter(ID != "height")
pheno18$Date<- as.Date(pheno18$Date)

colnames(pheno19)
pheno19<- pheno19 %>% 
  tidylog::select(-GRYLD,-year) %>% 
  filter(ID != "height")
pheno19$Date<- as.Date(pheno19$Date)

pheno17 %>% 
  ggplot(aes(x = Date, y = value, colour = Plot_ID)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  labs(title = "VI by line over time 2016/2017") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 20),
        legend.position = "none")

pheno18 %>% 
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  labs(title = "VI by line over time 2017/2018") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 20),
        legend.position = "none")

pheno18 %>% 
  filter(Date > "2018-04-1") %>% 
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  labs(title = "VI by line over time 2017/2018") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 20),
        legend.position = "none")

pheno19 %>% 
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  labs(title = "VI by line over time 2018/2019") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 20),
        legend.position = "none")

pheno19 %>% 
  filter(Date > "2019-04-1") %>% 
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap( ~ID, scales = "free",ncol = 2) + 
  theme_bw() +
  labs(title = "VI by line over time 2018/2019") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 20),
        legend.position = "none")

#Anova of linear reg
pheno17$Date<- as.factor(pheno17$Date)
nested17<- pheno17 %>%
  tidylog::select(-Plot_ID) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested17, "./Phenotype_Database/ANOVA_VIbyDateVariety17.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

pheno18$Date<- as.factor(pheno18$Date)
nested18All<- pheno18 %>%
  tidylog::select(-entity_id) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested18All, "./Phenotype_Database/ANOVA_VIbyDateVariety18.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

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

write.table(nested18AfterV, 
            "./Phenotype_Database/ANOVA_VIbyDateVariety18_afterVern.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

pheno19$Date<- as.factor(pheno19$Date)
nested19All<- pheno19 %>%
  tidylog::select(-entity_id) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested19All, "./Phenotype_Database/ANOVA_VIbyDateVariety19.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")

nested19AfterV<- pheno19 %>%
  tidylog::select(-entity_id) %>%
  filter(Date != "2019-01-03") %>% 
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested19AfterV, 
            "./Phenotype_Database/ANOVA_VIbyDateVariety19_afterVern.txt",
            quote = F, row.names = F, col.names = T, sep = "\t")


nested17 %>% 
  tidylog::filter(ID !="RedEdge") %>% 
  tidylog::filter(ID !="NIR") %>% 
  ggplot(aes(x=p.value, fill=ID)) +
  geom_histogram()

