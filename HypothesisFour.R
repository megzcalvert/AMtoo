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
  filter(ID != "height") %>% 
  tidylog::select(-entity_id, -Variety, -year) %>% 
  group_by(Date, ID) %>% 
  nest() %>% 
  mutate(correlation = map(data, ~ cor.test(.x$GRYLD,.x$value)),
         tidyCor = map(correlation,glance)) %>% 
  unnest(tidyCor) %>% 
  tidylog::select(-data,-correlation)

nested17$Date<- as.Date(nested17$Date)
nested18$Date<- as.Date(nested18$Date)

nested17 %>% 
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) + 
  geom_hline(yintercept = 0, linetype = 2, 
             colour = "darkgrey") +
  theme_bw() +
  scale_x_date(date_breaks = "1 week", 
               date_labels = "%d%b") +
  labs(title = "Correlation with CI 2017") +
  ylab("Pearson correlation co-efficient")

nested18 %>% 
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) + 
  geom_hline(yintercept = 0, linetype = 2, 
             colour = "darkgrey") +
  theme_bw() +
  scale_x_date(date_breaks = "1 week", 
               date_labels = "%d%b") +
  labs(title = "Correlation with CI 2018") +
  ylab("Pearson correlation co-efficient")

write.table(nested17, "./Phenotype_Database/Correlation_VI_2017.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(nested18, "./Phenotype_Database/Correlation_VI_2018.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
