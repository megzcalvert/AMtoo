rm(list = objects()); ls()

library(readr)
library(ggbiplot)
library(data.table)
library(tidyverse)
library(janitor)
library(GGally)
require(lubridate)
library(car)
library(tidylog)
library(broom)
library(readxl)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

hddt17<- fread("./Phenotype_Database/HDDT2017.txt")
hddt17$hddt17<- as.Date(hddt17$hddt17)
hddt18<- fread("./Phenotype_Database/HDDT2018.txt")
hddt18$hddt18<- as.Date(hddt18$hddt18)

htpPheno<- c("GNDVI","GRVI","height","NDRE","NDVI","Nir","RE")

path <- "./Phenotype_Database/2017_Ashland_AM3_traits_UAS/2017_ASH_AM_vis.xlsx"
htp17<- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path) 

htp17long<- map2_df(htp17, names(htp17), ~ mutate(.x, ID = .y)) %>%
  group_by(ID) %>% 
  gather(key = Date,value = value,`20170331`:`20170609`)
htp17long$Date<- as.Date(htp17long$Date,"%Y%m%d")
str(htp17long)

htp17long %>%
  ggplot(aes(x = Date, y = value, color = Plot_ID)) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~ID, scales = "free") +
  theme(legend.position = "none") + 
  scale_x_date(breaks = "1 week",date_labels = "%b %d")
  
