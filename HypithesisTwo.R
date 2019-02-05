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

#####################################################################
#### 2017 HTP vegetaion indices over time and relation to hddt ####

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
  facet_wrap(~ID, scales = "free", ncol = 2) + 
  scale_x_date(breaks = "1 week",date_labels = "%b %d") +
  theme(legend.position = "none",
        plot.background = element_rect(colour = "white",
                                       fill = "white"),
        panel.grid.major = element_line(color = "#d9d9d9", 
                                        linetype = 2,
                                        size = 0.5),
        panel.grid.minor = element_line(NULL),
        panel.background = element_rect(fill = "white", 
                                        colour = "black"),
        strip.background = element_rect(fill = "white", 
                                        colour = "black")) +
  geom_vline(color = "#41ae76",
             aes(xintercept = hddt17), hddt17)

#####################################################################
#### 2018 HTP vegetaion indices over time and relation to hddt ####

htpFileLoad<- function(htp, ...) {
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
    
  }
  return(t)
}

htp18<- htpFileLoad(htpPheno)

