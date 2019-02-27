rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(RMySQL)
library(janitor)
require(lubridate)
library(tidylog)
library(broom)
library(car)
library(readxl)
library(lme4)
library(Hmisc)
library(psych)


# Connect to database
wheatgenetics = dbConnect(MySQL( ),
                          user=rstudioapi::askForPassword("Database user"),
                          dbname='wheatgenetics', host='beocat.cis.ksu.edu',
                          password = 
                            rstudioapi::askForPassword("Database password"),
                          port = 6306) 

#SQL Query to get all AM Panel phenotype data

pheno_query <- "SELECT phenotype.entity_id,
phenotype.trait_id,
phenotype.phenotype_value,
phenotype.phenotype_date,
phenotype.phenotype_person,
plot.plot_name AS 'Variety',
plot.block,
plot.rep,
plot.range,
plot.column
FROM wheatgenetics.phenotype LEFT JOIN wheatgenetics.plot 
ON plot.plot_id = phenotype.entity_id
WHERE wheatgenetics.phenotype.entity_id LIKE '18ASH30%' OR
wheatgenetics.phenotype.entity_id LIKE '17ASH1%';"

#run the query to get plot information

pheno <- dbGetQuery(wheatgenetics, pheno_query)

#save original data
getwd( ) #get working directory set if needed
setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

saveRDS(pheno, "./Phenotype_Database/Pheno1718.RDS") 

dbDisconnect(wheatgenetics) #disconnect from database

#### Data summaries of long format from database

## clean workspace
rm(list = objects()); ls()

getwd()

sessionInfo() # useful infos **reproducible research**

#Read in original data
pheno_long<- readRDS("./Phenotype_Database/Pheno1718.RDS")

str(pheno_long) # structure
head(pheno_long, 10)

iniNames<- tabyl(pheno_long$Variety)
names(iniNames)[1:2] <- c("Variety", "Count")

## Long Format Summaries with modified names
pheno_long$Variety<- tolower(pheno_long$Variety)
pheno_long$Variety<- str_replace_all(pheno_long$Variety, " ", "_")
pheno_long$Variety<- str_replace_all(pheno_long$Variety, "-", "_")
pheno_long$Variety<- str_replace_all(pheno_long$Variety, "'", "")

#Function to add a column based on a portion of text in another column
ff = function(x, patterns, replacements = patterns, fill = NA, ...)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  
  return(ans)
}

#Adding a year column based on the entity_id
pheno_long$year<- ff(pheno_long$entity_id, 
                     c("17ASH","18ASH"), 
                     c("17","18"),
                     "NA", ignore.case = TRUE)

modNames<- tabyl(pheno_long$Variety)
names(modNames)[1:2] <- c("Variety", "Count")
dataMaid::summarize(pheno_long[ , c("Variety", "trait_id", "rep", "year")])
knitr::kable(tabyl(pheno_long$trait_id), 
             caption = "Number of observations per trait")
knitr::kable(tabyl(pheno_long$rep), 
             caption = "Number of observations of reps")
knitr::kable(tabyl(pheno_long$year), 
             caption = "Number of observations per year")

### Removing lines that do not fit the assumptions
pheno_long<- pheno_long %>% 
  filter(Variety != "blank") %>% 
  filter(!is.na(Variety)) %>% 
  filter(rep != 0) %>%
  filter(!is.na(rep))

finalNames<- tabyl(pheno_long$Variety)

str(pheno_long)
knitr::kable(tabyl(pheno_long$trait_id), 
             caption = "Final number of observations per traits")
knitr::kable(tabyl(pheno_long$rep), 
             caption = "Final number of observations per rep")
knitr::kable(tabyl(pheno_long$year), 
             caption = "Final number of observations per year")

### Extract PCTHEAD and add indicator

rm(modNames,finalNames,iniNames)
colnames(pheno_long)

pheno<- pheno_long %>% 
  filter(trait_id == "PCTHEAD")
pheno$phenotype_value<- as.numeric(pheno$phenotype_value)

pheno17<- pheno %>% 
  filter(year == "17") %>% 
  tidylog::select(entity_id, trait_id, phenotype_value, phenotype_date, Variety,
                  block, rep, range, column) %>% 
  arrange(phenotype_date) %>% 
  mutate(Indicator = if_else(phenotype_value >= 50, 
                             .data$phenotype_date, "NA")) %>% 
  glimpse()

pheno17$phenotype_date<- as.Date(pheno17$phenotype_date)

pheno18<- pheno %>% 
  filter(year == "18") %>% 
  tidylog::select(entity_id, trait_id, phenotype_value, phenotype_date, Variety,
                  block, rep, range, column) %>% 
  arrange(phenotype_date) %>%
  mutate(Indicator = if_else(phenotype_value >= 50, 
                             .data$phenotype_date, "NA")) %>%  
  glimpse()

pheno18$phenotype_date<- as.Date(pheno18$phenotype_date)

write.table(pheno17, "./Phenotype_Database/IndicatorPCTHEAD_17.txt", quote = F,
            sep = "\t", na = "NA", row.names = F,col.names = T)
write.table(pheno18, "./Phenotype_Database/IndicatorPCTHEAD_18.txt", quote = F,
            sep = "\t", na = "NA", row.names = F,col.names = T)
pheno17<- fread("./Phenotype_Database/IndicatorPCTHEAD_17.txt")
pheno18<- fread("./Phenotype_Database/IndicatorPCTHEAD_18.txt")

## Read in HTP

##### 2017 HTP VI data load ####
path <- "./Phenotype_Database/2017_Ashland_AM3_traits_UAS/2017_ASH_AM_vis.xlsx"
htp17<- path %>% 
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = path) 

plotData17<- pheno17 %>% 
  tidylog::select(entity_id, Variety, block, rep, range, column) %>% 
  distinct()

htp17long<- map2_df(htp17, names(htp17), ~ mutate(.x, ID = .y)) %>%
  group_by(ID) %>% 
  gather(key = Date,value = value,`20170331`:`20170609`) %>% 
  glimpse()

htp17long$Date<- as.Date(htp17long$Date,"%Y%m%d")
pheno17$phenotype_date<- as.Date(pheno17$phenotype_date)

htp17long<- htp17long %>% 
  rename(entity_id = Plot_ID, trait_id = ID, phenotype_value = value,
         phenotype_date = Date) %>% 
  left_join(plotData17) %>% 
  mutate(Indicator = NA) %>% 
  glimpse() 
  
htpPct17<- htp17long %>% 
  bind_rows(pheno17) %>% 
  arrange(entity_id,phenotype_date,trait_id) %>% 
  tidylog::select(entity_id,Variety,phenotype_date,trait_id,phenotype_value,
                  Indicator, block,rep,range,column) %>% 
  group_by(entity_id) %>% 
  glimpse()

hddt17<- htpPct17 %>% 
  tidylog::select(entity_id,Indicator) %>% 
  tidylog::filter(!is.na(Indicator)) %>% 
  distinct() %>% 
  group_by(entity_id) %>% 
  summarise(first(Indicator)) %>% 
  rename(HDDT = `first(Indicator)`)

htpPct17 <- htpPct17 %>% 
  left_join(hddt17) %>% 
  mutate(Ind = if_else(phenotype_date < HDDT, "1", "2", "NA")) %>% 
  tidylog::select(-Indicator,-HDDT)



##### 2018 HTP VI data load ####

htpPheno<- c("GNDVI","GRVI","height","NDRE","NDVI","Nir","RE")

plotData18<- pheno18 %>% 
  tidylog::select(entity_id, Variety, block, rep, range, column) %>% 
  distinct()

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

htp18Wide<- htpFileLoad(htpPheno,plotData18)

htp18long<- htp18Wide %>% 
  gather(key = "trait_id", value = "phenotype_value",
         `20171120_GNDVI`:`20180613_RE`) %>% 
  separate(trait_id, c("phenotype_date","trait_id"), sep = "_") 

#htp18long$phenotype_date<- as.Date(htp18long$phenotype_date, format = "%Y%m%d")

htpPct18<- htp18long %>% 
  bind_rows(pheno18) %>% 
  arrange(entity_id,phenotype_date,trait_id) %>% 
  tidylog::filter(phenotype_date > "2018-04-01") %>% 
  tidylog::filter(trait_id != "height") %>% 
  tidylog::select(entity_id,Variety,phenotype_date,trait_id,phenotype_value,
                  Indicator, block,rep,range,column) %>% 
  group_by(entity_id) %>% 
  glimpse()

rm(htp17,pheno,pheno_long,plotData17,plotData18,htpPheno,path,htp18Wide)


