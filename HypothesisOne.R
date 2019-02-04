rm(list = objects()); ls()

library(readr)
library(ggbiplot)
library(data.table)
library(tidyverse)
library(RMySQL)
library(janitor)
library(GGally)
require(lubridate)
library(car)

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
  dplyr::filter(Variety != "blank") %>% 
  dplyr::filter(!is.na(Variety)) %>% 
  dplyr::filter(rep != 0)

str(pheno_long)
knitr::kable(tabyl(pheno_long$trait_id), 
             caption = "Final number of observations per traits")
knitr::kable(tabyl(pheno_long$rep), 
             caption = "Final number of observations per rep")
knitr::kable(tabyl(pheno_long$year), 
             caption = "Final number of observations per year")

### Conversion to Wide format

#Add columns for awns and heading date
pheno_long$awns<- ifelse(pheno_long$phenotype_value == "Awned", 
                         pheno_long$phenotype_value,
                         NA)
str(pheno_long)

awns<- pheno_long[ , c(1,11)]

dataMaid::summarize(awns)

awns<- as.data.frame(na.omit(awns))

#Removing confounding phenotypes  
pheno_long$phenotype_date = as.Date(pheno_long$phenotype_date) ## set the date as a Date type
pheno_long<- subset(pheno_long, pheno_long$trait_id != "AWNS" )

pheno_pcthead<- subset(pheno_long, pheno_long$trait_id == "PCTHEAD" )
pheno_long<- subset(pheno_long, pheno_long$trait_id != "PCTHEAD" )

pheno_pcthead$phenotype_date = as.Date(pheno_pcthead$phenotype_date) ## set the date as a Date type
pheno_pcthead = pheno_pcthead[!is.na(pheno_pcthead$phenotype_date), ] ## remove some rows that are missing the phenotype date
pheno_pcthead = pheno_pcthead[order(pheno_pcthead$entity_id, pheno_pcthead$phenotype_date), ]  ## sort the rows by plot and date

runLogReg = function(dates, pctHEAD){
  res = data.frame(phi1=NA, phi2=NA, phi3=NA, hday=NA, hdate=NA)
  
  if(sum(!is.na(pctHEAD))==0){return(list(res=res, dates=dates, pctHEAD = pctHEAD, pred=NA))} ## skip if no phenotype data
  if(max(pctHEAD, na.rm=TRUE)<60){return(list(res=res, dates=dates, pctHEAD = pctHEAD, pred=NA))} ## skip plots that never reach 50% heading
  
  days = yday(dates) ## convert to day of year
  
  ## add some dummy variables for 0% and 100% at 10, 20, 30 days before/after the phenotyping range ##
  days = c(days, min(days)-c(10,20,30), max(days)+c(10,20,30))
  pctH = c(pctHEAD, c(0,0,0,100,100,100))
  
  ##find initial starting values for phi2 and phi3,  fix phi1 at 100
  phi1 = 100
  phi2 = coef(lm(logit(pctH/phi1)~days))[1]
  phi3 = coef(lm(logit(pctH/phi1)~days))[2]
  
  heading_model<-try(nls(pctH~100/(1+exp(-(phi2+phi3*days))), start=list(phi2=phi2,phi3=phi3), trace=FALSE), silent=TRUE)
  
  if(class(heading_model)=="try-error"){return(list(res=res, dates=dates, pctHEAD = pctHEAD, pred=NA))}  ## skip if doesn't converge
  
  #update model coefficients for predictions
  phi2 = coef(heading_model)[1]
  phi3 = coef(heading_model)[2]
  
  ## get predicted value
  pred = getPred(phi1, phi2, phi3)
  
  hday = pred$x[which.min(abs(pred$y-50))] ## get heading day of year
  hdate = as.Date(paste(year(dates[1]),"-01-01", sep="")) + hday ## convert to date
  
  res = data.frame(phi1, phi2, phi3, hday, hdate) ## set values
  
  return(list(res=res, dates=dates, pctHEAD = pctHEAD, pred=pred))
  
}

getPred = function(phi1, phi2, phi3){
  
  x<-c((80*10):(150*10))/10  #construct a range of x values for all days of year
  y<-phi1/(1+exp(-(phi2+phi3*x))) #predicted y values
  pred<-data.frame(x,y) #create the prediction data frame
  
  return(pred)
}

#################### TEST NON-LINEAR FIT #############################

## heading date calculation results
plots = unique(pheno_pcthead$entity_id) ## get list of plots
##plots = plots[grepl("17ASH", plots)] ## work with only 17ASH plots
head.dates = data.frame(plot_id = plots, day=NA, phi1=NA, phi2=NA, phi3=NA)  ## data frame to store results

##first example
test.plot='18ASH30001'

sample=pheno_pcthead[pheno_pcthead$entity_id==test.plot,]
sample$phenotype_value<- as.numeric(sample$phenotype_value)

test = runLogReg(dates=sample$phenotype_date, pctHEAD=sample$phenotype_value)
hddt<- as.data.frame(test$res)

plots = unique(pheno_pcthead$entity_id)

#head.dates = rep(list(list(vis=list(), length=length(plots)))) ## data frame to store results
#names(head.dates) = plots

for (i in plots){
  ##for (i in 681:length(plots)){  ## quick fix for running only 17ASH plots
  ##for (i in 1000:1010){
  ## subsample the data for the plot for analysis
  print(paste("Working on:", i))
  ## run for visual scores
  sample=pheno_pcthead[pheno_pcthead$entity_id == i,]
  sample$phenotype_value<- as.numeric(sample$phenotype_value)
  #sample=pheno_pcthead[pheno_pcthead$entity_id==plots[i],]
  vis.logistic = runLogReg(dates=sample$phenotype_date, pctHEAD=sample$phenotype_value)
  n<- as.data.frame(vis.logistic$res)
  hddt<- rbind(hddt, n)
  #head.dates[[plots[i]]] = list(vis=vis.logistic)
  ##print(head.dates[[plots[i]]]$vis$res)
  
}

hddt<- hddt[2:nrow(hddt),]
hddt<- cbind(plots,hddt)
plotDates<- hddt[,c(1,5)]

pheno_long$phenotype_value<-as.numeric(as.character(pheno_long$phenotype_value))

write.table(pheno_long, 
            file="./Phenotype_Database/Pheno_Long1718.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

pheno_long<- fread("./Phenotype_Database/Pheno_Long1718.txt", header = T)

#Changing from long to wide format

pheno_Wide <- mutate(dcast(pheno_long,  
                           entity_id + year + Variety + rep + block + 
                             range + column ~ trait_id, 
                           value.var = "phenotype_value", 
                           fun.aggregate = median, 
                           na.rm = TRUE),
                     Variety = factor(Variety),
                     rep = factor(rep))

dataMaid::summarize(pheno_Wide[ , c("Variety", "year", "rep")])

pheno_Wide <- left_join(pheno_Wide, awns)
pheno <- left_join(pheno_Wide, plotDates, by = c("entity_id" = "plots"))

head(pheno)

knitr::kable(tabyl(pheno$year), caption = "Number of observations per year")

names(pheno)[names(pheno)=="Variety"] <- "Taxa"

#To be used in the gapit_rrBLUP_AM.R file
write.table(pheno, "./Phenotype_Database/Pheno_1718.txt", sep = "\t") 

rm(awns, data, gndvi, grvi,hHTP,iniNames,modNames,ndre,ndvi,nir,
   pheno_long,pheno_Wide,redE,vis.logistic,test,sample,hddt,head.dates,
   n,pheno_pcthead,plotDates, i, plots, test.plot)


