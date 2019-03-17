rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(RMySQL)
library(janitor)
library(GGally)
require(lubridate)
library(car)
library(tidylog)
library(broom)
library(gvlma)
library(MASS)

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

### Conversion to Wide format

#Removing confounding phenotypes  - Jesse's Code
######################################################
## set the date as a Date type
pheno_long$phenotype_date = as.Date(pheno_long$phenotype_date) 
pheno_long<- subset(pheno_long, pheno_long$trait_id != "AWNS" )

pheno_pcthead<- subset(pheno_long, pheno_long$trait_id == "PCTHEAD" )
pheno_long<- subset(pheno_long, pheno_long$trait_id != "PCTHEAD" )

## set the date as a Date type
pheno_pcthead$phenotype_date = 
  as.Date(pheno_pcthead$phenotype_date) 
## remove some rows that are missing the phenotype date
pheno_pcthead = 
  pheno_pcthead[!is.na(pheno_pcthead$phenotype_date), ] 
## sort the rows by plot and date
pheno_pcthead = pheno_pcthead[order(pheno_pcthead$entity_id, 
                                    pheno_pcthead$phenotype_date), ]  

runLogReg = function(dates, pctHEAD){
  res = data.frame(phi1=NA, phi2=NA, phi3=NA, hday=NA, hdate=NA)
  ## skip if no phenotype data
  if(sum(!is.na(pctHEAD))==0){return(list(res=res, dates=dates, 
                                          pctHEAD = pctHEAD, pred=NA))} 
  ## skip plots that never reach 50% heading
  if(max(pctHEAD, na.rm=TRUE)<60){return(list(res=res, dates=dates, 
                                              pctHEAD = pctHEAD, pred=NA))} 
  ## convert to day of year
  days = yday(dates) 
  
  ## add some dummy variables for 0% and 100% at 10, 20, 30 
  # days before/after the phenotyping range ##
  days = c(days, min(days)-c(10,20,30), max(days)+c(10,20,30))
  pctH = c(pctHEAD, c(0,0,0,100,100,100))
  
  ##find initial starting values for phi2 and phi3,  fix phi1 at 100
  phi1 = 100
  phi2 = coef(lm(logit(pctH/phi1)~days))[1]
  phi3 = coef(lm(logit(pctH/phi1)~days))[2]
  
  heading_model<-try(nls(pctH~100/(1+exp(-(phi2+phi3*days))), 
                         start=list(phi2=phi2,phi3=phi3), trace=FALSE), 
                     silent=TRUE)
  ## skip if doesn't converge
  if(class(heading_model)=="try-error"){
    return(list(res=res, dates=dates, pctHEAD = pctHEAD, pred=NA))}  
  
  #update model coefficients for predictions
  phi2 = coef(heading_model)[1]
  phi3 = coef(heading_model)[2]
  
  ## get predicted value
  pred = getPred(phi1, phi2, phi3)
  ## get heading day of year
  hday = pred$x[which.min(abs(pred$y-50))] 
  ## convert to date
  hdate = as.Date(paste(year(dates[1]),"-01-01", sep="")) + hday 
  
  ## set values
  res = data.frame(phi1, phi2, phi3, hday, hdate) 
  
  return(list(res=res, dates=dates, pctHEAD = pctHEAD, pred=pred))
  
}

getPred = function(phi1, phi2, phi3){
  #construct a range of x values for all days of year
  x<-c((80*10):(150*10))/10  
  #predicted y values
  y<-phi1/(1+exp(-(phi2+phi3*x))) 
  #create the prediction data frame
  pred<-data.frame(x,y) 
  
  return(pred)
}

#################### TEST NON-LINEAR FIT #############################

## heading date calculation results
## get list of plots
plots = unique(pheno_pcthead$entity_id) 
## work with only 17ASH plots
plots17 = plots[grepl("17ASH", plots)] 
## data frame to store results
head.dates = data.frame(plot_id = plots17, day=NA, phi1=NA, phi2=NA, phi3=NA) 

##first example
test.plot='17ASH10002'

sample=pheno_pcthead[pheno_pcthead$entity_id==test.plot,]
sample$phenotype_value<- as.numeric(sample$phenotype_value)

test = runLogReg(dates=sample$phenotype_date, pctHEAD=sample$phenotype_value)
hddt17<- as.data.frame(test$res)

## data frame to store results
#head.dates = rep(list(list(vis=list(), length=length(plots)))) 
#names(head.dates) = plots

for (i in plots17){
  ## quick fix for running only 17ASH plots
  ##for (i in 681:length(plots)){  
  ##for (i in 1000:1010){
  ## subsample the data for the plot for analysis
  print(paste("Working on:", i))
  ## run for visual scores
  sample=pheno_pcthead[pheno_pcthead$entity_id == i,]
  sample$phenotype_value<- as.numeric(sample$phenotype_value)
  #sample=pheno_pcthead[pheno_pcthead$entity_id==plots[i],]
  vis.logistic = runLogReg(dates=sample$phenotype_date, 
                           pctHEAD=sample$phenotype_value)
  n<- as.data.frame(vis.logistic$res)
  hddt17<- rbind(hddt17, n)
  #head.dates[[plots[i]]] = list(vis=vis.logistic)
  ##print(head.dates[[plots[i]]]$vis$res)
  
}

hddt17<- hddt17[2:nrow(hddt17),]
hddt17<- cbind(plots17,hddt17)

plotDates17<- hddt17[,c(1,5:6)]
names(plotDates17)<- c("plots17","numHddt17","hddt17")

plots18 = plots[grepl("18ASH", plots)] 
## data frame to store results
head.dates = data.frame(plot_id = plots18, day=NA, phi1=NA, phi2=NA, phi3=NA) 

##first example
test.plot='18ASH30001'

sample=pheno_pcthead[pheno_pcthead$entity_id==test.plot,]
sample$phenotype_value<- as.numeric(sample$phenotype_value)

test = runLogReg(dates=sample$phenotype_date, pctHEAD=sample$phenotype_value)
hddt18<- as.data.frame(test$res)

## data frame to store results
#head.dates = rep(list(list(vis=list(), length=length(plots)))) 
#names(head.dates) = plots

for (i in plots18){
  ## quick fix for running only 17ASH plots
  ##for (i in 681:length(plots)){  
  ##for (i in 1000:1010){
  ## subsample the data for the plot for analysis
  print(paste("Working on:", i))
  ## run for visual scores
  sample=pheno_pcthead[pheno_pcthead$entity_id == i,]
  sample$phenotype_value<- as.numeric(sample$phenotype_value)
  #sample=pheno_pcthead[pheno_pcthead$entity_id==plots[i],]
  vis.logistic = runLogReg(dates=sample$phenotype_date, 
                           pctHEAD=sample$phenotype_value)
  n<- as.data.frame(vis.logistic$res)
  hddt18<- rbind(hddt18, n)
  #head.dates[[plots[i]]] = list(vis=vis.logistic)
  ##print(head.dates[[plots[i]]]$vis$res)
  
}

hddt18<- hddt18[2:nrow(hddt18),]
hddt18<- cbind(plots18,hddt18)

plotDates18<- hddt18[,c(1,5:6)]
names(plotDates18)<- c("plots18","numHddt18","hddt18")
#####################################################################

pheno_long$phenotype_value<-as.numeric(as.character(pheno_long$phenotype_value))

write.table(plotDates17, "./Phenotype_Database/HDDT2017.txt",
            col.names = T, row.names = F, sep = "\t",quote = F)
write.table(plotDates18, "./Phenotype_Database/HDDT2018.txt",
            col.names = T, row.names = F, sep = "\t",quote = F)
write.table(pheno_long, 
            file="./Phenotype_Database/Pheno_Long1718.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

rm(hddt17,hddt18,head.dates,iniNames,modNames,n,pheno_long,pheno_pcthead,
   sample,test,vis.logistic,i,plots,plots17,plots18,getPred,runLogReg,test.plot)

pheno_long<- fread("./Phenotype_Database/Pheno_Long1718.txt", header = T)

#Selecting pheno of choice

unique(pheno_long$trait_id)

pheno17<- pheno_long %>% 
  filter(year == "17") %>%
  filter(trait_id == "GRYLD") %>% 
  dplyr::arrange(entity_id) %>%
  left_join(plotDates17, 
            by = c("entity_id" = "plots17")) %>% 
  tidylog::select(entity_id,phenotype_value,hddt17,numHddt17)

pheno18<- pheno_long %>%
  filter(year == "18") %>%
  filter(trait_id == "GRYLD") %>%
  dplyr::arrange(entity_id) %>%
  left_join(plotDates18, 
            by = c("entity_id" = "plots18")) %>%
  tidylog::select(entity_id,phenotype_value,hddt18,numHddt18)

str(pheno17)

# LOESS - Locally Weighted Scatterplot Smoothing
# popular for regression through time lines

ggplot(pheno17,aes(x = hddt17, y = phenotype_value)) +
  geom_point() +
  geom_smooth(method = "lm",color = "#af8dc3",alpha = 0.25) +
  geom_smooth(method = "loess",color = "#7fbf7b",alpha = 0.25) +
  scale_x_date(date_breaks = "3 day",date_labels = "%b %d") +
  theme_bw() +
  labs(y = "GRYLD")

ggplot(pheno18,aes(x = hddt18, y = phenotype_value)) +
  geom_point() +
  geom_smooth(method = "lm",color = "#af8dc3",alpha = 0.25) +
  geom_smooth(method = "loess",color = "#7fbf7b",alpha = 0.25) +
  scale_x_date(date_breaks = "3 day",date_labels = "%b %d") +
  theme_bw() +
  labs(y = "GRYLD")

cor.test(pheno17$numHddt17,pheno17$phenotype_value)
cor.test(pheno18$numHddt18,pheno18$phenotype_value)


linReg2017<- lm(phenotype_value ~ numHddt17, data = pheno17)
summary(linReg2017)
linReg2017tidy<- tidy(linReg2017)
linReg2017augmented<- augment(linReg2017)
linReg2017glance<- glance(linReg2017)

outlierTest(linReg2017)
plot(linReg2017)

linReg2018<- lm(phenotype_value ~ numHddt18, data = pheno18)
summary(linReg2018)
linReg2018tidy<- tidy(linReg2018)
linReg2018augmented<-augment(linReg2018)
linReg2018glance<- glance(linReg2018)

outlierTest(linReg2018)
plot(linReg2018)


# Influential Observations
# added variable plots
avPlot(linReg2017,variable = "numHddt17")
# Cook's D plot
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(pheno17)-length(linReg2017$coefficients)-2)) 
plot(linReg2017, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(linReg2017, id.method="identify", 
              main="Influence Plot", 
              sub="Circle size is proportial to Cook's Distance" )

avPlot(linReg2018,variable = "numHddt18")
# Cook's D plot
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(pheno18)-length(linReg2018$coefficients)-2)) 
plot(linReg2018, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(linReg2018, id.method="identify", 
              main="Influence Plot", 
              sub="Circle size is proportial to Cook's Distance" )

# Non-normality

sresid <- studres(linReg2017) 
hist(sresid, freq=FALSE, 
     main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)

sresid <- studres(linReg2018) 
hist(sresid, freq=FALSE, 
     main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)


# Evaluate homoscedasticity
# non-constant error variance test
ncvTest(linReg2017)
# plot studentized residuals vs. fitted values 
spreadLevelPlot(linReg2017)

ncvTest(linReg2018)
# plot studentized residuals vs. fitted values 
spreadLevelPlot(linReg2018)

# Evaluate Nonlinearity
# component + residual plot 
crPlots(linReg2017)
crPlots(linReg2018)

# Test for Autocorrelated Errors
durbinWatsonTest(linReg2017)
durbinWatsonTest(linReg2018)


gvmodel <- gvlma(linReg2017) 
summary(gvmodel)

gvmodel <- gvlma(linReg2018) 
summary(gvmodel)

