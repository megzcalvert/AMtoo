rm(list = objects()); ls()

library(readr)
library(tidyverse)
library(RMySQL)
library(lme4)
library(janitor)
library(readr)
library(data.table)
library(ggplot2)
library(GGally)
library(ggbiplot)
library(qqman)
library(rrBLUP)
require(car)
require(ggplot2)
require(nlme)
require(plyr)
require(lubridate)

# Connect to database
wheatgenetics = dbConnect(MySQL( ),user="megzcalvert",
                          
                          dbname='wheatgenetics', host='beocat.cis.ksu.edu',
                          password = "calvertm" , port = 6306) 

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

FROM wheatgenetics.phenotype LEFT JOIN wheatgenetics.plot ON plot.plot_id = phenotype.entity_id

WHERE wheatgenetics.phenotype.entity_id LIKE '18ASH30%' OR

wheatgenetics.phenotype.entity_id LIKE '17ASH1%';"

#run the query to get plot information

pheno <- dbGetQuery(wheatgenetics, pheno_query)

#save original data
print(getwd( )) #get working directory set if needed
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/R/AMPanel/")

saveRDS(pheno, "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/Pheno1718.RDS") #save original plot information 


dbDisconnect(wheatgenetics) #disconnect from database


#### Data summaries of long format from database

## clean workspace
rm(list = objects()); ls()

getwd()

sessionInfo() # useful infos **reproducible research**

#Read in original data
pheno_long<- readRDS("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/Pheno1718.RDS")

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
                     c("15ASH", "16ASH", "17ASH","18ASH"), 
                     c("15", "16", "17","18"),
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
modNames
modNames <- modNames %>% 
  filter(Variety != "blank")

modNames  

head(pheno_long)
str(pheno_long)
pheno_long<- semi_join(pheno_long, modNames, by = c('Variety' = 'Variety'))
str(pheno_long)
pheno_long<- subset(pheno_long, pheno_long$rep != 0)
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
            file="~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/Pheno_Long1718.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

pheno_long<- fread("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/Pheno_Long1718.txt", header = T)

#Changing from long to wide format

pheno_Wide <- mutate(dcast(pheno_long,  
                           entity_id + year + Variety + rep + block + range + column ~ trait_id, 
                           value.var = "phenotype_value", fun.aggregate = median, 
                           na.rm = TRUE),
                     Variety = factor(Variety),
                     rep = factor(rep))

dataMaid::summarize(pheno_Wide[ , c("Variety", "year", "rep")])

pheno_Wide <- left_join(pheno_Wide, awns)
pheno <- left_join(pheno_Wide, plotDates, by = c("entity_id" = "plots"))

head(pheno)

knitr::kable(tabyl(pheno$year), caption = "Number of observations per year")

names(pheno)[names(pheno)=="Variety"] <- "Taxa"

write.table(pheno, "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/Pheno_1718.txt", sep = "\t") #To be used in the gapit_rrBLUP_AM.R file

rm(awns, data, gndvi, grvi,hHTP,iniNames,modNames,ndre,ndvi,nir,pheno_long,pheno_Wide,redE,vis.logistic,test,sample,hddt,head.dates,n,pheno_pcthead,
   plotDates, i, plots, test.plot)

## Visual Phenotypic Summaries

pheno<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/Pheno_1718.txt", 
                   sep = "\t", header = TRUE, stringsAsFactors = TRUE)

myvars <- names(pheno) %in% c("awns") 
pd<- pheno[!myvars]
plotList<- list(0)

histFacet.plot <- function(x, results, info, ...) {
  
  md <- names(x) %in% c("rn","Taxa","year","rep","block","column",
                        "range", "entity_id")
  traits <- names(x[ , !md])
  plotList = list()
  for (i in traits) {
    thisPlot <- ggplot(data = x, aes_string(x = i)) + 
      geom_histogram(colour="black", fill="white") + 
      facet_grid(x$year ~ .) +
      theme_bw() +
      xlab(paste0(i)) +
      ylab("Frequency") +
      theme(panel.grid.major = element_blank()) +
      theme(panel.grid.minor = element_blank()) +
      theme(axis.text = element_text(size = 15)) +
      theme(axis.title = element_text(size = 15)) +
      theme(strip.text = element_text(size = 15)) 
    
    plotList[[i]] = thisPlot
    print(i)
    
  }
  return(plotList)
}

raw1718<- histFacet.plot(pd,'~/Dropbox/Research_Poland_Lab/AM Panel/Figures/Hist/',
                         "_raw_1718")
ggpubr::ggarrange(plotlist = raw1718, ncol = 2, nrow = 2)

cpg_dot<- ggplot(data = pheno, aes(x = pheno$PTHT, y = pheno$GRWT)) +
  geom_point() + 
  geom_smooth(method = lm) +
  facet_grid(pheno$year ~ .) +
  xlab("Plant Height (cm)") + 
  ylab("Grain weight") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15)) 

cpg_dot

cmg_dot<- ggplot(data = pheno, aes(x = pheno$MOIST, y = pheno$GRWT)) +
  geom_point() + 
  geom_smooth(method = lm) +
  facet_grid(pheno$year ~ .) +
  xlab("Moisture") + 
  ylab("Grain weight") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15)) 

cmg_dot

cpm_dot<- ggplot(data = pheno, aes(x = pheno$PTHT, y = pheno$MOIST)) +
  geom_point() + 
  geom_smooth(method = lm) +
  facet_grid(pheno$year ~ .) +
  xlab("Plant Height (cm)") + 
  ylab("Moisture") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15)) 

cpm_dot

#### Heritability

#Heritability is calcultated with lme4

pheno<- pheno[ ,c(2:ncol(pheno))] 

snpChip <- read_delim("~/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/90KsnpChipHapMap/AMsnpChipImputed.hmp.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

snpChip<- snpChip %>% 
  clean_names()

chipLines<- as.data.frame(unique(colnames(snpChip[, 13:ncol(snpChip)])))
names(chipLines)[1] <- "snpChip"

remove_outliers <- function(x, na.rm = T, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  caps <- quantile(x, probs = c(.05, .95), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  #abs(scale(x)) >= 3
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

pheno<- pheno %>% 
  select(Taxa,year,rep,block,GRWT,MOIST,PTHT,TESTWT)
pheno<- as.data.frame(semi_join(pheno, chipLines, by = c("Taxa" = "snpChip")))
boxplot.stats(pheno[which(pheno$year == "17"),6])$out
boxplot.stats(pheno[which(pheno$year == "17"),5])$out
boxplot.stats(pheno[which(pheno$year == "17"),7])$out
boxplot.stats(pheno[which(pheno$year == "17"),8])$out
boxplot.stats(pheno[which(pheno$year == "18"),6])$out
boxplot.stats(pheno[which(pheno$year == "18"),5])$out
boxplot.stats(pheno[which(pheno$year == "18"),7])$out
boxplot.stats(pheno[which(pheno$year == "18"),8])$out
pheno[which(pheno$year == '17'),c(5:ncol(pheno))]<- as.data.frame(
  lapply(pheno[which(pheno$year == '17'),c(5:ncol(pheno))], remove_outliers))
pheno[which(pheno$year == '18'),c(5:ncol(pheno))]<- as.data.frame(
  lapply(pheno[which(pheno$year == '18'),c(5:ncol(pheno))], remove_outliers))
boxplot.stats(pheno[which(pheno$year == "17"),6])$out
boxplot.stats(pheno[which(pheno$year == "17"),5])$out
boxplot.stats(pheno[which(pheno$year == "17"),7])$out
boxplot.stats(pheno[which(pheno$year == "17"),8])$out
boxplot.stats(pheno[which(pheno$year == "18"),6])$out
boxplot.stats(pheno[which(pheno$year == "18"),5])$out
boxplot.stats(pheno[which(pheno$year == "18"),7])$out
boxplot.stats(pheno[which(pheno$year == "18"),8])$out

histFacet.plot(pheno,'~/Dropbox/Research_Poland_Lab/AM Panel/Figures/Hist/',
               'clean_2018')

cpg_dot<- ggplot(data = pheno, aes(x = pheno$PTHT, y = pheno$GRWT)) +
  geom_point() + 
  geom_smooth(method = lm) +
  facet_grid(pheno$year ~ .) +
  xlab("Plant Height (cm)") + 
  ylab("Grain weight") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15)) 

cpg_dot

cmg_dot<- ggplot(data = pheno, aes(x = pheno$MOIST, y = pheno$GRWT)) +
  geom_point() + 
  geom_smooth(method = lm) +
  facet_grid(pheno$year ~ .) +
  xlab("Moisture") + 
  ylab("Grain weight") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15)) 

cmg_dot

cpm_dot<- ggplot(data = pheno, aes(x = pheno$PTHT, y = pheno$MOIST)) +
  geom_point() + 
  geom_smooth(method = lm) +
  facet_grid(pheno$year ~ .) +
  xlab("Plant Height (cm)") + 
  ylab("Moisture") +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15)) 

cpm_dot

pheno$Taxa <- as.factor(pheno$Taxa)
pheno$rep <- as.factor(pheno$rep)
pheno$block <- as.factor(pheno$block)

calcH2yr <- function(dat, fill = NA, ...) {
  y=length(table(dat$year))
  r=length(table(dat$rep))
  effectvars <- names(dat) %in% c("block", "rep", "Taxa", "year", "column", 
                                  "row", "experiment_id")
  t <- colnames(dat[ , !effectvars])
  for (i in t) {
    
    print(paste("Working on trait", i))
    h = as.data.frame(VarCorr(lmer(paste0(i, "~(1|Taxa) + (1|year) + 
                                          (1| year:rep) + (1|year:rep:block) + 
                                          (1|year:Taxa)"), data = dat)))
    H2= h[2,4] / (h[2,4] + (h[1,4] / y) + (h[6,4] / (y*r)))
    print(H2)
  }
  
}

calcH2yr(pheno)

blues.yrb <- function(traits, dat = ".") {
  b<- as.data.frame(fixef(lmer(paste0(traits, "~ 0 + Taxa + (1|year) + 
(1| year:rep) + (1|year:rep:block) + (1|year:Taxa)"), 
                               data = dat)))
}

effectvars <- names(pheno) %in% c("block", "rep", "Taxa", "year", "column", 
                                   "row", "experiment_id")

traits <- colnames(pheno[ , !effectvars])

bluesA<- lapply(traits, blues.yrb, dat = pheno)

bluesA <- map2(bluesA, traits, ~ set_names(..1, ..2) %>%
                 rownames_to_column(var = "Taxa")) %>%
  reduce(full_join)

bluesA$Taxa<- sub("Taxa","",bluesA$Taxa)

write.table(bluesA, 
            "/Users/megzcalvert/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/cleanPheno1718Blues.txt",
            quote = FALSE, row.names = F, col.names = T, sep = "\t", 
            fileEncoding = "UTF-8")


#### Comparing names to those in genotype file

phenoLines<- as.data.frame(unique(pheno[ , 1]))
names(phenoLines)[1] <- "Phenotype"

missingGenotype <- anti_join(phenoLines, chipLines, by = c("Phenotype" = "snpChip")) 
missingPhenotype <- anti_join(chipLines, phenoLines, by = c("snpChip" = "Phenotype"))
nrow(missingGenotype)
nrow(missingPhenotype)

#### Initial GWAS analysis

snpChip <- read_delim("~/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/90KsnpChipHapMap/AMsnpChipImputed.hmp.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
snpChip<- snpChip %>% 
  clean_names()

missAmbiguous = c('0', '+', '-')
hetCodes = c('R','Y','S','W','K','M','B','D','H','V')
hapgeno=as.matrix(snpChip[,13:ncol(snpChip)])
hapgeno[hapgeno %in% missAmbiguous]=NA
hapgeno[hapgeno=='N']=NA
hapgeno[hapgeno %in% hetCodes]='H'
snpChip=cbind(snpChip[,1:12], hapgeno)
rm(hapgeno)

write.table(snpChip, file="~/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/SelectedImputedBeagle.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="~/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/SelectedImputedBeagle.txt", 
                header=TRUE, check.names=F, sep = "\t")

snpChip[snpChip == snpChip$allele_a] = -1
snpChip[snpChip == snpChip$allele_b] = 1
snpChip[snpChip == "H"] = 0
snpChip[snpChip == "C"] = NA
snpChip[snpChip == "A"] = NA
snpChip[snpChip == "T"] = NA
snpChip[snpChip == "G"] = NA
snpChip[snpChip == "-"] = NA
snpChip[snpChip == "."] = NA

snpChip<- snpChip[ ,c(1,4,5,13:311)]

write.table(snpChip, file="~/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/SelectedImputedBeagleNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="~/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/SelectedImputedBeagleNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")
snpMatrix<- t(snpChip[ , c(-1, -2, -3)])

pcaMethods::checkData(snpMatrix)  #Check PCA assumptions

pcaAM<- pcaMethods::pca(snpMatrix, nPcs = 10) #SVD PCA

sumPCA<- as.data.frame(summary(pcaAM))
knitr::kable(sumPCA)

lineInfo <- fread("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/LineDetailsAMPanel.txt", 
                  header = T, check.names = F, sep = '\t', data.table = F)
lineInfo$Name <- tolower(lineInfo$Name)
lineInfo$Name<- str_replace_all(lineInfo$Name, " ", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "-", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "'", "")

program<- lineInfo[,c("Name","Program")]

Scores<- as.data.frame(pcaMethods::scores(pcaAM)) 
Scores<- setDT(Scores, keep.rownames = TRUE)
Scores<- left_join(Scores, program, by = c("rn" = "Name"))

pca.plot <- function(x, p, q, results, ...) {
  
  plots<-ggplot(data = x, aes_string(x = p, y = q)) +
    geom_point(position = "jitter",aes(colour = factor(Program))) +
    theme_bw() +
    labs(title = paste0("PCA Plot ", p, " and ", q), 
         x = paste(p, "R2 = ",(round(sumPCA[1,paste0(p)],3))*100,"%"), 
         y = paste(q, "R2 = ",(round(sumPCA[1,paste0(q)],3))*100,"%")) +
    theme(legend.title=element_blank())
  ggsave(paste0("Biplot", p, "and", q,".pdf"), 
         path=paste(results, sep=''))
  print(plots)
  
}

pca.plot(Scores, "PC1", "PC2", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA/')
pca.plot(Scores, "PC1", "PC3", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC2", "PC3", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC1", "PC4", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC2", "PC4", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')
pca.plot(Scores, "PC3", "PC4", 
         '~/Dropbox/Research_Poland_Lab/AM Panel/Figures/PCA')


rm(all.opts,blues,chipLines,chrSum,cmg_dot,cor,cordist,cpg_dot,cpm_dot,geno,
   geno.opts,geno.optsDesc,genoOb,her2017,Her2017,lineInfo,misc.opts,
   misc.optsDesc,missingGenotype,missingPhenotype,opts,opts.desc,output,
   pcaAM,pd,pheno,pheno.opts,pheno.optsDesc,phenoDat,phenoLines,phenoOb,
   plotopts,plotopts.optsDesc,program,res2,resultsSingleUnP18,resultsVarSel18,
   Scores,snpMatrix,sumPCA,toplot,towrite,effectvars,hetCodes,missAmbiguous,
   myvars,numPhenos,resDir,traits)


blues<- read.table("/Users/megzcalvert/Dropbox/Research_Poland_Lab/AM Panel/Genotype_Database/cleanPheno1718Blues.txt", 
                   sep = "\t", 
                   header = TRUE, 
                   stringsAsFactors = TRUE)



gwaBlue<- rrBLUP::GWAS(pheno = blues,
                       geno = snpChip, 
                       fixed = NULL,
                       K = NULL,
                       n.PC = 0,
                       min.MAF = 0.05,
                       P3D = F, 
                       plot = F)

blues <- na.omit(blues)

filterChip<- snpChip[,4:ncol(snpChip)]
hist(rowSums(filterChip))
filterChip<- filterChip[which(rowSums(filterChip) != 299 | rowSums(filterChip) != -299), ]

gwaBlueHddt <- rrBLUP::GWAS(pheno = blues, geno = snpChip, fixed = "hday", P3D = F, plot = F)

write.table(gwaBlue, file = "~/Dropbox/Research_Poland_Lab/AM Panel/R/rrBlup/gwaBLUES1718.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
#write.table(gwaBlueHddt, file = "~/Dropbox/Research_Poland_Lab/AM Panel/AM_Panel/AM_Panel/AM_Panel/OriginalData/gwaBLUESHDDT.txt", sep = "\t",
#            quote = FALSE, row.names = FALSE, col.names = TRUE)

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("rs_number", "pos", "chrom")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","Figures","qqman_Figures",
                        paste("rrBLUP_2018_", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 70, height = 25)
    qqman::manhattan(dat,
                     main = paste("rrBLUP ",i),
                     chr = "chrom",
                     bp = "pos",
                     p = "P",
                     snp = "rs_number",
                     col = c("blue","grey40","black"),
                     chrlabs = c("1A","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     logp = F,
                     cex.axis = 2.5,
                     cex.lab = 3,
                     cex.main= 4,
                     cex = 3)
    dev.off()
    
  }
}

qqman.plot(gwaBlue)

bluesNA<- na.omit(blues)
bluesNA$Taxa<- toupper(bluesNA$Taxa)

write.table(bluesNA, 
            "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res2<- Hmisc::rcorr(phenoDat)
cor<- flattenCorrMatrix(res2$r, res2$P)

cordist<- ggplot(data = cor, aes(cor)) +
  geom_histogram(binwidth = 0.1, fill = "white", colour = "black") +
  theme_bw()

cordist

write.table(cor, "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/BluesCorr2018.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)



