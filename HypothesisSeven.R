rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(asreml)
library(lme4)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

pheno17<- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18<- fread("./Phenotype_Database/pheno18_htpLong.txt")
phenoLong<- fread("./Phenotype_Database/Pheno_Long1718.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(phenoLong)

phenoLong<- phenoLong %>% 
  dplyr::rename(Plot_ID = entity_id) %>% 
  tidylog::select(Plot_ID,Variety,block,rep,range,column)

pheno17$Date<- as.Date(pheno17$Date,format = "%Y-%m-%d")
pheno17$Date<- format(pheno17$Date, "%Y%m%d")
pheno18$Date<- as.Date(pheno18$Date,format = "%Y-%m-%d")
pheno18$Date<- format(pheno18$Date, "%Y%m%d")

pheno17<- pheno17 %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20170331:RedEdge_20170609) %>% 
  tidylog::inner_join(phenoLong) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse()

pheno18<- pheno18 %>% 
  dplyr::rename(Plot_ID = entity_id)  %>% 
  unite("ID",c("ID","Date")) %>% 
  spread(key = ID, value = value) %>% 
  tidylog::select(Plot_ID,Variety,GRYLD,
                  GNDVI_20171120:RE_20180613) %>% 
  tidylog::inner_join(phenoLong) %>% 
  glimpse() %>% 
  distinct() %>% 
  glimpse()

###############################################################################
####                  ASREML to calculate heritability                  ####

pheno17$Plot_ID<- as.factor(pheno17$Plot_ID)
pheno17$Variety<- as.factor(pheno17$Variety)
pheno17$block<- as.factor(pheno17$block)
pheno17$rep<- as.factor(pheno17$rep)
pheno17$range<- as.factor(pheno17$range)
pheno17$column<- as.factor(pheno17$column)

pheno18$Plot_ID<- as.factor(pheno18$Plot_ID)
pheno18$Variety<- as.factor(pheno18$Variety)
pheno18$block<- as.factor(pheno18$block)
pheno18$rep<- as.factor(pheno18$rep)
pheno18$range<- as.factor(pheno18$range)
pheno18$column<- as.factor(pheno18$column)

#2017 trial
# Without the auto-reggression correlation

asreml.license.status()

t17<- asreml(fixed = GRVI_20171120 ~ 1,
             random = ~Variety + rep + rep:block,
             data = pheno18)
plot(t17)
coef(t17)$random
fitted(t17)
summary(t17)
resid(t17)

h<- as.data.frame(summary(t17)$varcomp)
h

h2<-as.data.frame(h[3,1] / (h[3,1] + (h[4,1]/2)))
h2

## 2017 all
effectvars <- names(pheno17) %in% c("block", "rep", "Variety", "year", "column", 
                                "range", "Plot_ID")
traits <- colnames(pheno17[ , !effectvars])
H2_2017<- data.frame(traits)
H2_2017$Heritability<- NA
fieldInfo<- pheno17 %>% 
  tidylog::select(Variety, rep, block, column, range)
ntraits<- 1:nrow(H2_2017)

for (i in ntraits) {
  print(paste("Working on trait", H2_2017[i,1]))
  j<- H2_2017[i,1]
  
  data<- cbind(fieldInfo, pheno17[,paste(j)])
  names(data)<- c("Variety","rep","block","column","range","Trait")
  
  t17<- asreml(fixed = Trait ~ 1,
               random = ~Variety + rep + rep:block,
               data = data)
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock17_",H2_2017[i,1],".pdf"))
  plot(t17)
  dev.off()
  h<- as.data.frame(summary(t17)$varcomp)
  print(paste("Creating Data Frame", j))
  print(h)
  
  h2<- (h[3,1] / (h[3,1] + (h[4,1]/2)))
  h2
  H2_2017[i,2]<- h2
}

dev.off()
# lme4 comparison 2017
calcH2r <- function(dat, fill = NA, ...) {
  
  r=length(table(dat$rep))
  
  effectvars <- names(dat) %in% c("block", "rep", "Variety", "year", "column", 
                                  "range", "Plot_ID")
  
  t <- colnames(dat[ , !effectvars])
  
  for (i in t) {
    
    print(paste("Working on trait", i))
    h = as.data.frame(VarCorr(
      lmer(paste0(i, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)))
    H2= h[1,4] / (h[1,4] + (h[4,4] / r)) 
    print(H2)
  }
  
}

calcH2r(pheno17)

## 2018
dev.off()

t17<- asreml(fixed = GRYLD ~ 1,
             random = ~Variety + rep + rep:block,
             data = pheno18)
plot(t17)

h<- as.data.frame(summary(t17)$varcomp)
h

h2<-as.data.frame(h[3,1] / (h[3,1] + (h[4,1]/2)))
h2

## 2018 all
effectvars <- names(pheno18) %in% c("block", "rep", "Variety", "year", "column", 
                                    "range", "Plot_ID")
traits <- colnames(pheno18[ , !effectvars])
H2_2018<- data.frame(traits)
H2_2018$Heritability<- NA
fieldInfo<- pheno18 %>% 
  tidylog::select(Variety, rep, block, column, range)
ntraits<- 1:nrow(H2_2018)

for (i in ntraits) {
  print(paste("Working on trait", H2_2018[i,1]))
  j<- H2_2018[i,1]
  print(paste("Creating Data Frame", j))
  data<- cbind(fieldInfo, pheno18[,paste(j)])
  names(data)<- c("Variety","rep","block","column","range","Trait")
  
  t17<- asreml(fixed = Trait ~ 1,
               random = ~Variety + rep + rep:block,
               data = data)
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock18_",H2_2018[i,1],".pdf"))
  plot(t17)
  dev.off()
  h<- as.data.frame(summary(t17)$varcomp)
  print(h)
  
  h2<- (h[3,1] / (h[3,1] + (h[4,1]/2)))
  h2
  H2_2018[i,2]<- h2
}

# lme4 comparison 2018
calcH2r <- function(dat, fill = NA, ...) {
  
  r=length(table(dat$rep))
  
  effectvars <- names(dat) %in% c("block", "rep", "Variety", "year", "column", 
                                  "range", "Plot_ID")
  
  t <- colnames(dat[ , !effectvars])
  
  for (i in t) {
    
    print(paste("Working on trait", i))
    h = as.data.frame(VarCorr(
      lmer(paste0(i, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)))
    H2= h[1,4] / (h[1,4] + (h[4,4] / r)) 
    print(H2)
  }
  
}

calcH2r(pheno18)

H2_2017<- H2_2017 %>% 
  separate(traits, c("Trait","Date"), sep = "_")
H2_2017$Date<- as.Date(H2_2017$Date, format = "%Y%m%d")

H2_2017 %>% 
  tidylog::filter(Trait != "GRYLD") %>% 
  ggplot(aes(x = Date, y = Heritability)) +
  geom_point() + 
  facet_wrap(~Trait, scales = "free") +
  theme_bw() + 
  labs(title = "Broad-sense heritability of VI over Date 2016/2017") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 14))

H2_2018<- H2_2018 %>% 
  separate(traits, c("Trait","Date"), sep = "_")
H2_2018$Date<- as.Date(H2_2018$Date, format = "%Y%m%d")

H2_2018 %>% 
  tidylog::filter(Trait != "GRYLD") %>% 
  tidylog::filter(Trait != "height") %>% 
  tidylog::filter(Date > as.Date("20180101", format = "%Y%m%d")) %>% 
  ggplot(aes(x = Date, y = Heritability)) +
  geom_point() + 
  facet_wrap(~Trait, scales = "free") +
  theme_bw() + 
  scale_x_date(date_breaks = "10 days", date_labels = "%b%d") +
  labs(title = "Broad-sense heritability of VI over Date 2017/2018") +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 16), 
        title = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 14))

###### Examine weird residuals #####

t18<- asreml(fixed = GNDVI_20171215 ~ 1,
             random = ~Variety + rep + rep:block,
             data = pheno18)
plot(t18)

residuals<- setDT(as.data.frame(t18[["residuals"]]), keep.rownames = T)
names(residuals)<- c("plot","residual")
fits<- as.data.frame(fitted.asreml(t18))
names(fits)<- "fitted"
plotInfo<- pheno18 %>% 
  tidylog::select(block,rep,range,column) %>% 
  bind_cols(residuals) %>% 
  bind_cols(fits)

plotInfo %>% 
  ggplot(aes(x = fitted, y = residual, colour = rep)) +
  geom_point() + 
  theme_bw() +
  labs(x = "fitted", y = "Residual", 
       title = "Residual plots for GNDVI_20171215")

coef(t18)$random

h2<-as.data.frame(h[3,1] / (h[3,1] + (h[4,1]/2)))
h2
