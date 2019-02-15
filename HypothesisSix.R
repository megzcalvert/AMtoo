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