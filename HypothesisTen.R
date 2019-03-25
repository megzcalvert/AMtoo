rm(list = objects()); ls()

library(MatrixModels)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(RMySQL)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

snpChip <- read_delim(
  "./Genotype_Database/90KsnpChipHapMap/AMsnpChipImputed.hmp.txt", 
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

write.table(snpChip, file="./Genotype_Database/SelectedImputedBeagle.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./Genotype_Database/SelectedImputedBeagle.txt", 
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

write.table(snpChip, file="./Genotype_Database/SelectedImputedBeagleNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./Genotype_Database/SelectedImputedBeagleNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")
snpMatrix<- t(snpChip[ , c(-1, -2, -3)])
snpMatrix %>% glimpse()
snpMatrix[1:10,1:10]

snpMatrix<- setDT(as.data.frame(snpMatrix),keep.rownames = T)
snpMatrix[1:10,1:10]

##### Phenotypes ####
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

