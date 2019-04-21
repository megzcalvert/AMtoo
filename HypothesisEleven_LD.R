rm(list = objects()); ls()

library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(beepr)
library(rrBLUP)
library(Hmisc)
library(MVN)


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

snpChip[snpChip == snpChip$allele_a] = 0
snpChip[snpChip == snpChip$allele_b] = 2
snpChip[snpChip == "H"] = 1
snpChip[snpChip == "C"] = NA
snpChip[snpChip == "A"] = NA
snpChip[snpChip == "T"] = NA
snpChip[snpChip == "G"] = NA
snpChip[snpChip == "-"] = NA
snpChip[snpChip == "."] = NA

snpChip[1:5,1:13]
snpChip<- snpChip[ ,c(1,4,5,13:311)]

write.table(snpChip, file="./Genotype_Database/SelectedImputedBeagleNumeric.txt",
            col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

snpChip = fread(file="./Genotype_Database/SelectedImputedBeagleNumeric.txt", 
                header=TRUE, check.names=F, sep = "\t")

chrSum<- plyr::count(snpChip, vars = "chrom")

##### LD analysis ####

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

chrSum
snpChip[1:5,1:5]
str(snpChip)

chr <- 1

pos<- snpChip %>% 
  tidylog::select(rs_number,chrom,pos) %>% 
  tidylog::filter(chrom == chr) %>% 
  glimpse()

markers<- snpChip %>% 
  tidylog::filter(chrom == chr) %>% 
  tidylog::select(-rs_number,-chrom,-pos) %>% 
  t() %>% 
  as.matrix()

rownames(markers) <- c()
colnames(markers) <- 1:ncol(markers)

c<- rcorr(markers)
fltC<- flattenCorrMatrix(c$r,c$P)
fltC$row<- as.integer(fltC$row)
fltC$column<- as.integer(fltC$column)
fltC<- fltC %>% 
  inner_join(pos, by = c("row" = "rs_number")) %>% 
  dplyr::rename(Chr1 = chrom, Pos1 = pos) %>% 
  inner_join(pos, by = c("column" = "rs_number")) %>% 
  dplyr::rename(Chr2 = chrom, Pos2 = pos) %>% 
  mutate(Dist = abs(Pos1 - Pos2)) %>% 
  glimpse()
fltC %>% 
  ggplot(aes(x = cor)) +
  geom_density() +
  theme_bw() +
  labs(title = paste("Distribution of correlation between markers"),
       subtitle = paste("Chromosome ",chr))
fltC %>% 
  ggplot(aes(x = Dist, y = abs(cor), colour = p)) +
  geom_point() +
  theme_bw() + 
  geom_smooth() +
  scale_color_gradient2(low = "#762a83",
                        mid = "#f7f7f7",
                        high = "#5aae61",
                        midpoint = 0.25) +
  labs(title = "Correlation of markers over physical distance",
       subtitle = paste("Chromosome ",chr))


