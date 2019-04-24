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
library(RColorBrewer)
library(reshape2)


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
snpChip$rs_number<- as.character(snpChip$rs_number)
chromosomes<- 1:21

for (i in chromosomes) {
  chr <- i
  print(paste("working on ", chr))
  
  pos<- snpChip %>% 
    tidylog::select(rs_number,chrom,pos) %>% 
    tidylog::filter(chrom == chr) %>% 
    glimpse()
  
  markers<- snpChip %>% 
    tidylog::filter(chrom == chr) %>% 
    tidylog::select(-rs_number,-chrom,-pos) %>% 
    t() %>% 
    as.matrix()
  str(markers)
  
  rownames(markers) <- c()
  colnames(markers) <- pos$rs_number
  
  c<- rcorr(markers)
  fltC<- flattenCorrMatrix(c$r,c$P)
  str(fltC)
  fltC$row<- as.character(fltC$row)
  fltC$column<- as.character(fltC$column)
  
  fltC<- fltC %>% 
    inner_join(pos, by = c("row" = "rs_number")) %>% 
    dplyr::rename(Chr1 = chrom, Pos1 = pos) %>% 
    inner_join(pos, by = c("column" = "rs_number")) %>% 
    dplyr::rename(Chr2 = chrom, Pos2 = pos) %>% 
    mutate(Dist = abs(Pos1 - Pos2))
  
  write.table(fltC,paste("./Genotype_Database/MarkerCorrelations_",chr,".txt"),
              sep = "\t",quote = F,col.names = T,row.names = F)
  
  png(paste("./Figures/LD/MarkerCorrelationDensity_",chr,".png"),
      width = 1200, height = 1000, units = "px")
  #pdf(paste("./Figures/LD/MarkerCorrelationDensity_",chr,".pdf"))
  d <- ggplot(data = fltC, aes(x = cor)) +
    geom_density() +
    theme_bw() +
    labs(title = paste("Distribution of correlation between markers"),
         subtitle = paste("Chromosome ",chr))
  print(d)
  png(paste("./Figures/LD/MarkerCorrelationOverDistance_",chr,".png"),
      width = 1200, height = 1000, units = "px")
  #pdf(paste("./Figures/LD/MarkerCorrelationOverDistance_",chr,".pdf"))
  ld<- ggplot(data = fltC, aes(x = Dist, y = abs(cor), colour = p)) +
    geom_point() +
    theme_bw() + 
    geom_smooth() +
    scale_color_gradient2(low = "#762a83",
                          mid = "#f7f7f7",
                          high = "#5aae61",
                          midpoint = 0.25) +
    labs(title = "Correlation of markers over physical distance",
         subtitle = paste("Chromosome ",chr))
  print(ld)
  
  dev.off()
}

graphics.off()

chr<- 21

pos<- snpChip %>% 
  tidylog::select(rs_number,chrom,pos) %>% 
  tidylog::filter(chrom == chr) %>% 
  glimpse()

markers<- snpChip %>% 
  tidylog::filter(chrom == chr) %>% 
  tidylog::select(-rs_number,-chrom,-pos) %>% 
  t() %>% 
  as.matrix()
str(markers)

rownames(markers) <- c()
colnames(markers) <- pos$rs_number

c<- cor(markers)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri<- get_upper_tri(c)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap
png(paste("./Figures/LD/MarkerCorrelationBlocks_",chr,".png"),
    width = 1200, height = 1000, units = "px")

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#00441b", high = "#40004b", mid = "#f7f7f7", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_bw() + 
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(1,1500,250) ) +
  scale_y_continuous(breaks = seq(1,1500,250)) +
  coord_fixed() +
  labs(title = "Marker Correlations by position",
       subtitle = paste("Chromosome ",chr))
dev.off()


