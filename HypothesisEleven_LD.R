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

snpChip$chrom[snpChip$chrom == 1] = "1A"
snpChip$chrom[snpChip$chrom == 2] = "1B"
snpChip$chrom[snpChip$chrom == 3] = "1D"
snpChip$chrom[snpChip$chrom == 4] = "2A"
snpChip$chrom[snpChip$chrom == 5] = "2B"
snpChip$chrom[snpChip$chrom == 6] = "2D"
snpChip$chrom[snpChip$chrom == 7] = "3A"
snpChip$chrom[snpChip$chrom == 8] = "3B"
snpChip$chrom[snpChip$chrom == 9] = "3D"
snpChip$chrom[snpChip$chrom == 10] = "4A"
snpChip$chrom[snpChip$chrom == 11] = "4B"
snpChip$chrom[snpChip$chrom == 12] = "4D"
snpChip$chrom[snpChip$chrom == 13] = "5A"
snpChip$chrom[snpChip$chrom == 14] = "5B"
snpChip$chrom[snpChip$chrom == 15] = "5D"
snpChip$chrom[snpChip$chrom == 16] = "6A"
snpChip$chrom[snpChip$chrom == 17] = "6B"
snpChip$chrom[snpChip$chrom == 18] = "6D"
snpChip$chrom[snpChip$chrom == 19] = "7A"
snpChip$chrom[snpChip$chrom == 20] = "7B"
snpChip$chrom[snpChip$chrom == 21] = "7D"

chromosomes<- unique(snpChip$chrom)

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
  
  #png(paste0("./Figures/LD/MarkerCorrelationDensity_",chr,".png"),
      #width = 40, height = 30, units = "cm",
      #res = 320)
  #pdf(paste("./Figures/LD/MarkerCorrelationDensity_",chr,".pdf"))
  d <- ggplot(data = fltC, aes(x = cor)) +
    geom_density() +
    theme_bw() +
    labs(title = paste("Distribution of correlation between markers"),
         subtitle = paste("Chromosome ",chr))
  print(d)
  png(paste0("./Figures/LD/MarkerCorrelationOverDistance_",chr,".png"),
      width = 12.5, height = 5, units = "cm",
      res = 320)
  ld<- ggplot(data = fltC, aes(x = Dist, y = cor^2, colour = p)) +
    geom_point(size = 0.15) +
    theme_bw() + 
    geom_smooth() +
    scale_color_gradient2(low = "#762a83",
                          mid = "#f7f7f7",
                          high = "#5aae61",
                          midpoint = 0.25,
                          name = "p-value") +
    labs(subtitle = paste("Chromosome ",chr),
         x = "Distance (bp)",
         y = expression(paste("r"^2)))
  print(ld)
  
  dev.off()
}

graphics.off()

##### Heatmap type block things #####
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

chr<- 1

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
upper_tri<- get_upper_tri(c)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

str(melted_cormat)

png(paste0("./Figures/LD/MarkerCorrelationBlocks_",chr,".png"),
    width = 1200, height = 1000, units = "px")

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value^2))+
  geom_tile(color = "white",width = 1, height = 1) +
  #geom_point(size = 0.25) +
  scale_fill_gradient(low = "#f7f7f7", high = "#40004b", 
                      limit = c(0,1), space = "Lab", 
                      name="Pearson\nCorrelation") +
  theme_bw() + 
  theme(axis.text = element_blank(),
        panel.grid = element_blank()) +
  #coord_cartesian(xlim = c(0,594102056),ylim = c(0,594102056)) +
  #coord_fixed() +
  labs(title = "Marker Correlations by position",
       subtitle = paste("Chromosome ",chr))

dev.off()

melted_cormat$Var1<- as.character(melted_cormat$Var1)
melted_cormat$Var2<- as.character(melted_cormat$Var2)
melted_cormat<- melted_cormat %>%
  inner_join(pos, by = c("Var1" = "rs_number")) %>%
  dplyr::rename(Chr1 = chrom, Pos1 = pos) %>%
  inner_join(pos, by = c("Var2" = "rs_number")) %>%
  dplyr::rename(Chr2 = chrom, Pos2 = pos) %>%
  glimpse()


