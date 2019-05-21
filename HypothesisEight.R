rm(list = objects()); ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(plotly)
library(rrBLUP)

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

pcaMethods::checkData(snpMatrix)  #Check PCA assumptions

pcaAM<- pcaMethods::pca(snpMatrix, nPcs = 15) #SVD PCA

sumPCA<- as.data.frame(summary(pcaAM))

lineInfo <- fread(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/LineDetailsAMPanel.txt", 
  header = T, check.names = F, sep = '\t', data.table = F)
lineInfo$Name <- tolower(lineInfo$Name)
lineInfo$Name<- str_replace_all(lineInfo$Name, " ", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "-", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "'", "")

program<- lineInfo[,c("Name","Program")]

Scores<- as.data.frame(pcaMethods::scores(pcaAM)) 
Scores<- setDT(Scores, keep.rownames = TRUE)
Scores<- left_join(Scores, program, by = c("rn" = "Name"))
Scores$Program[is.na(Scores$Program)] = "NA"
Scores$Program<- as.factor(Scores$Program)

pca.plot <- function(x, p, q, ...) {
  
  plots<-ggplot(data = x, aes_string(x = p, y = q)) +
    geom_point(position = "jitter",aes(colour = Program)) +
    theme_bw() +
    labs(title = paste0("PCA Plot ", p, " and ", q), 
         x = paste(p, "R2 = ",(round(sumPCA[1,paste0(p)],3))*100,"%"), 
         y = paste(q, "R2 = ",(round(sumPCA[1,paste0(q)],3))*100,"%")) +
    scale_color_manual(values = c("#e6194B", "#3cb44b", "#a9a9a9", "#4363d8", 
                                  "#f58231","#000000", "#42d4f4", "#f032e6", 
                                  "#fabebe", 
                                  "#469990", "#e6beff", "#9A6324", 
                                  "#800000", "#000075"))
  print(plots)
  
}

pca.plot(Scores, "PC1", "PC2")
pca.plot(Scores, "PC1", "PC3")
pca.plot(Scores, "PC2", "PC3")
pca.plot(Scores, "PC1", "PC4")
pca.plot(Scores, "PC2", "PC4")
pca.plot(Scores, "PC3", "PC4")


p <- plot_ly(Scores, x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~Program, size = 2) %>%
  add_markers(size = 2) %>%
  layout(scene = list(xaxis = list(title = 'PC1 R2 = 5.2%'),
                      yaxis = list(title = 'PC2 R2 = 4.5%'),
                      zaxis = list(title = 'PC3 R2 = 3.7%')))
p

ggplot(data = Scores, aes_string(x = "PC2", y = "PC3")) +
  geom_point(position = "jitter",aes(colour = Program)) +
  theme_bw() +
  labs(title = "PCA Plot PC2 and PC3", 
       x = expression(paste("PC2 ", "R"^{2}, "= 4.5%")), 
       y =  expression(paste("PC3 ", "R"^{2}, "= 3.7%"))) +
  scale_color_manual(values = c("#e6194B", "#3cb44b", "#767676", "#4363d8", 
                                "#f58231","#000000", "#42d4f4", "#f032e6", 
                                "#fabebe", "#469990", "#e6beff", "#9A6324", 
                                "#800000", "#000075"))

