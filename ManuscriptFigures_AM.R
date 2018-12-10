library(tidyverse)
library(data.table)
library(myprettyreport)
library(patchwork)
library(sysfonts)
library(ggpubr)
library(cowplot)
library(janitor)
library(extrafont)
library(sysfonts)
library(ggbiplot)
library(GGally)
library(qqman)
library(ggrepel)


quartzFonts(CMUBright = c("CMUBright-Roman", 
                          "CMUBright-Bold",  
                          "CMUBright-Oblique", 
                          "CMUBright-BoldOblique"))

font_paths()
font_add(family = "CMUBright",
         regular = "/Library/Fonts/cmunbmr.ttf")
font_families()

#######Heritability plots for HTP
her17<- fread(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/2017heritability.txt")

head(her17)
#Separate trait into date and HTP
her17<- separate(her17, "Trait",c("Trait","Date"), sep = "_")
str(her17)
#Adjust date format in case we want to use date as plotting factor
her17$Date<- as.Date(her17$Date, format = "%Y%m%d")
#Remove non-HTP traits
# her17<- her17[which(her17$Trait != c("BYDV","GRWT","MOIST","PTHT","TESTWT",
#                                      "AWNS","HDDT")),]
#Change names so everything fits
her17$Trait[her17$Trait == "RedEdge"] <- "RE"

her18<- fread(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/2018heritability.txt")

head(her18)

#Separate trait into date and HTP
her18<- separate(her18, "Trait",c("Date","Trait"), sep = "_")
str(her18)

#Date was placed first in these files with an X, remove X
her18$Date<- sub('.', '', her18$Date)
str(her18)

#Adjust date format in case we want to use date as plotting factor
her18$Date<- as.Date(her18$Date, format = "%Y%m%d")
her18<- filter(her18, Date != c("2018-05-16","2018-05-29"))
str(her18)

#Change names so everything fits
her18$Trait[her18$Trait == "Nir"] <- "NIR"

str(her17)

herPlot17<- ggplot(data = her17,
                   aes(x = Date,
                       y = heritability,
                       colour = Trait,
                       group = Trait)) +
  geom_line(data = filter(her17, Date != "2017-06-23"),
    size = 3) +
  geom_point(data = filter(her17, Date != "2017-06-23"),
             size = 4 ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(face = "bold",
                                  size = 24),
        axis.title = element_text(size = 20),
        plot.tag = element_text(size = 26),
        text = element_text(family = "CMUBright",
                            colour = "black")) +
  labs(title = "Broad-Sense Heritability for HTP Traits for 2017 Growing Season",
       tag = "A",
       x = "Date",
       y = "Broad-Sense Heritability") +
  scale_color_manual(values = c("#762a83","#af8dc3","#dbc0de","#c7c7c7",
                                "#b1e0a7","#7fbf7b","#1b7837")) +
  scale_x_date(date_labels = "%b %d",
               date_breaks = "2 weeks") +
  geom_text_repel(data = filter(her17, Date == "2017-06-23"),
            aes(x = Date,
                y = heritability,
                label = Trait),
            nudge_x = 1, 
            colour = "black") +
  geom_point(data = filter(her17, Date == "2017-06-23"),
             colour = "#adaeae")


herPlot17

herPlot18<- ggplot(data = her18,
                   aes(x = Date,
                     y = heritability,
                     colour = Trait,
                     group = Trait)) +
  geom_line(data = filter(her18, Date != "2018-06-27"), size = 3) +
  geom_point(data = filter(her18, Date != "2018-06-27"),size = 4 ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(face = "bold",
                                  size = 24),
        axis.title = element_text(size = 20),
        plot.tag = element_text(size = 26),
        text = element_text(family = "CMUBright",
                            colour = "black")) +
  labs(title = "Broad-Sense Heritability for HTP Traits for 2018 Growing Season",
       tag = "B",
       x = "Date",
       y = "Broad-Sense Heritability") +
  scale_color_manual(values = c("#762a83","#af8dc3","#dbc0de","#c7c7c7",
                                "#b1e0a7","#7fbf7b","#1b7837")) +
  scale_x_date(date_labels = "%b %d",
               date_breaks = "2 weeks") +
  geom_text_repel(data = filter(her18, Date == "2018-06-27"),
            aes(x = Date,
                y = heritability,
                label = Trait),
            nudge_x = 1, 
            colour = "black") +
  geom_point(data = filter(her18, Date == "2018-06-27"),
             colour = "#adaeae")


herPlot18

ggarrange(herPlot17,herPlot18,
          common.legend = T,
          legend = "right",
          ncol = 1,
          nrow = 2)


######## PCA of genetic markers

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

Scores[is.na(Scores)] <- "Unknown"

pcaPlt1<-ggplot(data = Scores, aes(x = PC1, y = PC2, colour = Program)) +
  geom_point(alpha = 1,position = "jitter") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(face = "bold",
                                  size = 14),
        axis.title = element_text(size = 12),
        plot.tag = element_text(size = 16),
        text = element_text(family = "CMUBright",
                            colour = "black")) +
  labs(title = "PCA Plot of PC1 and PC2",
       tag = "A",
       x = paste("PC1 ", "R2 = ",(round(sumPCA[1,"PC1"],3))*100,"%"),
       y =  paste("PC2 ", "R2 = ",(round(sumPCA[1,"PC2"],3))*100,"%")) + 
  scale_color_manual(name="Breeding\nProgram",
                     values = c('#40004b','#762a83','#9970ab',"#af8dc3",
                                '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
                                "#bbe4b2",'#9cd795','#5aae61','#1b7837',
                                "#808080",'#00441b'))
pcaPlt1

pcaPlt2<-ggplot(data = Scores, aes(x = PC1, y = PC3, colour = Program)) +
  geom_point(alpha = 1,position = "jitter") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(face = "bold",
                                  size = 14),
        axis.title = element_text(size = 12),
        plot.tag = element_text(size = 16),
        text = element_text(family = "CMUBright",
                            colour = "black")) +
  labs(title = "PCA Plot of PC1 and PC3",
       tag = "B",
       x = paste("PC1 ", "R2 = ",(round(sumPCA[1,"PC1"],3))*100,"%"),
       y =  paste("PC3 ", "R2 = ",(round(sumPCA[1,"PC3"],3))*100,"%")) + 
  scale_color_manual(name="Breeding\nProgram",
                     values = c('#40004b','#762a83','#9970ab',"#af8dc3",
                                '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
                                "#bbe4b2",'#9cd795','#5aae61','#1b7837',
                                "#808080",'#00441b'))

pcaPlt2

pcaPlt3<-ggplot(data = Scores, aes(x = PC2, y = PC3, colour = Program)) +
  geom_point(alpha = 1,position = "jitter") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(face = "bold",
                                  size = 14),
        axis.title = element_text(size = 12),
        plot.tag = element_text(size = 16),
        text = element_text(family = "CMUBright",
                            colour = "black")) +
  labs(title = "PCA Plot of PC2 and PC3",
       tag = "C",
       x = paste("PC2 ", "R2 = ",(round(sumPCA[1,"PC2"],3))*100,"%"),
       y =  paste("PC3 ", "R2 = ",(round(sumPCA[1,"PC3"],3))*100,"%")) + 
  scale_color_manual(name="Breeding\nProgram",
                     values = c('#40004b','#762a83','#9970ab',"#af8dc3",
                                '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
                                "#bbe4b2",'#9cd795','#5aae61','#1b7837',
                                "#808080",'#00441b'))

pcaPlt3

ggarrange(pcaPlt1, pcaPlt2, pcaPlt3,
          nrow = 2,
          ncol = 1,
          common.legend = T,
          legend = "right")

prow <- plot_grid( pcaPlt1 + theme(legend.position="none"),
                   pcaPlt2 + theme(legend.position="none"),
                   pcaPlt3 + theme(legend.position="none"),
                   nrow = 2,
                   ncol = 2
)
prow

legend <- get_legend(pcaPlt1)

p <- plot_grid( prow, legend, rel_widths = c(3, .3))
p


################ PCA of blues pheno

blues2018 <- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
                        header = T, sep = "\t")
blues2018 <- as.data.frame( blues2018[, !names(blues2018) %in% c("NormGRWT","TESTWT","GRYLD")])

blues2018 <- rename(blues2018, 
                    replace = c("GRWT" = "GRWT_2018",
                                "MOIST" = "MOSIT_2018",
                                "PTHT" = "PTHT_2018",
                                "SPNAREA" = "SPNAREA_2018",
                                "awns" = "AWNS_2018",
                                "hday" = "HDDT_2018"))

lineInfo <- fread(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/LineDetailsAMPanel.txt", 
  header = T, check.names = F, sep = '\t', data.table = F)
lineInfo$Name <- toupper(lineInfo$Name)
lineInfo$Name<- str_replace_all(lineInfo$Name, " ", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "-", "_")
lineInfo$Name<- str_replace_all(lineInfo$Name, "'", "")

program<- lineInfo[,c(1,4)]

b2018<- blues2018[,2:ncol(blues2018)]

blues2018<- left_join(blues2018, program, by = c("Taxa" = "Name"))

b.program<- blues2018[,97]

b28.pca<- prcomp(b2018, scale. = T)

summary(b28.pca)

ggscreeplot(b28.pca)

bluesPCA12018<- ggbiplot2(
  b28.pca,
  choices = 1:2,
  scale = 1,
  groups = b.program,
  varname.size = 2, 
  varname.adjust = 1.25,
  color = "grey", 
  linetype = "dashed", 
  alpha_arrow = 0.5,
  textcolor = "grey",
  pointsize = 1
) + 
  scale_color_manual(
    name="Breeding\nProgram",
    values = c('#40004b','#762a83','#9970ab',"#af8dc3",
               '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
               "#bbe4b2",'#9cd795','#5aae61','#1b7837',
               "#808080",'#00441b')
  ) +
  labs(
    title = "PCA of Phenotypic BLUEs 2018",
    tag = "D",
    x = "Standardized PC1 (33.4% explained variance)",
    y = "Standardized PC2 (10.5% explained variance)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    legend.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 14,colour = "black"
    ),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.tag = element_text(size = 16,
                            colour = "black"
    ),
    text = element_text(family = "CMUBright")
  )
bluesPCA12018

bluesPCA22018<- ggbiplot2(
  b28.pca,
  choices = c(1,3),
  scale = 1.25,
  groups = b.program,
  varname.size = 2, 
  varname.adjust = 1.25,
  color = "grey", 
  linetype = "dashed", 
  alpha_arrow = 0.5,
  textcolor = "grey",
  pointsize = 1
) + 
  scale_color_manual(
    name="Breeding\nProgram",
    values = c('#40004b','#762a83','#9970ab',"#af8dc3",
               '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
               "#bbe4b2",'#9cd795','#5aae61','#1b7837',
               "#808080",'#00441b')
  )  +
  labs(
    title = "PCA of Phenotypic BLUEs 2018",
    tag = "E",
    x = "Standardized PC1 (33.4% explained variance)",
    y = "Standardized PC3 (9.0% explained variance)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    legend.text = element_text(
      size = 10,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = 14,
      colour = "black"
    ),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.tag = element_text(
      size = 16,
      colour = "black"
    ),
    text = element_text(family = "CMUBright")
  )
bluesPCA22018

bluesPCA32018<- ggbiplot2(
  b28.pca,
  choices = 2:3,
  scale = 1.5,
  groups = b.program,
  varname.size = 2, 
  varname.adjust = 1.25,
  color = "grey", 
  linetype = "dashed", 
  alpha_arrow = 0.5,
  textcolor = "grey",
  pointsize = 1
) + 
  scale_color_manual(
    name="Breeding\nProgram",
    values = c('#40004b','#762a83','#9970ab',"#af8dc3",
               '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
               "#bbe4b2",'#9cd795','#5aae61','#1b7837',
               "#808080",'#00441b')
  ) +
  labs(
    title = "PCA of Phenotypic BLUEs 2018",
    tag = "F",
    x = "Standardized PC2 (10.5% explained variance)",
    y = "Standardized PC3 (9.0% explained variance)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(
      size = 11,
      colour = "black"),
    legend.text = element_text(
      size = 10,
      colour = "black"),
    legend.title = element_text(
      size = 12,
      colour = "black"),
    plot.title = element_text(
      face = "bold",
      size = 14,
      colour = "black"),
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    plot.tag = element_text(
      size = 16,
      colour = "black"
    ),
    text = element_text(family = "CMUBright")
  )
bluesPCA32018

blues2017 <- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
                        header = T, sep = "\t")
blues2017 <- as.data.frame( blues2017[, !names(blues2017) %in% 
                                        c("NormGRWT","TESTWT","GRYLD")])

blues2017 <- rename(blues2017, 
                    replace = c("GRWT" = "GRWT_2017",
                                "MOIST" = "MOSIT_2017",
                                "PTHT" = "PTHT_2017",
                                "BYDV" = "BYDV_2017",
                                "awns" = "AWNS_2017",
                                "hday" = "HDDT_2017"))

b2017<- blues2017[,2:ncol(blues2017)]

blues2017<- left_join(blues2017, program, by = c("Taxa" = "Name"))

b.program<- blues2017[,73]

b27.pca<- prcomp(b2017, scale. = T)

summary(b27.pca)

ggscreeplot(b27.pca)

bluesPCA12017<- ggbiplot2(
  b27.pca,
  choices = 1:2,
  scale = 1,
  groups = b.program,
  varname.size = 2, 
  varname.adjust = 1.25,
  color = "grey", 
  linetype = "dashed", 
  alpha_arrow = 0.5,
  textcolor = "grey",
  pointsize = 1) + 
  scale_color_manual(
    name="Breeding\nProgram",
    values = c('#40004b','#762a83','#9970ab',"#af8dc3",
               '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
               "#bbe4b2",'#9cd795','#5aae61','#1b7837',
               "#808080",'#00441b')
  ) +
  labs(
    title = "PCA of Phenotypic BLUEs 2017",
    tag = "A",
    x = "Standardized PC1 (30.7% explained variance)",
    y = "Standardized PC2 (18.9% explained variance)"
  ) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(size = 12),
        plot.tag = element_text(size = 16),
        text = element_text(
          family = "CMUBright",
          colour = "black"
        ))
bluesPCA12017

bluesPCA22017<- ggbiplot2(b27.pca,
                          choices = c(1,3),
                          scale = 1.25,
                          groups = b.program,
                          varname.size = 2, 
                          varname.adjust = 1.25,
                          color = "grey", 
                          linetype = "dashed", 
                          alpha_arrow = 0.5,
                          textcolor = "grey") + 
  scale_color_manual(name="Breeding\nProgram",
                     values = c('#40004b','#762a83','#9970ab',"#af8dc3",
                                '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
                                "#bbe4b2",'#9cd795','#5aae61','#1b7837',
                                "#808080",'#00441b'))  +
  labs(
    title = "PCA of Phenotypic BLUEs 2017",
    tag = "B",
    x = "Standardized PC1 (33.4% explained variance)",
    y = "Standardized PC3 (11.2% explained variance)"
  ) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(
          face = "bold",
          size = 14),
        axis.title = element_text(size = 12),
        plot.tag = element_text(size = 16),
        text = element_text(family = "CMUBright",
                            colour = "black"))
bluesPCA22017

bluesPCA32017<- ggbiplot2(b27.pca,
                          choices = 2:3,
                          scale = 1.5,
                          groups = b.program,
                          varname.size = 2, 
                          varname.adjust = 1.25,
                          color = "grey", 
                          linetype = "dashed", 
                          alpha_arrow = 0.5,
                          textcolor = "grey") + 
  scale_color_manual(name="Breeding\nProgram",
                     values = c('#40004b','#762a83','#9970ab',"#af8dc3",
                                '#c2a5cf','#cfacd4',"#dbc0de",'#c7c7c7',
                                "#bbe4b2",'#9cd795','#5aae61','#1b7837',
                                "#808080",'#00441b')) +
  labs(
    title = "PCA of Phenotypic BLUEs 2017",
    tag = "C",
    x = "Standardized PC2 (18.9% explained variance)",
    y = "Standardized PC3 (11.2% explained variance)"
  ) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(
          face = "bold",
          size = 14
        ),
        axis.title = element_text(size = 12),
        plot.tag = element_text(size = 16),
        text = element_text(family = "CMUBright",
                            colour = "black"))
bluesPCA32017

prow <- ggarrange(
  bluesPCA12017,
  bluesPCA22017,
  bluesPCA32017,
  bluesPCA12018,
  bluesPCA22018,
  bluesPCA32018,
  nrow = 4,
  ncol = 2,
  align = "hv",
  legend = "right",
  common.legend = T) 
prow

##### BLUEs Summary

#Correlation plots with GGally

head(blues2017)
head(blues2018)

correlations<- function(x, ...) {
  
  md<- names(x) %in% c("rn","Taxa","year","rep","block","column",
                       "range", "entity_id","HDDT","Program")
  corDa<-x[ , !md]
  nums <- sapply(corDa, is.numeric)
  corDa<- corDa[ , nums]
  ggcorr(corDa, 
         hjust = 1,
         size = 3,
         label_color = "grey",
         low = "#40004b",
         mid = "#c7c7c7",
         high = "#00441b",
         midpoint = 0,
         nbreaks = 8,
         layout.exp = 10,
         family = "CMUBright") 
}

cor17<-correlations(blues2017)
cor17
cor18<-correlations(blues2018, layout.exp = 11)
cor18

ggarrange(cor17,cor18,
          ncol = 1,
          nrow = 2,
          common.legend = T,
          legend = "right")

distBlues2017<- ggpairs(blues2017,
                        mapping = NULL,
                        columns = c(2:5,7),
                        title = "2017 Phenotypic BLUES",
                        upper = list(continuous = wrap("cor", size = 5, 
                                                       color = "black",
                                                       family = "CMUBright")),
                        columnLabels = c("BYDV","GRWT","MOIST","PTHT","HDDT"),
                        lower = list(continuous = "smooth"),
                        diag = list(continuous = "densityDiag", 
                                    discrete = "barDiag"),
                        # params = NULL,
                        #xlab = NULL,
                        #ylab = NULL,
                        axisLabels = c("show") #, ,show "internal", "none"
                        
) +
  theme(text = element_text(family = "CMUBright",size = 14),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        panel.border = element_rect(fill = NA, linetype = "solid")) 

distBlues2017

distBlues2018<- ggpairs(blues2018,
                        mapping = NULL,
                        columns = c(2:5,7),
                        title = "2018 Phenotypic BLUES",
                        upper = list(continuous = wrap("cor", size = 5, 
                                                       color = "black",
                                                       family = "CMUBright")),
                        columnLabels = c("GRWT","MOIST","PTHT",
                                         "SPNAREA","HDDT"),
                        lower = list(continuous = "smooth"),
                        diag = list(continuous = "densityDiag", 
                                    discrete = "barDiag"),
                        # params = NULL,
                        #xlab = NULL,
                        #ylab = NULL,
                        axisLabels = c("show") #, ,show "internal", "none"
                        
) +
  theme(text = element_text(family = "CMUBright",size = 14),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) 

distBlues2018


######## BLUES distribution histograms

histFacet.plot <- function(x, results, info, ...) {
  
  md <- names(x) %in% c("rn","Taxa","year","rep","block","column",
                        "range", "entity_id")
  traits <- names(x[ , !md])
  plotList = list()
  for (i in traits) {
    thisPlot <- ggplot(data = x, aes_string(x = i)) + 
      geom_histogram(colour="black", 
                     fill="white",
                     size = 1) + 
      theme_bw() +
      xlab(paste(i)) +
      ylab("Frequency") +
      theme(text = element_text(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 24,
                                     colour = "black"),
            axis.title = element_text(size = 30),
            strip.text = element_text(size = 16),
            axis.line = element_line(colour = "black",
                                     size = 1)) 
    
    
    plotList[[i]] = thisPlot
    print(i)
    
  }
  return(plotList)
}


distHist2017<- histFacet.plot(blues2017[,1:72],' ',"")
distHist2018<- histFacet.plot(blues2018[,1:96],"","")

distHist<- c(distHist2017,distHist2018)

ggarrange(plotlist = distHist,
          ncol = 3,
          nrow = 3) %>%
  ggexport(filename = "~/Dropbox/Research_Poland_Lab/AM Panel/AMPanel_Manuscript/Supplementary/SupplementaryFigure1_Distributions.pdf",
           width = 50,
           height = 50)

##### Manhattan and QQ plots

dataRRblup<- fread("~/Dropbox/Research_Poland_Lab/AM Panel/R/rrBlup/gwaBLUES_all_rrBlup.txt", header = T)
head(dataRRblup)

<<<<<<< HEAD
dataMultiphen<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Univariate/resultsSingle52018.wide.txt",
=======
<<<<<<< HEAD
dataMultiphen<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Univariate/resultsSingle52018.wide.txt",
=======
dataMultiphen<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Univariate/resultsSingle52017.wide.txt",
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
                           header = T, sep = "")

dataMultiphen<- separate(dataMultiphen, "rsid", c("chrom","pos"), sep = "_")
dataMultiphen$chrom <- as.numeric(dataMultiphen$chrom)
dataMultiphen$pos <- as.numeric(dataMultiphen$pos)
dataMultiphen<- dataMultiphen[which(dataMultiphen$label == "pval"), 
                              c(5,3:4,6:ncol(dataMultiphen))]
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
colnames(dataMultiphen)[colnames(dataMultiphen)=="GRWT"] <- "GRWT_2018"
colnames(dataMultiphen)[colnames(dataMultiphen)=="GRYLD"] <- "GRYLD_2018"
colnames(dataMultiphen)[colnames(dataMultiphen)=="awns"] <- "awns_2018"
colnames(dataMultiphen)[colnames(dataMultiphen)=="hday"] <- "HDDT_2018"
colnames(dataMultiphen)[colnames(dataMultiphen)=="MOIST"] <- "MOIST_2018"
colnames(dataMultiphen)[colnames(dataMultiphen)=="PTHT"] <- "PTHT_2018"
colnames(dataMultiphen)[colnames(dataMultiphen)=="TESTWT"] <- "TESTWT_2018"
<<<<<<< HEAD
=======
=======
colnames(dataMultiphen)[colnames(dataMultiphen)=="GRWT"] <- "GRWT_2017"
colnames(dataMultiphen)[colnames(dataMultiphen)=="GRYLD"] <- "GRYLD_2017"
colnames(dataMultiphen)[colnames(dataMultiphen)=="awns"] <- "awns_2017"
colnames(dataMultiphen)[colnames(dataMultiphen)=="hday"] <- "HDDT_2017"
colnames(dataMultiphen)[colnames(dataMultiphen)=="MOIST"] <- "MOIST_2017"
colnames(dataMultiphen)[colnames(dataMultiphen)=="PTHT"] <- "PTHT_2017"
colnames(dataMultiphen)[colnames(dataMultiphen)=="TESTWT"] <- "TESTWT_2017"
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100

#traits<- colnames(pheno[,3:ncol(pheno)])
chrBP<- dataRRblup[,c("rs_number","chrom","pos")]

######### Reading in files as a list of data frames
#Linear
<<<<<<< HEAD
fileNames<- list.files(path = "/Users/megzcalvert/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate",
=======
<<<<<<< HEAD
fileNames<- list.files(path = "/Users/megzcalvert/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate",
                       full.names = T,
                       pattern = ".wide.txt$")

#Names From Linear
traitNames<- basename(fileNames) %>%
  str_remove_all(c(".wide.txt"))
=======
fileNames<- list.files(path = "/Users/megzcalvert/Dropbox/Research_Poland_Lab/AM Panel/plink/amPanel/amPMult/",
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
                       full.names = T,
                       pattern = ".mqfam.total$")

#Names From Linear
traitNames<- basename(fileNames) %>%
<<<<<<< HEAD
  str_remove_all(c(".wide.txt"))
=======
  str_remove_all(c(".mqfam.total"))
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100

## File loading function
load.file<- function (filename) {
  d<- fread(file = filename,header = TRUE,check.names = F,data.table = F)
  d
}

#Read in data
data<- lapply(fileNames, load.file)
#perms<- lapply(filePerms, load.file)

#Adjust names and join
names(data)<- traitNames


## Sort data by SNPs or whatever GAPIT needs this
sortedData<- lapply(data, function(df) {
  df[order(df$SNP),]
})


##### Extract variable of interest from each dataframe in list and place into 1 dataframe

<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
circData <- map_dfc(data,`[`, c("label","rsid","SNP","JointModel")) #`[`, c("label","rsid","SNP","JointModel")) #Pick whichever columns, MultiPhen: (data,`[`, c("label","rsid","SNP","JointModel") or just pvalue

circData <- circData[which(circData$label == "pval"),]
circData <- circData %>%
  separate("rsid", c("chrom","pos"), sep = "_") 
bp<- circData[,2:4]
str(bp)
bp$chrom<- as.numeric(bp$chrom)
bp$pos<- as.numeric(bp$pos)
str(bp)
circData <- circData %>%
  select(starts_with("JointModel"))
names(circData) <- traitNames

<<<<<<< HEAD
## Combine with position information obtained previously
Circos<- cbind(bp,circData)
=======
=======
circData <- map_df(sortedData,"P") #`[`, c("label","rsid","SNP","JointModel")) #Pick whichever columns, MultiPhen: (data,`[`, c("label","rsid","SNP","JointModel") or just pvalue
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a

# colnames(circData)[colnames(circData)=="GRWT"] <- "GRWT_2018"
# colnames(circData)[colnames(circData)=="GRYLD"] <- "GRYLD_2018"
# colnames(circData)[colnames(circData)=="awns"] <- "awns_2018"
# colnames(circData)[colnames(circData)=="hday"] <- "HDDT_2018"
# colnames(circData)[colnames(circData)=="MOIST"] <- "MOIST_2018"
# colnames(circData)[colnames(circData)=="PTHT"] <- "PTHT_2018"
# colnames(circData)[colnames(circData)=="TESTWT"] <- "TESTWT_2018"

## Combine with position information obtained previously
<<<<<<< HEAD
Circos<- cbind(bp,circData)
=======
Circos<- cbind(chrBP,circData)
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
Circos10<- Circos
Circos10[,4:ncol(Circos10)] <- -log10(Circos10[,4:ncol(Circos10)])


write.table(Circos10, file = "~/Dropbox/Research_Poland_Lab/AM Panel/AMPanel_Manuscript/Supplementary/PLINK/results_Multivariate-logPValue_PLINK.txt",
            quote = F, row.names = F, col.names = T) 

<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
Circos <- Circos  %>%
  select(-results2017BluBckSelec5,-results2017BluBckSelecVIF5,
         -results2017bluLasSelecML5,
         -results2018bluBckSelec5,-results2018bluBckSelecVIF5)
<<<<<<< HEAD

=======
=======
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
#### Manhattan and QQPlot 
#qqman plots are backup, trying more customisable soulutions

qqman.plot <- function(x, prog, ...) {
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
  md <- names(x) %in% c("SNP", "pos", "chrom")
  traits <- colnames(x[,4:ncol(x)])
  info<- x[,c("SNP", "pos", "chrom")]
  
  plotList = list()
  
  mypath <- file.path(
    "~/Dropbox/Research_Poland_Lab/AM Panel/AMPanel_Manuscript/Supplementary",
                      paste("ManhattanPlots_",prog,".pdf", 
                            sep = ""))
  pdf(file = mypath,
      onefile = T,
      paper = "a4r",
      family = "CMU Sans Serif",
      pointsize = 10)
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    par(mar=c(5,6,4,3)+0.1,
        mfrow = c(2,1))
    thisPlot<- qqman::manhattan(dat,
                                main = paste(prog,i),
                                chr = "chrom",
                                bp = "pos",
                                p = "P",
                                snp = "SNP",
                                col = c("#1b7837","#762a83","#c2a5cf"),
                                chrlabs = c("1A","1B","1D",
                                            "2A","2B","2D",
                                            "3A","3B","3D",
                                            "4A","4B","4D",
                                            "5A","5B","5D",
                                            "6A","6B","6D",
                                            "7A","7B","7D"),
                                genomewideline = -log10(0.05 / nrow(dat)),
                                suggestiveline = F,
                                logp = T,
                                #ylim = c(0,20),
                                cex.axis = 1,
                                cex.lab = 1.5,
                                cex.main= 2,
                                cex = 1)
<<<<<<< HEAD
=======
=======
  md <- names(x) %in% c("rs_number", "pos", "chrom")
  traits <- colnames(x[,4:ncol(x)])
  info<- x[,c("rs_number", "pos", "chrom")]
  
  plotList = list()
  pdf(file = 
        "~/Dropbox/Research_Poland_Lab/AM Panel/AMPanel_Manuscript/Supplementary/PLINK/AssociationPlots_PLINK.pdf",
      onefile = T)
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    par(mar=c(5,6,4,3)+0.1)
    thisPlot<- qqman::manhattan(dat,
                           main = paste(prog,i),
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
                           logp = T,
                           #ylim = c(0,20),
                           cex.axis = 1.5,
                           cex.lab = 2,
                           cex.main= 2.5,
                           cex = 1.5)
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
    plotList[[paste0(i,"_manhattan")]] = thisPlot
    print(paste(i,"_manhattan"))
    
    thisPlot<- qqman::qq(dat$P,
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
                         main = paste(prog,i),
                         cex.axis = 1,
                         cex.lab = 1.5,
                         cex.main= 2,
                         cex = 1)
<<<<<<< HEAD
=======
=======
                   main = paste(prog,i),
                   cex.axis = 1.2,
                   cex.lab = 1.5,
                   cex.main= 2,
                   cex = 1.5)
>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100
    plotList[[paste0(i,"_qq")]] = thisPlot
    print(paste(i,"_qq"))
    
  }
  return(plotList)
  dev.off()
}

<<<<<<< HEAD
plinkMulti<- qqman.plot(Circos,'MV-Multiphen')

dev.off()
=======
<<<<<<< HEAD
plinkMulti<- qqman.plot(Circos,'MV-Multiphen')

dev.off()
=======
plinkMulti<- qqman.plot(Circos,'PLINK')

>>>>>>> 392711150022b98d9adc03dd5f5ce4f857f8c58a
>>>>>>> 542c14909c22f425c5f365ce10b1149d97628100

