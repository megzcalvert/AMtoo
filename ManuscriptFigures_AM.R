library(tidyverse)
library(data.table)
library(myprettyreport)
library(data.table)
library(patchwork)
library(sysfonts)
library(ggpubr)
library(cowplot)
library(janitor)
library(extrafont)
library(sysfonts)
library(ggbiplot)
library(GGally)

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
her17<- her17[which(her17$Trait != c("BYDV","GRWT","MOIST","PTHT","TESTWT",
                                     "AWNS","HDDT")),]
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
#Remove non-HTP traits
her18<- her18[which(her18$Trait != c("GRWT","GRYLD","MOIST","PTHT","SPNAREA",
                                     "TESTWT","AWNS","HDDT")),]
#Change names so everything fits
her18$Trait[her18$Trait == "Nir"] <- "NIR"

herPlot17<- ggplot(data = her17,
                   aes(x = factor(Date),
                       y = heritability,
                       colour = Trait,
                       group = Trait)) +
  geom_line(size = 3) +
  geom_point(size = 4 ) +
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
  scale_x_discrete(labels = c("Season","01Nov","08Nov","11Nov","23Nov",
                              "28Nov","01Dec","31Mar","06Apr",
                              "12Apr","17Apr","23Apr","03May",
                              "05May","12May","23May","02Jun",
                              "06Jun","09Jun")) +
  geom_hline(colour = "#b2b2b2",
             yintercept = c(0.5393841, 0.8172119,0.6478161,0.7612984,0.7803741,
                            0.9544459,0.9401581),
             linetype = "dashed") +
  annotate("text", x = 1, y = 0.55, label = "BYDV",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.821, label = "GRWT",
           family = "CMUBright", size = 6)  +
  annotate("text", x = 1, y = 0.65, label = "MOIST",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.75, label = "PTHT",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.79, label = "TESTWT",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.96, label = "AWNS",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.93, label = "HDDT",
           family = "CMUBright", size = 6)

herPlot17

herPlot18<- ggplot(data = her18,
                   aes(
                     x = factor(Date),
                     y = heritability,
                     colour = Trait,
                     group = Trait
                   )) +
  geom_line(size = 3) +
  geom_point(size = 4 ) +
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
  scale_x_discrete(labels = c("Season","20Nov","27Nov","05Dec",
                              "15Dec","18Dec","04Apr",
                              "12Apr","19Apr","23Apr",
                              "04May","14May","16May",
                              "29May","06Jun","13Jun")) +
  geom_hline(colour = "#b2b2b2",
             yintercept = c(0.5889003, 0.4652914,0.9464354,0.5794162,0.6473628,
                            1,0.8817521),
             linetype = "dashed") +
  annotate("text", x = 1, y = 0.6, label = "GRWT",
           family = "CMUBright", size = 6)  +
  annotate("text", x = 1, y = 0.465, label = "MOIST",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.95, label = "PTHT",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.57, label = "SPNAREA",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.655, label = "TESTWT",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 1, label = "AWNS",
           family = "CMUBright", size = 6) +
  annotate("text", x = 1, y = 0.88, label = "HDDT",
           family = "CMUBright", size = 6)

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
      geom_histogram(colour="black", fill="white") + 
      theme_bw() +
      xlab(i) +
      ylab("Frequency") +
      theme(text = element_text(family = "CMUBright",
                                size = 14,
                                colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            strip.text = element_text(size = 15)) 
    
    plotList[[i]] = thisPlot
    print(i)
    
  }
  return(plotList)
}

distHist2017<- histFacet.plot(blues2017[,1:72],' ',"")

ggarrange(plotlist = distHist2017,
          ncol = 3,
          nrow = 6,
          align = "hv")
