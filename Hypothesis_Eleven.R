rm(list = objects())
ls()

library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(asreml)
library(beepr)
library(rrBLUP)
library(Hmisc)
library(MVN)


getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

snpChip <- read_delim(
  "./Genotype_Database/90KsnpChipHapMap/AMsnpChipImputed.hmp.txt",
  "\t",
  escape_double = FALSE, trim_ws = TRUE
)
snpChip <- snpChip %>%
  clean_names()

missAmbiguous <- c("0", "+", "-")
hetCodes <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")
hapgeno <- as.matrix(snpChip[, 13:ncol(snpChip)])
hapgeno[hapgeno %in% missAmbiguous] <- NA
hapgeno[hapgeno == "N"] <- NA
hapgeno[hapgeno %in% hetCodes] <- "H"
snpChip <- cbind(snpChip[, 1:12], hapgeno)
rm(hapgeno)

write.table(snpChip,
            file = "./Genotype_Database/SelectedImputedBeagle.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

snpChip <- fread(
  file = "./Genotype_Database/SelectedImputedBeagle.txt",
  header = TRUE, check.names = F, sep = "\t"
)

snpChip[snpChip == snpChip$allele_a] <- 0
snpChip[snpChip == snpChip$allele_b] <- 2
snpChip[snpChip == "H"] <- 1
snpChip[snpChip == "C"] <- NA
snpChip[snpChip == "A"] <- NA
snpChip[snpChip == "T"] <- NA
snpChip[snpChip == "G"] <- NA
snpChip[snpChip == "-"] <- NA
snpChip[snpChip == "."] <- NA

snpChip <- snpChip[, c(1, 4, 5, 13:311)]

write.table(snpChip,
            file = "./Genotype_Database/SelectedImputedBeagleNumeric.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

snpChip <- fread(
  file = "./Genotype_Database/SelectedImputedBeagleNumeric.txt",
  header = TRUE, check.names = F, sep = "\t"
)

chrSum <- plyr::count(snpChip, vars = "chrom")
snpMatrix <- t(snpChip[, c(-1, -2, -3)])
snpMatrix %>% glimpse()
snpMatrix[1:10, 1:10]

snpChip[1:10, 1:10]

myGD <- setDT(as.data.frame(snpMatrix), keep.rownames = T) %>% 
  arrange(rn)

myGD[1:5, 1:5]
colnames(myGD)[colnames(myGD) == "rn"] <- "Taxa"
myGM <- snpChip[, 1:3]
colnames(myGM) <- c("Name", "Chromosome", "Position")

##### Phenotypes ####

pheno17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoLong <- fread("./Phenotype_Database/Pheno_Long171819.txt")

dat17 <- fread("./Phenotype_Database/Hyp10BLUEs_17.txt")
dat18 <- fread("./Phenotype_Database/Hyp10BLUEs_18.txt")
dat19 <- fread("./Phenotype_Database/Hyp10BLUEs_19.txt")

snpChip[1:10, 1:10]
snpMatrix[1:10, 1:10]

snpLines <- setDT(
  as.data.frame(colnames(snpChip[, 3:ncol(snpChip)]))
)
colnames(snpLines) <- "rn"
dat17 <- semi_join(dat17, snpLines, by = "rn") %>% 
  arrange(rn)
colnames(dat17)[colnames(dat17) == "rn"] <- "Taxa"

dat18 <- semi_join(dat18, snpLines, by = "rn") %>% 
  arrange(rn)
colnames(dat18)[colnames(dat18) == "rn"] <- "Taxa"

dat19 <- semi_join(dat19, snpLines, by = "rn") %>% 
  arrange(rn)
colnames(dat19)[colnames(dat19) == "rn"] <- "Taxa"

write.table(dat17,
            file = "./Phenotype_Database/Hyp11Blue_2017.txt", quote = F,
            sep = "\t", row.names = F, col.names = T
)
write.table(dat18,
            file = "./Phenotype_Database/Hyp11Blue_2018.txt", quote = F,
            sep = "\t", row.names = F, col.names = T
)
write.table(dat19,
            file = "./Phenotype_Database/Hyp11Blue_2019.txt", quote = F,
            sep = "\t", row.names = F, col.names = T
)

numberedCols <- paste(1:(ncol(myGD) - 1), sep = ",")

myGM <- myGM %>%
  tidylog::mutate(V = "V") %>%
  unite(Name, c("V", "Name"), sep = "") %>%
  glimpse()

write.table(myGD, "./R/Gapit/HypothesisEleven/myGD.txt",
            sep = "\t", quote = F,
            col.names = T, row.names = F
)
write.table(myGM, "./R/Gapit/HypothesisEleven/myGM.txt",
            sep = "\t", quote = F,
            col.names = T, row.names = F
)

##### rrBLUP ####
library(rrBLUP)
library(tidyverse)

setwd("~/Dropbox/Research_Poland_Lab/AM Panel")

dat_awns <- read.table(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/ASREMLBlup_awns.txt",
  head = TRUE
)
dat_awns[1:10, ]
dat_awns <- dat_awns %>%
  semi_join(snpLines, by = c("Taxa" = "rn")) %>%
  arrange(var = "Taxa") %>%
  tidylog::select(Taxa, phenotype_value) %>% 
  group_by(Taxa) %>% 
  summarise(Awns = median(phenotype_value)) 

myGD_awns <- myGD %>% 
  semi_join(dat_awns)

snpMarkers <- t(myGD_awns) %>% 
  row_to_names(row_number = 1) %>% 
  as.matrix()

snpMarkers[snpMarkers == 0] <- -1
snpMarkers[snpMarkers == 1] <- 0
snpMarkers[snpMarkers == 2] <- 1

snpRR <- cbind(myGM, snpMarkers)

write.table(snpRR,
            file = "./Genotype_Database/Intermediate_awns_geno.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

snpRR <- fread(
  file = "./Genotype_Database/Intermediate_awns_geno.txt",
  header = TRUE, check.names = F, sep = "\t"
)

print(paste('Check that marker matrix and phenotypes align',  
            all(colnames(snpRR[,4:ncol(snpRR)]) == dat_awns$Taxa), sep = ' '))

res <- GWAS(
  pheno = dat_awns,
  geno = snpRR, 
  n.PC = 4,
  P3D = F
)

write.table(res,
            file = "./R/rrBlup/HypothesisEleven/awns_Allyrs_4PC.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T
)

effectvars <- c("Taxa")
traits <- dat17 %>%
  dplyr::select(!any_of(effectvars)) %>%
  colnames()

for (i in traits) {
  print(paste("working on ", i))
  
  dat <- dat17 %>%
    tidylog::select("Taxa", paste(i))
  
  png(
    filename =
      paste0("./R/rrBlup/HypothesisEleven/PC4_2017/", i, "_2017_4PC.png")
  )
  
  res <- GWAS(
    pheno = dat,
    geno = snpRR, n.PC = 4,
    P3D = F
  )
  dev.off()
  
  write.table(res,
              file = paste0(
                "./R/rrBlup/HypothesisEleven/PC4_2017/",
                i, "_2017_4PC.txt"
              ), quote = F, sep = "\t",
              row.names = F, col.names = T
  )
}

graphics.off()
dev.off()

traits <- dat18 %>%
  dplyr::select(!any_of(effectvars)) %>%
  colnames()

for (i in traits) {
  print(paste("working on ", i))
  
  dat <- dat18 %>%
    tidylog::select("Taxa", paste(i))
  
  png(
    filename =
      paste0("./R/rrBlup/HypothesisEleven/PC4_2018/", i, "_2018_4PC.png")
  )
  res <- GWAS(
    pheno = dat,
    geno = snpRR, n.PC = 4,
    P3D = F
  )
  dev.off()
  write.table(res,
              file = paste0(
                "./R/rrBlup/HypothesisEleven/PC4_2018/",
                i, "_2018_4PC.txt"
              ), quote = F, sep = "\t",
              row.names = F, col.names = T
  )
}

traits <- dat19 %>%
  dplyr::select(!any_of(effectvars)) %>%
  colnames()

for (i in traits) {
  print(paste("working on ", i))
  
  dat <- dat19 %>%
    tidylog::select("Taxa", paste(i))
  
  png(
    filename =
      paste0("./R/rrBlup/HypothesisEleven/PC4_2019/", i, "_2019_4PC.png")
  )
  res <- GWAS(
    pheno = dat,
    geno = snpRR, n.PC = 4,
    P3D = F
  )
  dev.off()
  write.table(res,
              file = paste0(
                "./R/rrBlup/HypothesisEleven/PC4_2019/",
                i, "_2019_4PC.txt"
              ), quote = F, sep = "\t",
              row.names = F, col.names = T
  )
}

beep(9)

###############################################################################
#####       Comparing the GAPIT and rrBLUP results on the same scale       ####

######### Reading in files as a list of data frames

load.file <- function(filename) {
  d <- fread(file = filename, header = TRUE, check.names = F, data.table = F)
  d
}

## rrBLUP results

fileNames <- list.files(
  path = "./R/rrBlup/HypothesisEleven/PC4_2017",
  full.names = T,
  pattern = "_2017_4PC.txt$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all("_2017_4PC.txt")

rrB4_17 <- lapply(fileNames, load.file)

names(rrB4_17) <- traitNames

fileNames <- list.files(
  path = "./R/rrBlup/HypothesisEleven/PC4_2018",
  full.names = T,
  pattern = "_2018_4PC.txt$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all("_2018_4PC.txt")

rrB4_18 <- lapply(fileNames, load.file)

names(rrB4_18) <- traitNames

fileNames <- list.files(
  path = "./R/rrBlup/HypothesisEleven/PC4_2019",
  full.names = T,
  pattern = "_2019_4PC.txt$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all("_2019_4PC.txt")

rrB4_19 <- lapply(fileNames, load.file)

names(rrB4_19) <- traitNames

snpPos <- rrB4_17$GNDVI_20170331[, 1:3]

rrB4_17F <- snpPos %>%
  bind_cols(map_dfr(rrB4_17, 4)) %>%
  mutate(
    Name = as.character(Name),
    Position = as.numeric(Position)
  ) %>%
  arrange(Chromosome, Position)

rrB4_18F <- snpPos %>%
  bind_cols(map_dfr(rrB4_18, 4)) %>%
  mutate(
    Name = as.character(Name),
    Position = as.numeric(Position)
  ) %>%
  arrange(Chromosome, Position)

rrB4_19F <- snpPos %>%
  bind_cols(map_dfr(rrB4_19, 4)) %>%
  mutate(
    Name = as.character(Name),
    Position = as.numeric(Position)
  ) %>%
  arrange(Chromosome, Position)

## GAPIT results
# 
# fileNames <- list.files(
#   path = "./R/Gapit/HypothesisEleven/PC4_01_2017",
#   full.names = T,
#   pattern = "GWAS.Results.csv$"
# )
# traitNames <- basename(fileNames) %>%
#   str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")
# 
# gap4_17 <- lapply(fileNames, load.file)
# 
# names(gap4_17) <- traitNames
# 
# gap4_17F <- gap4_17$GNDVI_20170331[, 1:3] %>%
#   arrange(Chromosome, Position) %>%
#   bind_cols(map_dfr(gap4_17, 4)) %>%
#   dplyr::rename(chrom = Chromosome, pos = Position) %>%
#   mutate(pos = as.numeric(pos))
# 
# fileNames <- list.files(
#   path = "./R/Gapit/HypothesisEleven/PC4_01_2018",
#   full.names = T,
#   pattern = "GWAS.Results.csv$"
# )
# traitNames <- basename(fileNames) %>%
#   str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")
# 
# gap4_18 <- lapply(fileNames, load.file)
# 
# names(gap4_18) <- traitNames
# 
# gap4_18F <- gap4_18$GNDVI_20171120[, 1:3] %>%
#   arrange(Chromosome, Position) %>%
#   bind_cols(map_dfr(gap4_18, 4)) %>%
#   dplyr::rename(chrom = Chromosome, pos = Position) %>%
#   mutate(pos = as.numeric(pos))
# 
# fileNames <- list.files(
#   path = "./R/Gapit/HypothesisEleven/PC4_2019",
#   full.names = T,
#   pattern = "GWAS.Results.csv$"
# )
# traitNames <- basename(fileNames) %>%
#   str_remove_all("GAPIT.MLM.|.GWAS.Results.csv")
# 
# gap4_19 <- lapply(fileNames, load.file)
# 
# names(gap4_19) <- traitNames
# 
# gap4_19F <- gap4_19$GNDVI_20190103[, 1:3] %>%
#   arrange(Chromosome, Position) %>%
#   bind_cols(map_dfr(gap4_19, 4)) %>%
#   dplyr::rename(chrom = Chromosome, pos = Position) %>%
#   mutate(pos = as.numeric(pos))

##### Comparison plots between GWAS methods and PCs ####

don_rrB4_17 <- rrB4_17F %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_17F, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

axisdf <- don_rrB4_17 %>%
  group_by(chrom) %>%
  dplyr::summarize(
    center = (max(BPcum) + min(BPcum)) / 2,
    maximum = max(BPcum)
  )

don_rrB4_18 <- rrB4_18F %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_18F, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

don_rrB4_19 <- rrB4_19F %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_19F, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

don_gap4_17 <- gap4_17F %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gap4_17F, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

don_gap4_18 <- gap4_18F %>%
  # Compute chromosome size
  group_by(chrom) %>%
  summarise(chr_len = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gap4_18F, ., by = c("chrom" = "chrom")) %>%
  # Add a cumulative position of each SNP
  arrange(chrom, pos) %>%
  mutate(BPcum = pos + tot)

ggplot(don_rrB4_18, aes(x = BPcum)) +
  # Show all points
  geom_point(aes(y = -log10(don_rrB4_18$GNDVI_20180613)),
             colour = "#4d9221",
             alpha = 0.5, size = 1
  ) +
  geom_point(aes(y = -log10(don_gap4_18$GNDVI_20180613)),
             colour = "#c51b7d",
             alpha = 0.5, size = 1
  ) +
  # Significance Threshold
  geom_hline(yintercept = -log10(0.05 / nrow(don_rrB4_17)), linetype = 2) +
  # geom_vline(xintercept = )
  # custom X axis:
  scale_x_continuous(label = axisdf$chrom, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0.05)) + # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 16),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16)
  ) +
  labs(
    title = "GWAS results GNDVI_20180613",
    subtitle = "rrBLUP results - green, GAPIT results - pink, Bonferroni corection alpha = 0.05",
    x = "Chromosome",
    y = "-log10(P)"
  )

ggplot(don_rrB4_19, aes(x = BPcum, colour = as.factor(don_rrB4_19$chrom))) +
  # Show all points
  geom_point(aes(y = -log10(GRYLD)),
             alpha = 0.5, size = 1
  ) +
  scale_color_manual(values = rep(c("#2ca25f", "#8856a7", "#43a2ca"), 22)) +
  # Significance Threshold
  geom_hline(yintercept = -log10(0.05 / nrow(don_rrB4_19)), linetype = 2) +
  geom_vline(xintercept = 6979839046, linetype = 3) +
  # custom X axis:
  scale_x_continuous(
    label = c(
      "1A", "1B", "1D",
      "2A", "2B", "2D",
      "3A", "3B", "3D",
      "4A", "4B", "4D",
      "5A", "5B", "5D",
      "6A", "6B", "6D",
      "7A", "7B", "7D"
    ),
    breaks = axisdf$center
  ) +
  scale_y_continuous(expand = c(0, 0.05)) + # remove space between plot area and x axis
  coord_cartesian(ylim = c(0, 6.5)) +
  labs(
    title = "GWAS results GRYLD 2018-2019",
    subtitle = "Bonferroni Threshold alpha = 0.05",
    x = "Chromosome",
    y = "-log10(P)"
  ) +
  annotate(geom = "text", x = 7000000000, y = 6, label = "Rht-1B")

##### Making this a function to plot all ####

myManhattans <- function(dat, traits, saveFileFigures, identifier,
                         colourOne, colourTwo, colourThree, ...) {
  dat <- dat %>%
    # Compute chromosome size
    group_by(Chromosome) %>%
    summarise(chr_len = max(Position)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    tidylog::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(dat, ., by = c("Chromosome" = "Chromosome")) %>%
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate(BPcum = Position + tot) %>%
    select(1:3, tot, BPcum, everything())
  
  axisdf <- dat %>%
    group_by(Chromosome) %>%
    dplyr::summarize(
      center = (max(BPcum) + min(BPcum)) / 2,
      maximum = max(BPcum)
    )
  plotList <- list()
  
  for (i in traits) {
    thisPlot <- ggplot(
      data = dat,
      mapping = aes(
        x = BPcum,
        colour = as.factor(Chromosome)
      )
    ) +
      # Show all points
      geom_point(aes_string(y = paste(i)),
                 alpha = 0.5, size = 1
      ) +
      scale_color_manual(values = rep(
        c(colourOne, colourTwo, colourThree),
        22
      )) +
      # Significance Threshold
      geom_hline(yintercept = -log10(0.05 / nrow(dat)), linetype = 2) +
      # custom X axis:
      scale_x_continuous(
        label = c(
          "1A", "1B", "1D",
          "2A", "2B", "2D",
          "3A", "3B", "3D",
          "4A", "4B", "4D",
          "5A", "5B", "5D",
          "6A", "6B", "6D",
          "7A", "7B", "7D"
        ),
        breaks = axisdf$center
      ) +
      scale_y_continuous(expand = c(0, 0.05)) + # remove space between plot area and x axis
      # Custom the theme:
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.subtitle = element_text(size = 16,
                                     hjust = 0)
      ) +
      labs(
        title = paste0("GWAS results ", i),
        subtitle = "Bonferroni Threshold alpha = 0.05",
        x = "Chromosome",
        y = "-log10(P)"
      )
    
    plotList[[i]] <- thisPlot
    print(i)
    ggpubr::ggexport(thisPlot,
                     filename = paste0(
                       saveFileFigures,
                       identifier, "_", i, ".png"
                     ),
                     width = 5000, height = 1500,
                     res = 360
    )
  }
  return(plotList)
}

traits <- colnames(rrB4_17F[, 4:ncol(rrB4_17F)])
manhattan17 <- myManhattans(
  dat = rrB4_17F, traits = traits,
  saveFileFigures = "./Figures/AssociationPlots/",
  identifier = "rrBlup4PC_17",
  colourOne = "#2ca25f", colourTwo = "#43a2ca",
  colourThree = "#8856a7"
)

rrB4_18F<- rrB4_18F %>% 
  rename(GRYLD_2018 = GRYLD)

traits <- colnames(rrB4_18F[, 4:ncol(rrB4_18F)])
manhattan18 <- myManhattans(
  dat = rrB4_18F, traits = traits,
  saveFileFigures = "./Figures/AssociationPlots/",
  identifier = "rrBlup4PC_18",
  colourOne = "#2ca25f", colourTwo = "#43a2ca",
  colourThree = "#8856a7"
)

traits <- colnames(rrB4_19F[, 4:ncol(rrB4_19F)])
manhattan19 <- myManhattans(
  dat = rrB4_19F, traits = traits,
  saveFileFigures = "./Figures/AssociationPlots/",
  identifier = "rrBlup4PC_19",
  colourOne = "#2ca25f", colourTwo = "#43a2ca",
  colourThree = "#8856a7"
)

# traits <- colnames(gap4_17F[, 4:ncol(gap4_17F)])
# manhattan19 <- myManhattans(
#   dat = gap4_17F, traits = traits,
#   saveFileFigures = "./Figures/AssociationPlots/",
#   identifier = "Gapit4PC_17",
#   colourOne = "#2ca25f", colourTwo = "#43a2ca",
#   colourThree = "#8856a7"
# )
# 
# traits <- colnames(gap4_18F[, 4:ncol(gap4_18F)])
# manhattan18 <- myManhattans(
#   dat = gap4_18F, traits = traits,
#   saveFileFigures = "./Figures/AssociationPlots/",
#   identifier = "Gapit4PC_18",
#   colourOne = "#2ca25f", colourTwo = "#43a2ca",
#   colourThree = "#8856a7"
# )
# 
# traits <- colnames(gap4_19F[, 4:ncol(gap4_19F)])
# manhattan19 <- myManhattans(
#   dat = gap4_19F, traits = traits,
#   saveFileFigures = "./Figures/AssociationPlots/",
#   identifier = "Gapit4PC_19",
#   colourOne = "#2ca25f", colourTwo = "#43a2ca",
#   colourThree = "#8856a7"
# )

#### Examining significant ones ####

snps <- setDT(as.data.frame(snpMatrix), keep.rownames = TRUE)
## GRVI_20170609
grvi20170609 <- rrB4_17F %>%
  select(Name, GRVI_20170609) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(GRVI_20170609 >= -log10(0.05 / nrow(rrB4_17F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-GRVI_20170609)

grvi20170609$Name

grvi20170609 <- setDT(as.data.frame(t(grvi20170609)), keep.rownames = TRUE)

grvi20170609 <- dat17 %>%
  select(Taxa, GRYLD, GRVI_20170609) %>%
  left_join(grvi20170609, by = c("Taxa" = "rn")) %>%
  rename(snp_10783 = V1, snp_10790 = V2) %>%
  pivot_longer(
    cols = snp_10783:snp_10790,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:GRVI_20170609,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

grvi20170609_plot <- ggplot(
  grvi20170609,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

grvi20170609_plot
ggpubr::ggexport(grvi20170609_plot,
                 filename = "./Figures/AssociationPlots/grvi20170609_effect.png",
                 width = 3000, height = 1250, res = 360)

## NDRE_20170609
ndre20170609 <- rrB4_17F %>%
  select(Name, NDRE_20170609) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(NDRE_20170609 >= -log10(0.05 / nrow(rrB4_17F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-NDRE_20170609)

ndre20170609$Name

ndre20170609 <- setDT(as.data.frame(t(ndre20170609)), keep.rownames = TRUE)

ndre20170609 <- dat17 %>%
  select(Taxa, GRYLD, NDRE_20170609) %>%
  left_join(ndre20170609, by = c("Taxa" = "rn")) %>%
  rename(snp_10783 = V1, snp_10790 = V2) %>%
  pivot_longer(
    cols = snp_10783:snp_10790,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:NDRE_20170609,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

ndre20170609_plot <- ggplot(
  ndre20170609,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

ndre20170609_plot

ggpubr::ggexport(ndre20170609_plot,
                 filename = "./Figures/AssociationPlots/ndre20170609_effect.png",
                 width = 3000, height = 1250, res = 360)

## NIR_20170503
nir20170503 <- rrB4_17F %>%
  select(Name, NIR_20170503) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(NIR_20170503 >= -log10(0.05 / nrow(rrB4_17F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-NIR_20170503)

nir20170503$Name

nir20170503 <- setDT(as.data.frame(t(nir20170503)), keep.rownames = TRUE)

nir20170503 <- dat17 %>%
  select(Taxa, GRYLD, NIR_20170503) %>%
  left_join(nir20170503, by = c("Taxa" = "rn")) %>%
  rename(snp_3880 = V1) %>%
  pivot_longer(
    cols = snp_3880,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:NIR_20170503,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

nir20170503_plot <- ggplot(
  nir20170503,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

nir20170503_plot
ggpubr::ggexport(nir20170503_plot,
                 filename = "./Figures/AssociationPlots/nir20170503_effect.png",
                 width = 3000, height = 1250, res = 360)

## GNDVI_20180613
gndvi20180613 <- rrB4_18F %>%
  select(Name, GNDVI_20180613) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(GNDVI_20180613 >= -log10(0.05 / nrow(rrB4_18F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-GNDVI_20180613)

gndvi20180613$Name

gndvi20180613 <- setDT(as.data.frame(t(gndvi20180613)), keep.rownames = TRUE)

gndvi20180613 <- dat18 %>%
  select(Taxa, GRYLD, GNDVI_20180613) %>%
  left_join(gndvi20180613, by = c("Taxa" = "rn")) %>%
  rename(
    snp_1197 = V1,
    snp_1198 = V2,
    snp_4200 = V3
  ) %>%
  pivot_longer(
    cols = snp_1197:snp_4200,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:GNDVI_20180613,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

gndvi20180613_plot <- ggplot(
  gndvi20180613,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

gndvi20180613_plot

ggpubr::ggexport(gndvi20180613_plot,
                 filename = "./Figures/AssociationPlots/gndvi20180613_effect.png",
                 width = 3800, height = 1250, res = 360)

## GNDVI_20180613
gryld2018 <- rrB4_18F %>%
  select(Name, GRYLD_2018) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(GRYLD_2018 >= -log10(0.05 / nrow(rrB4_18F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-GRYLD_2018)

gryld2018$Name

gryld2018 <- setDT(as.data.frame(t(gryld2018)), keep.rownames = TRUE)

gryld2018 <- dat18 %>%
  select(Taxa, GRYLD) %>%
  left_join(gryld2018, by = c("Taxa" = "rn")) %>%
  rename(
    snp_8026 = V1,
    snp_8043 = V2,
    snp_8057 = V3,
    snp_8064 = V4,
    snp_8065 = V5,
    snp_8406 = V6
  ) %>%
  pivot_longer(
    cols = snp_8026:snp_8406,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

gryld2018_plot <- ggplot(
  gryld2018,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  ggpubr::stat_compare_means(method = "t.test")

gryld2018_plot

ggpubr::ggexport(gryld2018_plot,
                 filename = "./Figures/AssociationPlots/gryld2018_effect.png",
                 width = 3750, height = 1250, res = 360)

## RE_20180423
re20180423 <- rrB4_18F %>%
  select(Name, RE_20180423) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(RE_20180423 >= -log10(0.05 / nrow(rrB4_18F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-RE_20180423)

re20180423$Name

re20180423 <- setDT(as.data.frame(t(re20180423)), keep.rownames = TRUE)

re20180423 <- dat18 %>%
  select(Taxa, GRYLD, RE_20180423) %>%
  left_join(re20180423, by = c("Taxa" = "rn")) %>%
  rename(
    snp_12550 = V1
  ) %>%
  pivot_longer(
    cols = snp_12550,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:RE_20180423,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

re20180423_plot <- ggplot(
  re20180423,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

re20180423_plot
ggpubr::ggexport(re20180423_plot,
                 filename = "./Figures/AssociationPlots/re20180423_effect.png",
                 width = 3000, height = 1250, res = 360)

## RE_20180606
re20180606 <- rrB4_18F %>%
  select(Name, RE_20180606) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(RE_20180606 >= -log10(0.05 / nrow(rrB4_18F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-RE_20180606)

re20180606$Name

re20180606 <- setDT(as.data.frame(t(re20180606)), keep.rownames = TRUE)

re20180606 <- dat18 %>%
  select(Taxa, GRYLD, RE_20180606) %>%
  left_join(re20180606, by = c("Taxa" = "rn")) %>%
  rename(
    snp_4960 = V1,
    snp_4961 = V2
  ) %>%
  pivot_longer(
    cols = snp_4960:snp_4961,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:RE_20180606,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

re20180606_plot <- ggplot(
  re20180606,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

re20180606_plot
ggpubr::ggexport(re20180606_plot,
                 filename = "./Figures/AssociationPlots/re20180606_effect.png",
                 width = 3000, height = 1250, res = 360)

## NDRE_20190624
ndre20190624 <- rrB4_19F %>%
  select(Name, NDRE_20190624) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(NDRE_20190624 >= -log10(0.05 / nrow(rrB4_19F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-NDRE_20190624)

ndre20190624$Name

ndre20190624 <- setDT(as.data.frame(t(ndre20190624)), keep.rownames = TRUE)

ndre20190624 <- dat19 %>%
  select(Taxa, GRYLD, NDRE_20190624) %>%
  left_join(ndre20190624, by = c("Taxa" = "rn")) %>%
  rename(
    snp_4751 = V1,
    snp_4754 = V2
  ) %>%
  pivot_longer(
    cols = snp_4751:snp_4754,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:NDRE_20190624,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

ndre20190624_plot <- ggplot(
  ndre20190624,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

ndre20190624_plot
ggpubr::ggexport(ndre20190624_plot,
                 filename = "./Figures/AssociationPlots/ndre20190624_effect.png",
                 width = 3250, height = 1250, res = 360)

## NDVI_20190424
ndvi20190424 <- rrB4_19F %>%
  select(Name, NDVI_20190424) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(NDVI_20190424 >= -log10(0.05 / nrow(rrB4_19F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-NDVI_20190424)

ndvi20190424$Name

ndvi20190424 <- setDT(as.data.frame(t(ndvi20190424)), keep.rownames = TRUE)

ndvi20190424 <- dat19 %>%
  select(Taxa, GRYLD, NDVI_20190424) %>%
  left_join(ndvi20190424, by = c("Taxa" = "rn")) %>%
  rename(
    snp_9109 = V1
  ) %>%
  pivot_longer(
    cols = snp_9109,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:NDVI_20190424,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

ndvi20190424_plot <- ggplot(
  ndvi20190424,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

ndvi20190424_plot

ggpubr::ggexport(ndvi20190424_plot,
                 filename = "./Figures/AssociationPlots/ndvi20190424_effect.png",
                 width = 3000, height = 1250, res = 360)

## NDVI_20190624
ndvi20190624 <- rrB4_19F %>%
  select(Name, NDVI_20190624) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(NDVI_20190624 >= -log10(0.05 / nrow(rrB4_19F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-NDVI_20190624)

ndvi20190624$Name

ndvi20190624 <- setDT(as.data.frame(t(ndvi20190624)), keep.rownames = TRUE)

ndvi20190624 <- dat19 %>%
  select(Taxa, GRYLD, NDVI_20190624) %>%
  left_join(ndvi20190624, by = c("Taxa" = "rn")) %>%
  rename(
    snp_4751 = V1,
    snp_4754 = V2
  ) %>%
  pivot_longer(
    cols = snp_4751:snp_4754,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:NDVI_20190624,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

ndvi20190624_plot <- ggplot(
  ndvi20190624,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

ndvi20190624_plot
ggpubr::ggexport(ndvi20190624_plot,
                 filename = "./Figures/AssociationPlots/ndvi20190624_effect.png",
                 width = 3250, height = 1250, res = 360)

## NIR_20190530
nir20190530 <- rrB4_19F %>%
  select(Name, Nir_20190530) %>%
  tidylog::mutate(Name = str_remove(Name, "V"),
                  Name = as.numeric(Name)) %>%
  filter(Nir_20190530 >= -log10(0.05 / nrow(rrB4_19F))) %>%
  left_join(snpChip, by = c("Name" = "rs_number")) %>%
  select(-Nir_20190530)

nir20190530$Name

nir20190530 <- setDT(as.data.frame(t(nir20190530)), keep.rownames = TRUE)

nir20190530 <- dat19 %>%
  select(Taxa, GRYLD, Nir_20190530) %>%
  left_join(nir20190530, by = c("Taxa" = "rn")) %>%
  rename(
    snp_8026 = V1,
    snp_8043 = V2,
    snp_8057 = V3,
    snp_8064 = V4,
    snp_8065 = V5,
    snp_8406 = V6
  ) %>%
  pivot_longer(
    cols = snp_8026:snp_8406,
    names_to = "SNP",
    values_to = "State"
  ) %>%
  tidylog::mutate(State = as.factor(State)) %>%
  pivot_longer(
    cols = GRYLD:Nir_20190530,
    names_to = "phenotype",
    values_to = "phenotype_value"
  )

nir20190530_plot <- ggplot(
  nir20190530,
  aes(
    x = SNP,
    y = phenotype_value,
    colour = State
  )
) +
  geom_boxplot() +
  facet_wrap(~phenotype, scales = "free") +
  ggpubr::stat_compare_means(method = "t.test")

nir20190530_plot
ggpubr::ggexport(nir20190530_plot,
                 filename = "./Figures/AssociationPlots/nir20190530_effect.png",
                 width = 7000, height = 1250, res = 360)




