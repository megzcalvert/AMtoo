rm(list = objects())
ls()

library(MatrixModels)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(asreml)
library(beepr)
library(rrBLUP)
library(Hmisc)
library(broom)
library(MVN)
library(MASS)
library(car)
library(ape)
library(ClassDiscovery)
library(ggdendro)
library(pvclust)
library(viridis)
library(RColorBrewer)

##### Set up work space ####
#### Theme set
custom_theme <- theme_minimal() %+replace%
  theme(
    axis.title = element_text(
      colour = "black",
      size = rel(2)
    ),
    axis.title.x = element_text(
      vjust = 0,
      margin = margin(
        t = 0, r = 0.25,
        b = 0, l = 0,
        unit = "cm"
      )
    ),
    axis.title.y = element_text(
      vjust = 1,
      angle = 90,
      margin = margin(
        t = 0, r = 0.25,
        b = 0, l = 0.1,
        unit = "cm"
      )
    ),
    axis.text = element_text(
      colour = "black",
      size = rel(1.5)
    ),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(3, "pt"),
    axis.line = element_line(
      color = "black",
      size = 0.5
    ),
    legend.key.size = unit(1.5, "lines"),
    legend.margin = margin(
      t = 0, r = 0.75,
      b = 0, l = 0.75,
      unit = "cm"
    ),
    legend.text = element_text(size = rel(2)),
    legend.title = element_text(size = rel(1.5)),
    panel.grid.major = element_line(
      colour = "#969696",
      linetype = 3
    ),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(
      size = rel(2),
      margin = margin(
        t = 0.1, r = 0.1,
        b = 0.1, l = 0.1,
        unit = "cm"
      )
    ),
    plot.margin = margin(
      t = 0.5, r = 0.5,
      b = 0.5, l = 0,
      unit = "cm"
    ),
    plot.title = element_text(
      colour = "black",
      size = rel(2.25),
      vjust = 0,
      hjust = 0,
      margin = margin(
        t = 0.25, r = 0.25,
        b = 0.5, l = 0.25,
        unit = "cm"
      )
    ),
    plot.subtitle = element_text(
      colour = "black",
      size = rel(1.75)
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 1
    ),
    strip.text = element_text(
      colour = "black",
      size = rel(1.25)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")
set.seed(1642)

#### Load data ####
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

snpChip[snpChip == snpChip$allele_a] <- 1
snpChip[snpChip == snpChip$allele_b] <- -1
snpChip[snpChip == "H"] <- 0
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

snpMatrix <- setDT(as.data.frame(snpMatrix), keep.rownames = T) %>% 
  arrange(rn)
snpMatrix[1:10, 1:10]

##### Phenotypes ####

pheno17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoLong <- fread("./Phenotype_Database/Pheno_Long171819.txt")
hddt17<- fread("./Phenotype_Database/HDDT2017.txt")
hddt18<- fread("./Phenotype_Database/HDDT2018.txt")
hddt19<- fread("./Phenotype_Database/HDDT2019.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(pheno19)
glimpse(phenoLong)

phenoLong <- phenoLong %>%
  dplyr::rename(Plot_ID = entity_id) %>%
  tidylog::select(Plot_ID, Variety, block, rep, range, column) %>% 
  arrange(Variety)

pheno17 <- pheno17 %>%
  mutate(
    Date = as.Date(Date, format = "%Y-%m-%d"),
    Date = format(Date, "%Y%m%d")
  ) %>%
  unite("ID", c("ID", "Date")) %>%
  spread(key = ID, value = value) %>%
  tidylog::select(
    Plot_ID, Variety, GRYLD,
    GNDVI_20170331:RedEdge_20170609
  ) %>%
  tidylog::inner_join(phenoLong) %>%
  tidylog::inner_join(hddt17, by = c("Plot_ID" = "plots17")) %>% 
  tidylog::select(-hddt17) %>% 
  distinct() %>%
  mutate(
    Plot_ID = as.factor(Plot_ID),
    Variety = as.factor(Variety),
    block = as.factor(block),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  ) %>% 
  arrange(Variety)

pheno18 <- pheno18 %>%
  mutate(
    Date = as.Date(Date, format = "%Y-%m-%d"),
    Date = format(Date, "%Y%m%d")
  ) %>%
  dplyr::rename(Plot_ID = entity_id) %>%
  unite("ID", c("ID", "Date")) %>%
  spread(key = ID, value = value) %>%
  tidylog::select(
    Plot_ID, Variety, GRYLD,
    GNDVI_20171120:RE_20180613
  ) %>%
  tidylog::inner_join(phenoLong) %>%
  tidylog::inner_join(hddt18, by = c("Plot_ID" = "plots18")) %>% 
  tidylog::select(-starts_with("height_"),
                  -hddt18) %>%
  distinct() %>%
  mutate(
    Plot_ID = as.factor(Plot_ID),
    Variety = as.factor(Variety),
    block = as.factor(block),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  ) %>% 
  arrange(Variety)

pheno19 <- pheno19 %>%
  mutate(
    Date = as.Date(Date, format = "%Y-%m-%d"),
    Date = format(Date, "%Y%m%d")
  ) %>%
  dplyr::rename(Plot_ID = entity_id) %>%
  unite("ID", c("ID", "Date")) %>%
  spread(key = ID, value = value) %>%
  tidylog::select(
    Plot_ID, Variety, GRYLD,
    GNDVI_20190103:RE_20190624
  ) %>%
  tidylog::inner_join(phenoLong) %>%
  tidylog::inner_join(hddt19, by = c("Plot_ID" = "plots19")) %>% 
  tidylog::select(-starts_with("height_"),
                  -hddt19) %>%
  distinct() %>%
  mutate(
    Plot_ID = as.factor(Plot_ID),
    Variety = as.factor(Variety),
    block = as.factor(block),
    rep = as.factor(rep),
    range = as.factor(range),
    column = as.factor(column)
  ) %>% 
  arrange(Variety)

asreml.license.status()

##### Trialing something ####

set.seed(1962)

pheno17 <- pheno17 %>%
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>% 
  arrange(Variety)

snpMatrix17 <- snpMatrix %>%
  semi_join(pheno17, by = c("rn" = "Variety")) %>% 
  arrange(rn)
snpMatrix17[1:5, 1:5]

snpLines <- snpMatrix17$rn
snpMatrix17 <- snpMatrix17 %>%
  column_to_rownames(var = "rn")

dim(pheno17)
dim(snpMatrix17)

mean(pheno17$GRYLD)

dat<- pheno17 %>% 
  tidylog::select(Variety,rep,block,GRYLD)

t17 <- asreml(
  fixed = GRYLD ~ 0 + Variety,
  random = ~ rep + rep:block,
  data = dat
)
summary(t17)
plot(t17)
blues <- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "Variety_")
dat17 <- blues %>%
  dplyr::rename(GRYLD = effect) %>%
  glimpse()

##### Making it a function for all of the VI's 2017

calculatingBlues <- function(dat, saveFile, joinFile, ...) {
  drop <- c(
    "block", "rep", "Variety", "year",
    "column", "range", "Plot_ID", "GRYLD"
  )
  effectvars <- dat %>% 
    dplyr::select(!any_of(drop))
  
  traits <- colnames(effectvars)
  traits
  fieldInfo <- dat %>%
    tidylog::select(Variety, rep, block, column, range)
  
  for (i in traits) {
    print(paste("Working on trait", i))
    j <- i
    
    traitData<- dat %>% 
      pull(i)
    
    data <- cbind(fieldInfo, traitData)
    names(data) <- c("Variety", "rep", "block", "column", "range", "Trait")
    print(colnames(data))
    
    t17 <- asreml(
      fixed = Trait ~ 0 + Variety,
      random = ~ rep + rep:block,
      data = data
    )
    pdf(paste0(
      saveFile,
      i, ".pdf"
    ))
    plot(t17)
    
    blues <- setDT(as.data.frame(coef(t17)$fixed), keep.rownames = T)
    blues$rn <- str_remove(blues$rn, "Variety_")
    colnames(blues)[colnames(blues) == "effect"] <- paste(i)
    joinFile <- blues %>%
      inner_join(joinFile)
    graphics.off()
  }
  return(joinFile)
}

blues17 <- calculatingBlues(
  dat = pheno17,
  saveFile = "./Figures/AsremlPlots/ASREML_Blues17_",
  joinFile = dat17
)

##### Making it a function for all of the VI's 2018
dat<- pheno18 %>% 
  tidylog::select(GRYLD, Variety, rep, block)

t18 <- asreml(
  fixed = GRYLD ~ 0 + Variety,
  random = ~ rep + rep:block,
  data = dat
)
plot(t18)
blues <- setDT(as.data.frame(coef(t18)$fixed), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "Variety_")
dat18 <- blues %>%
  tidylog::rename(GRYLD = effect) %>%
  glimpse()

blues18 <- calculatingBlues(
  dat = pheno18,
  saveFile = "./Figures/AsremlPlots/ASREML_Blues18_",
  joinFile = dat18
)

##### Making it a function for all of the VI's 2019
dat<- pheno19 %>% 
  tidylog::select(GRYLD, Variety, rep, block)

t19 <- asreml(
  fixed = GRYLD ~ 0 + Variety,
  random = ~ rep + rep:block,
  data = dat
)
summary(t19)
plot(t19)
blues <- setDT(as.data.frame(coef(t19)$fixed), keep.rownames = T)
blues$rn <- str_remove(blues$rn, "Variety_")
dat19 <- blues %>%
  tidylog::rename(GRYLD = effect) %>%
  glimpse()

blues19 <- calculatingBlues(
  dat = pheno19,
  saveFile = "./Figures/AsremlPlots/ASREML_Blues19_",
  joinFile = dat19
)

write.table(blues17, "./Phenotype_Database/Hyp10BLUEs_17.txt",
            quote = F,
            sep = "\t", row.names = FALSE, col.names = T
)
write.table(blues18, "./Phenotype_Database/Hyp10BLUEs_18.txt",
            quote = F,
            sep = "\t", row.names = FALSE, col.names = T
)
write.table(blues19, "./Phenotype_Database/Hyp10BLUEs_19.txt",
            quote = F,
            sep = "\t", row.names = FALSE, col.names = T
)
dat17 <- fread("./Phenotype_Database/Hyp10BLUEs_17.txt", header = T)
dat18 <- fread("./Phenotype_Database/Hyp10BLUEs_18.txt", header = T)
dat19 <- fread("./Phenotype_Database/Hyp10BLUEs_19.txt", header = T)

##### Select in 2016/2017 and see what happens in 2017/2018 ####
# Using NDVI_20170512 for the selection
dat17 %>%
  ggplot(aes(x = NDVI_20170512)) +
  geom_histogram(colour = "black", fill = "white") +
  geom_vline(xintercept = mean(dat17$NDVI_20170512), linetype = 2)

ndvi0512_top <- dat17 %>%
  dplyr::arrange(desc(NDVI_20170512)) %>% # use desc() if you want highest to lowest
  tidylog::select(rn, GRYLD, NDVI_20170512) %>%
  glimpse()

0.2 * nrow(dat17)
ndvi0512_top <- ndvi0512_top[1:round(0.2 * nrow(dat17)), ]

dat17 %>%
  ggplot(aes(x = NDVI_20170512)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_vline(xintercept = mean(dat17$NDVI_20170512), linetype = 2) +
  geom_vline(
    xintercept = mean(ndvi0512_top$NDVI_20170512), colour = "blue",
    linetype = 2
  ) +
  geom_vline(xintercept = min(ndvi0512_top$NDVI_20170512), colour = "blue") +
  labs(
    title = "Distribution of NDVI_20170512",
    subtitle = "All lines - black, top 20% of NDVI_20170512 - blue"
  )

dat17 %>%
  ggplot(aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_histogram(
    data = ndvi0512_top, aes(x = GRYLD), colour = "blue",
    bins = 100
  ) +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2) +
  geom_vline(
    xintercept = mean(ndvi0512_top$GRYLD), linetype = 2,
    colour = "blue"
  ) +
  labs(
    title = "Distribution of GRYLD in 2016/2017",
    subtitle = "All lines - black, lines selected from NDVI_20170512 - blue"
  )

ndviSelect <- ndvi0512_top %>%
  left_join(dat18, by = "rn") %>%
  left_join(dat19, by = "rn") %>%
  tidylog::select(rn, GRYLD.x, NDVI_20170512, GRYLD.y, GRYLD) %>%
  glimpse()

dat18 %>%
  ggplot(aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_histogram(
    data = ndviSelect, aes(x = GRYLD.y), colour = "blue",
    bins = 100
  ) +
  geom_vline(xintercept = mean(dat18$GRYLD), linetype = 2) +
  geom_vline(
    xintercept = mean(ndviSelect$GRYLD.y), linetype = 2,
    colour = "blue"
  ) +
  geom_vline(
    xintercept = mean(dat17$GRYLD), linetype = 2,
    colour = "#737373"
  ) +
  labs(
    title = "Distribution of GRYLD in 2017/2018",
    subtitle = "All lines - black, lines selected from NDVI_20170512 - blue, mean from 2016/2017 - grey"
  )

dat19 %>%
  ggplot(aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_histogram(
    data = ndviSelect, aes(x = GRYLD), colour = "blue",
    bins = 100
  ) +
  geom_vline(xintercept = mean(dat19$GRYLD), linetype = 2) +
  geom_vline(
    xintercept = mean(ndviSelect$GRYLD), linetype = 2,
    colour = "blue"
  ) +
  geom_vline(
    xintercept = mean(dat17$GRYLD), linetype = 2,
    colour = "#737373"
  ) +
  labs(
    title = "Distribution of GRYLD in 2018/2019",
    subtitle = "All lines - black, lines selected from NDVI_20170512 - blue, mean from 2016/2017 - grey"
  )

mean(ndviSelect$GRYLD.x)
mean(ndviSelect$GRYLD.y)
mean(ndviSelect$GRYLD)
mean(dat17$GRYLD)
mean(dat18$GRYLD)
mean(dat19$GRYLD)

t.test(dat18$GRYLD, ndviSelect$GRYLD.y)
t.test(dat18$GRYLD, dat17$GRYLD)
t.test(dat19$GRYLD, ndviSelect$GRYLD)

# Check RedEdge because it's opposite
dat17 %>%
  ggplot(aes(x = RedEdge_20170512)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_vline(xintercept = mean(dat17$RedEdge_20170512), linetype = 2)

RE0512_top <- dat17 %>%
  dplyr::arrange(RedEdge_20170512) %>% # use desc() if you want highest to lowest
  tidylog::select(rn, GRYLD, RedEdge_20170512) %>%
  glimpse()
0.2 * nrow(dat17)
RE0512_top <- RE0512_top[1:round(0.2 * nrow(dat17)), ]

dat17 %>%
  ggplot(aes(x = RedEdge_20170512)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_vline(xintercept = mean(dat17$RedEdge_20170512), linetype = 2) +
  geom_vline(
    xintercept = mean(RE0512_top$RedEdge_20170512), colour = "blue",
    linetype = 2
  ) +
  geom_vline(xintercept = max(RE0512_top$RedEdge_20170512), colour = "blue") +
  labs(
    title = "Distribution of RedEdge_20170512",
    subtitle = "All lines - black, lowest 5% of RedEdge_20170512 - blue"
  )

dat17 %>%
  ggplot(aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_histogram(
    data = RE0512_top, aes(x = GRYLD), colour = "blue",
    bins = 100
  ) +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2) +
  geom_vline(
    xintercept = mean(RE0512_top$GRYLD), linetype = 2,
    colour = "blue"
  ) +
  labs(
    title = "Distribution of GRYLD in 2016/2017",
    subtitle = "All lines - black, lines selected from RedEdge_20170512 - blue"
  )

reSelect <- RE0512_top %>%
  left_join(dat18, by = "rn") %>%
  left_join(dat19, by = "rn") %>%
  tidylog::select(rn, GRYLD.x, RedEdge_20170512, GRYLD.y, GRYLD) %>%
  glimpse()

dat18 %>%
  ggplot(aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_histogram(
    data = reSelect, aes(x = GRYLD.y), colour = "blue",
    bins = 100
  ) +
  geom_vline(xintercept = mean(dat18$GRYLD), linetype = 2) +
  geom_vline(
    xintercept = mean(reSelect$GRYLD.y), linetype = 2,
    colour = "blue"
  ) +
  geom_vline(
    xintercept = mean(dat17$GRYLD), linetype = 2,
    colour = "#737373"
  ) +
  labs(
    title = "Distribution of GRYLD in 2017/2018",
    subtitle = "All lines - black, lines selected from RedEdge_20170512 - blue, mean from 2016/2017 - grey"
  )

dat19 %>%
  ggplot(aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white", bins = 100) +
  geom_histogram(
    data = reSelect, aes(x = GRYLD), colour = "blue",
    bins = 100
  ) +
  geom_vline(xintercept = mean(dat19$GRYLD), linetype = 2) +
  geom_vline(
    xintercept = mean(reSelect$GRYLD), linetype = 2,
    colour = "blue"
  ) +
  geom_vline(
    xintercept = mean(dat17$GRYLD), linetype = 2,
    colour = "#737373"
  ) +
  labs(
    title = "Distribution of GRYLD in 2018/2019",
    subtitle = "All lines - black, lines selected from RedEdge_20170512 - blue, mean from 2016/2017 - grey"
  )

mean(reSelect$GRYLD.x)
mean(reSelect$GRYLD.y)
mean(reSelect$GRYLD)
t.test(dat18$GRYLD, reSelect$GRYLD.y)
t.test(dat19$GRYLD, reSelect$GRYLD)
t.test(reSelect$GRYLD.y, ndviSelect$GRYLD.y)
t.test(reSelect$GRYLD.x, ndviSelect$GRYLD.x)

###############################################################
####                    rrBlup trial                      ####

snpMatrix[1:5, 1:5]

dat17 <- dat17 %>%
  semi_join(snpMatrix, by = "rn") %>% 
  arrange(rn)

snpMatrix17 <- snpMatrix %>%
  semi_join(dat17, by = "rn") %>% 
  arrange(rn) %>%
  column_to_rownames(var = "rn")

all(dat17$rn == rownames(snpMatrix17))

snpMatrix17 <- as.matrix(snpMatrix17)

dat18 <- dat18 %>%
  semi_join(snpMatrix, by = "rn") %>% 
  arrange(rn)

snpMatrix18 <- snpMatrix %>%
  semi_join(dat18, by = "rn")  %>% 
  arrange(rn) %>%
  column_to_rownames(var = "rn")

snpMatrix18 <- as.matrix(snpMatrix18)

dat19 <- dat19 %>%
  semi_join(snpMatrix, by = "rn") %>% 
  arrange(rn)

snpMatrix19 <- snpMatrix %>%
  semi_join(dat19, by = "rn")  %>% 
  arrange(rn) %>%
  column_to_rownames(var = "rn")

snpMatrix19 <- as.matrix(snpMatrix19)

# Predict marker effects

##### Determing marker effects for each trait 2017
gryldME <- mixed.solve(dat17$GRYLD, Z = snpMatrix17)
ndre14MayME <- mixed.solve(dat17$NDVI_20170512, Z = snpMatrix17)
tidy(cor.test(gryldME$u, ndre14MayME$u))
traitME_17 <- as.data.frame(gryldME$u)

colnames(traitME_17) <- "GRYLD"

traits <- dat17 %>%
  tidylog::select(-rn, -GRYLD) %>%
  colnames()

for (i in traits) {
  print(paste("Working on trait", i))
  y <- dat17[[i]]
  y
  meRes <- mixed.solve(y = y, Z = snpMatrix17)
  traitME_17[[i]] <- meRes$u
}

##### Determing marker effects for each trait 2018

gryldME <- mixed.solve(dat18$GRYLD, Z = snpMatrix18)
re20180613ME <- mixed.solve(dat18$RE_20180613, Z = snpMatrix18)
tidy(cor.test(gryldME$u, re20180613ME$u))

traitME_18 <- as.data.frame(gryldME$u)
colnames(traitME_18) <- "GRYLD"

traits <- dat18 %>%
  tidylog::select(-rn, -GRYLD) %>%
  colnames()

for (i in traits) {
  print(paste("Working on trait", i))
  y <- dat18[[i]]
  y
  meRes <- mixed.solve(y = y, Z = snpMatrix18)
  traitME_18[[i]] <- meRes$u
}

dev.list()
graphics.off()

##### Determing marker effects for each trait 2019

gryldME <- mixed.solve(dat19$GRYLD, Z = snpMatrix19)
re20190624ME <- mixed.solve(dat19$RE_20190624, Z = snpMatrix19)
tidy(cor.test(gryldME$u, re20190624ME$u))

traitME_19 <- as.data.frame(gryldME$u)
colnames(traitME_19) <- "GRYLD"

traits <- dat19 %>%
  tidylog::select(-rn, -GRYLD) %>%
  colnames()

for (i in traits) {
  print(paste("Working on trait", i))
  y <- dat19[[i]]
  y
  meRes <- mixed.solve(y = y, Z = snpMatrix19)
  traitME_19[[i]] <- meRes$u
}

dev.list()
graphics.off()

##### Looking at the marker effects distributions
# Trying something....

ggplot(data = pheno17, aes(x = GRYLD)) +
  geom_density(colour = "#008856") +
  geom_vline(
    xintercept = mean(pheno17$GRYLD), linetype = 2,
    colour = "#008856"
  ) +
  geom_density(data = pheno18, aes(x = GRYLD), colour = "#0067a5") +
  geom_vline(
    xintercept = mean(pheno18$GRYLD), linetype = 2,
    colour = "#0067a5"
  ) +
  geom_density(data = pheno19, aes(x = GRYLD), colour = "#604e97") +
  geom_vline(
    xintercept = mean(pheno19$GRYLD), linetype = 2,
    colour = "#604e97"
  ) +
  labs(
    title = "Distribution of GRYLD over 2016/2017, 2017/2018 and 2018/2019",
    subtitle = "2016/2017-green, 2017/2018-blue, 2018/2019-violet"
  )

ggplot(data = dat17, aes(x = GRYLD)) +
  geom_density(colour = "#008856") +
  geom_vline(xintercept = mean(dat17$GRYLD), linetype = 2, colour = "#008856") +
  geom_density(data = dat18, aes(x = GRYLD), colour = "#0067a5") +
  geom_vline(xintercept = mean(dat18$GRYLD), linetype = 2, colour = "#0067a5") +
  geom_density(data = dat19, aes(x = GRYLD), colour = "#604e97") +
  geom_vline(xintercept = mean(dat19$GRYLD), linetype = 2, colour = "#604e97") +
  labs(
    title = "Distribution of GRYLD BLUEs over 2016/2017, 2017/2018 and 2018/2019",
    subtitle = "2016/2017-green, 2017/2018-blue, 2018/2019-violet"
  )

ggplot(data = traitME_17, aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white") +
  labs(title = "Marker Effect distribution for GRYLD2016/2017")

ggplot(data = traitME_18, aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white") +
  labs(title = "Marker Effect distribution for GRYLD 2017/2018")

ggplot(data = traitME_19, aes(x = GRYLD)) +
  geom_histogram(colour = "black", fill = "white") +
  labs(title = "Marker Effect distribution for GRYLD 2018/2019")

ggplot(data = traitME_17, aes(x = GRYLD)) +
  geom_histogram(
    colour = "#008856",
    fill = NA,
    bins = 100
  ) +
  geom_histogram(
    data = traitME_18, aes(x = GRYLD), colour = "#0067a5",
    fill = NA, 
    bins = 100
  ) +
  geom_vline(xintercept = mean(traitME_17$GRYLD), colour = "#008856") +
  geom_vline(xintercept = mean(traitME_18$GRYLD), colour = "#0067a5") +
  geom_histogram(
    data = traitME_19, aes(x = GRYLD), colour = "#604e97",
    fill = NA, bins = 100
  ) +
  geom_vline(xintercept = mean(traitME_19$GRYLD), colour = "#604e97") +
  labs(
    title = "Marker Effect distribution for GRYLD 2016/2017, 2017/2018 and 2018/2019",
    subtitle = "2016/2017-green, 2017/2018-blue, 2018/2019-violet",
    y = "Frequency"
  )

write.table(traitME_17, "./Genotype_Database/traitMarkerEffects17.txt",
            quote = F, sep = "\t", row.names = F, col.names = T
)
write.table(traitME_18, "./Genotype_Database/traitMarkerEffects18.txt",
            quote = F, sep = "\t", row.names = F, col.names = T
)
write.table(traitME_19, "./Genotype_Database/traitMarkerEffects19.txt",
            quote = F, sep = "\t", row.names = F, col.names = T
)

traitME_17<- fread("./Genotype_Database/traitMarkerEffects17.txt")
traitME_18<- fread("./Genotype_Database/traitMarkerEffects18.txt")
traitME_19<- fread("./Genotype_Database/traitMarkerEffects19.txt")

mean(pheno17$GRYLD)
mean(pheno18$GRYLD)
mean(pheno19$GRYLD)
var(pheno17$GRYLD)
var(pheno18$GRYLD)
var(pheno19$GRYLD)

mean(dat17$GRYLD)
mean(dat18$GRYLD)
mean(dat19$GRYLD)
var(dat17$GRYLD)
var(dat18$GRYLD)
var(dat19$GRYLD)

mean(traitME_17$GRYLD)
mean(traitME_18$GRYLD)
mean(traitME_19$GRYLD)
var(traitME_17$GRYLD)
var(traitME_18$GRYLD)
var(traitME_19$GRYLD)

t.test(traitME_17$GRYLD, traitME_18$GRYLD)
t.test(traitME_17$GRYLD, traitME_19$GRYLD)
t.test(traitME_18$GRYLD, traitME_19$GRYLD)
t.test(pheno17$GRYLD, pheno18$GRYLD)
t.test(pheno17$GRYLD, pheno19$GRYLD)
t.test(pheno18$GRYLD, pheno19$GRYLD)
t.test(dat17$GRYLD, dat18$GRYLD)
t.test(dat17$GRYLD, dat19$GRYLD)
t.test(dat18$GRYLD, dat19$GRYLD)

ggplot() +
  geom_point(aes(x = mean(pheno17$GRYLD), y = mean(traitME_17$GRYLD)),
             colour = "#008856"
  ) +
  geom_point(aes(x = mean(pheno18$GRYLD), y = mean(traitME_18$GRYLD)),
             colour = "#0067a5"
  ) +
  geom_point(aes(x = mean(pheno19$GRYLD), y = mean(traitME_19$GRYLD)),
             colour = "#604e97"
  ) +
  labs(
    x = "mean GRYLD", y = "mean marker effect",
    title = "Comparison of mean GRYLD and mean marker effect",
    subtitle = "2016/2017-green, 2017/2018-blue, 2018/2019-violet"
  )

ggplot() +
  geom_point(aes(x = mean(pheno17$GRYLD), y = var(traitME_17$GRYLD)),
             colour = "#008856"
  ) +
  geom_point(aes(x = mean(pheno18$GRYLD), y = var(traitME_18$GRYLD)),
             colour = "#0067a5"
  ) +
  geom_point(aes(x = mean(pheno19$GRYLD), y = var(traitME_19$GRYLD)),
             colour = "#604e97"
  ) +
  labs(
    x = "mean GRYLD", y = "variance of marker effect",
    title = "Comparison of mean GRYLD and variance of marker effect",
    subtitle = "2016/2017-green, 2017/2018-blue, 2018/2019-violet"
  )

# Correlation between ME for GRYLD in different years
tidy(rcorr(traitME_17$GRYLD, traitME_18$GRYLD))
tidy(rcorr(traitME_17$GRYLD, traitME_19$GRYLD))
tidy(rcorr(traitME_18$GRYLD, traitME_19$GRYLD))

gryldME <- as.matrix(cbind(traitME_17$GRYLD, traitME_18$GRYLD, traitME_19$GRYLD))
colnames(gryldME) <- c("gryld17", "gryld18", "gryld19")

gryldME_chr1 <- gryldME[1:1149, ]
gryldME_chr2 <- gryldME[1150:2028, ]
gryldME_chr3 <- gryldME[2029:2511, ]
gryldME_chr4 <- gryldME[2512:3578, ]
gryldME_chr5 <- gryldME[3579:4677, ]
gryldME_chr6 <- gryldME[4678:5118, ]
gryldME_chr7 <- gryldME[5119:5969, ]
gryldME_chr8 <- gryldME[5970:6762, ]
gryldME_chr9 <- gryldME[6763:7062, ]
gryldME_chr10 <- gryldME[7063:7672, ]
gryldME_chr11 <- gryldME[7673:8096, ]
gryldME_chr12 <- gryldME[8097:8193, ]
gryldME_chr13 <- gryldME[8194:8859, ]
gryldME_chr14 <- gryldME[8860:9949, ]
gryldME_chr15 <- gryldME[9950:10508, ]
gryldME_chr16 <- gryldME[10509:11374, ]
gryldME_chr17 <- gryldME[11375:12216, ]
gryldME_chr18 <- gryldME[12217:12627, ]
gryldME_chr19 <- gryldME[12628:13526, ]
gryldME_chr20 <- gryldME[13527:14221, ]
gryldME_chr21 <- gryldME[14222:14523, ]

outCor <- zoo::rollapply(gryldME_chr1, 10,
                    by = 2,
                    function(x) c(cor(x)), by.column = F
)

outCor <- as.data.frame(outCor)
outCor <- outCor %>%
  tidylog::select(-V1, -V3, -V4) %>%
  dplyr::rename(Correlation = V2) %>%
  tidylog::mutate(Bin = row_number())

outCor %>%
  ggplot(aes(x = Bin, y = Correlation)) +
  geom_line() +
  geom_smooth() +
  labs(
    title = "Sliding window correlation of GRYLD marker effects",
    subtitle = "Chr1A window = 10, slide = 2"
  ) +
  coord_cartesian(ylim = c(-1, 1))

##### Correlation Matrix examination ####

corrMatrix_17 <- rcorr(as.matrix(traitME_17))
corrMatrix_18 <- rcorr(as.matrix(traitME_18))
corrMatrix_19 <- rcorr(as.matrix(traitME_19))

##### Distribution of Correlation Matrix

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

correlationME17 <- flattenCorrMatrix(
  corrMatrix_17$r,
  corrMatrix_17$P
)
correlationME18 <- flattenCorrMatrix(
  corrMatrix_18$r,
  corrMatrix_18$P
)
correlationME19 <- flattenCorrMatrix(
  corrMatrix_19$r,
  corrMatrix_19$P
)

phenoCorrMatrix_17 <- dat17 %>%
  tidylog::select(-rn)
phenoCorrMatrix_17 <- rcorr(as.matrix(phenoCorrMatrix_17))

correlationsPheno17 <- flattenCorrMatrix(
  phenoCorrMatrix_17$r,
  phenoCorrMatrix_17$P
)

phenoCorrMatrix_18 <- dat18 %>%
  tidylog::select(-rn)
phenoCorrMatrix_18 <- rcorr(as.matrix(phenoCorrMatrix_18))

correlationsPheno18 <- flattenCorrMatrix(
  phenoCorrMatrix_18$r,
  phenoCorrMatrix_18$P
)

phenoCorrMatrix_19 <- dat19 %>%
  tidylog::select(-rn)
phenoCorrMatrix_19 <- rcorr(as.matrix(phenoCorrMatrix_19))

correlationsPheno19 <- flattenCorrMatrix(
  phenoCorrMatrix_19$r,
  phenoCorrMatrix_19$P
)

## Distribution of all correlations for VI and GRYLD ME
correlationMarkerEffects <- ggplot() +
  geom_density(
    data = correlationME17, aes(x = cor), #binwidth = 0.05,
    fill = NA, colour = "#008856"
  ) +
  geom_density(
    data = correlationME18, aes(x = cor), #binwidth = 0.05,
    fill = NA, colour = "#0067a5"
  ) +
  geom_density(
    data = correlationME19, aes(x = cor), #binwidth = 0.05,
    fill = NA, colour = "#604e97"
  ) +
  labs(
    title =
      "Distribution of the correlations between marker effects",
    subtitle =
      "All VI and GRYLD, generated by rrBLUP 2016/2017-green, 2017/2018-blue, 2018/2019-violet"
  ) +
  xlim(-1, 1)
correlationMarkerEffects

## Distribution for only those related to GRYLD
correlationsGryld17 <- correlationsPheno17 %>%
  tidylog::filter(column == "GRYLD") %>%
  tidylog::select(-column)
unique(correlationsGryld17$row)

correlationsGryld18 <- correlationsPheno18 %>%
  tidylog::filter(column == "GRYLD") %>%
  tidylog::select(-column)
unique(correlationsGryld18$row)

correlationsGryld19 <- correlationsPheno19 %>%
  tidylog::filter(column == "GRYLD") %>%
  tidylog::select(-column)
unique(correlationsGryld19$row)

correlationME17 %>%
  tidylog::filter(row == "GRYLD") %>%
  left_join(correlationsGryld17, by = c("column" = "row")) %>%
  dplyr::rename(CorToGeno = cor.x, Pvalue = p.x, corToPheno = cor.y) %>%
  tidylog::select(-p.y) %>%
  glimpse() %>%
  ggplot(aes(x = CorToGeno)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill = "white") +
  labs(
    title =
      "Distribution of the correlations between marker effects",
    subtitle =
      "GRYLD correlations with VI, generated by rrBLUP 2016/2017 Season"
  ) +
  xlim(-1, 1)

correlationME18 %>%
  tidylog::filter(row == "GRYLD") %>%
  left_join(correlationsGryld18, by = c("column" = "row")) %>%
  dplyr::rename(CorToGeno = cor.x, Pvalue = p.x, corToPheno = cor.y) %>%
  tidylog::select(-p.y) %>%
  glimpse() %>%
  ggplot(aes(x = CorToGeno)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill = "white") +
  labs(
    title =
      "Distribution of the correlations between marker effects",
    subtitle =
      "GRYLD correlations with VI, generated by rrBLUP 2017/2018 Season"
  ) +
  xlim(-1, 1)

correlationME19 %>%
  tidylog::filter(row == "GRYLD") %>%
  left_join(correlationsGryld19, by = c("column" = "row")) %>%
  dplyr::rename(CorToGeno = cor.x, Pvalue = p.x, corToPheno = cor.y) %>%
  tidylog::select(-p.y) %>%
  glimpse() %>%
  ggplot(aes(x = CorToGeno)) +
  geom_histogram(binwidth = 0.025, colour = "black", fill = "white") +
  labs(
    title =
      "Distribution of the correlations between marker effects",
    subtitle =
      "GRYLD correlations with VI, generated by rrBLUP 2018/2019 Season"
  ) +
  xlim(-1, 1)

correlationME17 %>%
  tidylog::filter(row == "GRYLD") %>%
  left_join(correlationsGryld17, by = c("column" = "row")) %>%
  dplyr::rename(CorToGeno = cor.x, Pvalue = p.x, corToPheno = cor.y) %>%
  tidylog::select(-p.y) %>%
  separate(column, c("Trait", "Date"), sep = "_") %>%
  glimpse() %>%
  ggplot(aes(x = CorToGeno, y = corToPheno, colour = Date, shape = Trait)) +
  geom_point(size = 7) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856",
    "#e68fac", "#0067a5", "#f99379", "#604e97",
    "#f6a600", "#b3446c", "#dcd300", "#882d17",
    "#8db600", "#654522", "#e25822", "#2b3d26"
  )) +
  scale_shape_manual(values = c(0, 1, 2, 8, 11, 9)) +
  theme(
    aspect.ratio = 1:1
  ) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  labs(
    title =
      "Correlation between Marker Effects correlations and Phenotypic measurements correlations",
    subtitle =
      "GRYLD correlations with VI, generated by rrBLUP 2016/2017 Season",
    x = "Correlation to GRYLD marker effects matrix",
    y = "Correlation to GRYLD phenotypic measurements"
  )

correlationME18 %>%
  tidylog::filter(row == "GRYLD") %>%
  left_join(correlationsGryld18, by = c("column" = "row")) %>%
  dplyr::rename(CorToGeno = cor.x, Pvalue = p.x, corToPheno = cor.y) %>%
  tidylog::select(-p.y) %>%
  separate(column, c("Trait", "Date"), sep = "_") %>%
  glimpse() %>%
  ggplot(aes(x = CorToGeno, y = corToPheno, colour = Date, shape = Trait)) +
  geom_point(size = 6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856",
    "#e68fac", "#0067a5", "#604e97",
    "#f6a600", "#b3446c", "#dcd300", "#882d17",
    "#8db600", "#654522", "#e25822", "#2b3d26"
  )) +
  scale_shape_manual(values = c(0, 1, 2, 8, 11, 9)) +
  theme(
    aspect.ratio = 1:1
  ) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  labs(
    title =
      "Correlation between Marker Effects correlations and Phenotypic measurements correlations",
    subtitle =
      "GRYLD correlations with VI, generated by rrBLUP 2017/2018 Season",
    x = "Correlation to GRYLD marker effects matrix",
    y = "Correlation to GRYLD phenotypic measurements"
  )

correlationME19 %>%
  tidylog::filter(row == "GRYLD") %>%
  left_join(correlationsGryld19, by = c("column" = "row")) %>%
  dplyr::rename(CorToGeno = cor.x, Pvalue = p.x, corToPheno = cor.y) %>%
  tidylog::select(-p.y) %>%
  separate(column, c("Trait", "Date"), sep = "_") %>%
  glimpse() %>%
  ggplot(aes(x = CorToGeno, y = corToPheno, colour = Date, shape = Trait)) +
  geom_point(size = 6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_color_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856",
    "#e68fac", "#0067a5", "#604e97",
    "#f6a600", "#b3446c", "#dcd300", "#882d17",
    "#8db600", "#654522", "#e25822", "#2b3d26"
  )) +
  scale_shape_manual(values = c(0, 1, 2, 8, 11, 9)) +
  theme(
    aspect.ratio = 1:1,
  ) +
  xlim(-1, 1) +
  ylim(-1, 1) +
  labs(
    title =
      "Correlation between Marker Effects correlations and Phenotypic measurements correlations",
    subtitle =
      "GRYLD correlations with VI, generated by rrBLUP 2018/2019 Season",
    x = "Correlation to GRYLD marker effects matrix",
    y = "Correlation to GRYLD phenotypic measurements"
  )

##### Converting Correlation matrix to a distance matrix ####

distmat_17 <- distanceMatrix(as.matrix(traitME_17),
                             metric = "absolute pearson"
)
distmat_18 <- distanceMatrix(as.matrix(traitME_18),
                             metric = "absolute pearson"
)
distmat_19 <- distanceMatrix(as.matrix(traitME_19),
                             metric = "absolute pearson"
)

pcoord_17 <- pcoa(distmat_17)
pcoord_18 <- pcoa(distmat_18)
pcoord_19 <- pcoa(distmat_19)
biplot.pcoa(pcoord_17)
biplot.pcoa(pcoord_18)
biplot.pcoa(pcoord_19)

## Making better biplot
pcoOrd_17 <- setDT(as.data.frame(pcoord_17$vectors), keep.rownames = T)
pcoOrd_18 <- setDT(as.data.frame(pcoord_18$vectors), keep.rownames = T)
pcoOrd_19 <- setDT(as.data.frame(pcoord_19$vectors), keep.rownames = T)

pcoOrd_17 <- pcoOrd_17 %>%
  separate(rn, c("Trait", "date"), sep = "_")
pcoOrd_17[1, 2] <- "20170613"

write.table(pcoOrd_17, "./Genotype_Database/PCOanalysisGeneticDistance_17.txt",
            quote = F, sep = "\t", row.names = F, col.names = T
)

pcoOrd_17 %>%
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0, 1, 19, 2, 5, 6, 7)) +
  scale_color_manual(values = c(
    "#875692", "#f38400",
    "#be0032", "#008856",
    "#0067a5", "#604e97",
    "#f6a600", "#b3446c", "#222222"
  )) +
  labs(
    title =
      "Principal Coordinate analysis of genetic distance matrix",
    subtitle = "2016/2017 season"
  ) +
  coord_fixed(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))

pcoOrd_18 <- pcoOrd_18 %>%
  separate(rn, c("Trait", "date"), sep = "_")
pcoOrd_18[1, 2] <- "20180615"

write.table(pcoOrd_18, "./Genotype_Database/PCOanalysisGeneticDistance_18.txt",
            quote = F, sep = "\t", row.names = F, col.names = T
)

pcoOrd_18 %>%
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0, 1, 19, 2, 5, 6, 7)) +
  scale_color_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856",
    "#e68fac", "#0067a5", "#604e97",
    "#f6a600", "#b3446c", "#dcd300", "#882d17",
    "#8db600", "#654522", "#e25822", "#2b3d26"
  )) +
  labs(
    title = "Principal Coordinate analysis of genetic distance matrix",
    subtitle = "2017/2018 season"
  ) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    aspect.ratio = 1:1
  ) +
  coord_fixed(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))

group1_18 <- pcoOrd_18 %>%
  filter(Axis.1 > 0.04 & Axis.2 > -0.19 & Axis.2 < 0.25) %>%
  dplyr::select(Trait, date) %>%
  unite("trait", c("Trait", "date"), sep = "_")

group2_18 <- pcoOrd_18 %>%
  filter(Axis.1 < 0.04 & Axis.2 > 0) %>%
  dplyr::select(Trait, date) %>%
  unite("trait", c("Trait", "date"), sep = "_")

group3_18 <- pcoOrd_18 %>%
  filter(Axis.1 < 0.1 & Axis.2 < 0) %>%
  dplyr::select(Trait, date) %>%
  unite("trait", c("Trait", "date"), sep = "_")

group1_18[group1_18 == "GRYLD_20180615"] <- "GRYLD"

group1ME_18 <- traitME_18 %>%
  dplyr::select(group1_18$trait) %>%
  as.matrix()
group2ME_18 <- traitME_18 %>%
  dplyr::select(group2_18$trait) %>%
  as.matrix()
group3ME_18 <- traitME_18 %>%
  dplyr::select(group3_18$trait) %>%
  as.matrix()

group1Cor_18 <- cor(group1ME_18)
group2Cor_18 <- cor(group2ME_18)
group3Cor_18 <- cor(group3ME_18)

pcoOrd_19 <- pcoOrd_19 %>%
  separate(rn, c("Trait", "date"), sep = "_")
pcoOrd_19[1, 2] <- "20190702"

write.table(pcoOrd_19, "./Genotype_Database/PCOanalysisGeneticDistance_19.txt",
            quote = F, sep = "\t", row.names = F, col.names = T
)

pcoOrd_19 %>%
  ggplot(aes(x = `Axis.1`, y = `Axis.2`, colour = date, shape = Trait)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(0, 19, 2, 5, 6, 7)) +
  scale_color_manual(values = c(
    "#f3c300", "#875692", "#f38400", "#a1caf1",
    "#be0032", "#848482", "#008856",
    "#e68fac", "#0067a5", "#604e97",
    "#f6a600", "#b3446c", "#dcd300", "#882d17",
    "#8db600", "#654522", "#e25822", "#2b3d26"
  )) +
  labs(
    title = "Principal Coordinate analysis of genetic distance matrix",
    subtitle = "2018/2019 season"
  ) +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    aspect.ratio = 1:1
  ) +
  coord_fixed(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))

# mantel.test(group1Cor_18, group2Cor_18, graph = T)

## Hierarchical clustering

hClustering_17 <- hclust(distmat_17, method = "ward.D2")
plot(hClustering_17, hang = -1, no.margin = T)

hClustering_18 <- hclust(distmat_18, method = "ward.D2")
plot(hClustering_18, hang = -1, no.margin = T)

hClustering_19 <- hclust(distmat_19, method = "ward.D2")
plot(hClustering_19, hang = -1, no.margin = T)

colors <- c(
  "#762a83",
  "#1b7837",
  "#9970ab",
  "#5aae61",
  "#c2a5cf",
  "#a6dba0",
  "#e7d4e8",
  "#d9f0d3"
)
clus17 <- cutree(hClustering_17, 4)
plot(as.phylo(hClustering_17),
     tip.color = colors[clus17],
     label.offset = 0.01, cex = 0.7, no.margin = T
)

clus18 <- cutree(hClustering_18, 5)
plot(as.phylo(hClustering_18),
     tip.color = colors[clus18],
     label.offset = 0.01, cex = 0.7, no.margin = T
)

clus19 <- cutree(hClustering_19, 5)
plot(as.phylo(hClustering_19),
     tip.color = colors[clus19],
     label.offset = 0.01, cex = 0.7, no.margin = T
)

## Converting hclust to dendograms

hClustDen_17 <- as.dendrogram(hClustering_17)

ggdendrogram(hClustering_17, rotate = T)

hClustDen_18 <- as.dendrogram(hClustering_18)

ggdendrogram(hClustering_18, rotate = T)

hClustDen_19 <- as.dendrogram(hClustering_19)

ggdendrogram(hClustering_19, rotate = T)

hclustDenData_17 <- dendro_data(hClustDen_17)
hclustDenData_18 <- dendro_data(hClustDen_18)
hclustDenData_19 <- dendro_data(hClustDen_19)

p <- ggplot(hclustDenData_17$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(
    data = hclustDenData_17$labels, aes(x, y, label = label),
    hjust = 1, angle = 90, size = 3
  ) +
  ylim(-0.5, 1) 
p

p <- ggplot(hclustDenData_18$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(
    data = hclustDenData_18$labels, aes(x, y, label = label),
    hjust = 1, angle = 90, size = 2
  ) +
  ylim(-0.5, 1) 
p

p <- ggplot(hclustDenData_19$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(
    data = hclustDenData_19$labels, aes(x, y, label = label),
    hjust = 1, angle = 90, size = 2
  ) +
  ylim(-0.5, 1) 
p

## Neighbour-joining
nj_17 <- nj(distmat_17)
plot.phylo(nj_17,
           type = "phylogram", use.edge.length = F,
           node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
           edge.color = "black", edge.width = 1, edge.lty = 1, font = 1,
           cex = 0.75, adj = NULL, srt = 0, no.margin = T,
           root.edge = FALSE, label.offset = 0.25, underscore = FALSE,
           x.lim = NULL, y.lim = NULL, direction = "rightwards",
           lab4ut = NULL, tip.color = "black", plot = TRUE,
           rotate.tree = 0, open.angle = 0, node.depth = 1,
           align.tip.label = T
)


nj_18 <- nj(distmat_18)
plot.phylo(nj_18,
           type = "phylogram", use.edge.length = F,
           node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
           edge.color = "black", edge.width = 1, edge.lty = 1, font = 1,
           cex = 0.5, adj = NULL, srt = 0, no.margin = T,
           root.edge = FALSE, label.offset = 0.25, underscore = FALSE,
           x.lim = NULL, y.lim = NULL, direction = "rightwards",
           lab4ut = NULL, tip.color = "black", plot = TRUE,
           rotate.tree = 0, open.angle = 0, node.depth = 1,
           align.tip.label = T
)

nj_19 <- nj(distmat_19)
plot.phylo(nj_19,
           type = "phylogram", use.edge.length = F,
           node.pos = NULL, show.tip.label = TRUE, show.node.label = FALSE,
           edge.color = "black", edge.width = 1, edge.lty = 1, font = 1,
           cex = 0.5, adj = NULL, srt = 0, no.margin = T,
           root.edge = FALSE, label.offset = 0.25, underscore = FALSE,
           x.lim = NULL, y.lim = NULL, direction = "rightwards",
           lab4ut = NULL, tip.color = "black", plot = TRUE,
           rotate.tree = 0, open.angle = 0, node.depth = 1,
           align.tip.label = T
)

gplots::heatmap.2(as.matrix(distmat_17),
                  margins = c(8, 8), trace = "none",
                  dendrogram = "both",
                  density.info = "density",
                  col = "viridis",
                  main = "Hierarchichal Clustering of additive genetic effects 2016/2017 season"
)

gplots::heatmap.2(as.matrix(distmat_18),
                  margins = c(8, 8), trace = "none",
                  dendrogram = "both",
                  density.info = "density",
                  col = "viridis",
                  main = "Hierarchichal Clustering of additive genetic effects 2017/2018 season"
)

gplots::heatmap.2(as.matrix(distmat_19),
                  margins = c(8, 8), trace = "none",
                  dendrogram = "both",
                  density.info = "density",
                  col = "viridis",
                  main = "Hierarchichal Clustering of additive genetic effects 2018/2019 season"
)

## hierarchical clustering with significance
# Can take a long time

# pvClust_17_ward <- pvclust(as.matrix(dat17[, 2:ncol(dat17)]),
#   method.hclust = "ward.D2",
#   method.dist = "abscor",
#   use.cor = "pairwise.complete.obs",
#   nboot = 100000
# )
# beep(4)
#
# par(mar = c(0.5, 2, 2, 0.25))
#
# plot(pvClust_17_ward)
# pvrect(pvClust_17_ward, alpha = 0.95)
# par(mar = c(4, 4, 2, 0.25))
# x <- seplot(pvClust_17_ward, identify = TRUE)
# print(pvClust_17_ward, which = x)
#
# pvClust_17_ave <- pvclust(as.matrix(dat17[, 2:ncol(dat17)]),
#   method.hclust = "average",
#   method.dist = "abscor",
#   use.cor = "pairwise.complete.obs",
#   nboot = 100000
# )
# beep(5)
#
# par(mar = c(0.5, 2, 2, 0.25))
#
# plot(pvClust_17_ave)
# pvrect(pvClust_17_ave, alpha = 0.95)
# par(mar = c(4, 4, 2, 0.25))
# x <- seplot(pvClust_17_ave, identify = TRUE)
# print(pvClust_17_ave, which = x)
#
# pvClust_17_com <- pvclust(as.matrix(dat17[, 2:ncol(dat17)]),
#   method.hclust = "complete",
#   method.dist = "abscor",
#   use.cor = "pairwise.complete.obs",
#   nboot = 100000
# )
# beep(6)
#
# par(mar = c(0.5, 2, 2, 0.25))
#
# plot(pvClust_17_com)
# pvrect(pvClust_17_com, alpha = 0.95)
# par(mar = c(4, 4, 2, 0.25))
# x <- seplot(pvClust_17_com, identify = TRUE)
# print(pvClust_17_com, which = x)
#
# pvClust_18_ward <- pvclust(as.matrix(dat18[, 2:ncol(dat18)]),
#   method.hclust = "ward.D2",
#   method.dist = "abscor",
#   use.cor = "pairwise.complete.obs",
#   nboot = 100000
# )
# beep(7)
#
# par(mar = c(0.5, 2, 2, 0.25))
#
# plot(pvClust_18_ward)
# pvrect(pvClust_18_ward, alpha = 0.95)
# par(mar = c(4, 4, 2, 0.25))
# x <- seplot(pvClust_18_ward, identify = TRUE)
# print(pvClust_18_ward, which = x)
#
# pvClust_18_ave <- pvclust(as.matrix(dat18[, 2:ncol(dat18)]),
#   method.hclust = "average",
#   method.dist = "abscor",
#   use.cor = "pairwise.complete.obs",
#   nboot = 100000
# )
# beep(8)
#
# par(mar = c(0.5, 2, 2, 0.25))
#
# plot(pvClust_18_ave)
# pvrect(pvClust_18_ave, alpha = 0.95)
# par(mar = c(4, 4, 2, 0.25))
# x <- seplot(pvClust_18_ave, identify = TRUE)
# print(pvClust_18_ave, which = x)
#
# pvClust_18_com <- pvclust(as.matrix(dat18[, 2:ncol(dat18)]),
#   method.hclust = "complete",
#   method.dist = "abscor",
#   use.cor = "pairwise.complete.obs",
#   nboot = 100000
# )
# beep(8)
#
# par(mar = c(0.5, 2, 2, 0.25))
#
# plot(pvClust_18_com)
# pvrect(pvClust_18_com, alpha = 0.95)
# par(mar = c(4, 4, 2, 0.25))
# x <- seplot(pvClust_18_com, identify = TRUE)
# print(pvClust_18_com, which = x)
