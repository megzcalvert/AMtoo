rm(list = objects())
ls()

library(data.table)
library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(heritability)
library(rrBLUP)

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
    legend.key.size = unit(4, "lines"),
    # legend.background = element_rect(fill = NULL, colour = NULL),
    # legend.box = NULL,
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
      size = rel(3),
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
      size = rel(2),
      hjust = 0
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 1
    ),
    strip.text = element_text(
      colour = "black",
      size = rel(1)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")
set.seed(1642)
options(scipen = 999)
sessionInfo()

###############################################################################
#### Load data ####

pheno <- fread("./Phenotype_Database/Pheno_Long171819.txt")

unique(pheno$trait_id)

pheno <- pheno %>%
  filter(trait_id == "GRYLD") %>%
  tidylog::select(
    -phenotype_date, -phenotype_person,
    -block, -range, -column, -trait_id
  ) %>%
  rename(GRYLD = phenotype_value) %>%
  mutate(
    rep1 = if_else(rep == 1, 1, 0),
    rep2 = if_else(rep == 2, 1, 0)
  )

gryld17 <- pheno %>%
  filter(year == "17")

gryld18 <- pheno %>%
  filter(year == "18")

gryld19 <- pheno %>%
  filter(year == "19")

snpChip <- fread(
  file = "./Genotype_Database/SelectedImputedBeagleNumeric.txt",
  header = TRUE, check.names = F, sep = "\t"
)

snpMatrix <- t(snpChip[, c(-1, -2, -3)])
kinMat <- A.mat(snpMatrix)

snpMatrix <- setDT(as.data.frame(snpMatrix), keep.rownames = T)

## 2017
gryld17 <- gryld17 %>%
  semi_join(snpMatrix, by = c("Variety" = "rn"))

snpMatrix17 <- snpMatrix %>%
  semi_join(gryld17, by = c("rn" = "Variety"))

gryld_h2 <- marker_h2(
  data.vector = gryld17$GRYLD,
  geno.vector = gryld17$Variety,
  covariates = gryld17[, 6:7],
  K = kinMat
)

gryld_h2

phenoVI_17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
phenoVI_17 <- phenoVI_17 %>%
  tidylog::select(-GRYLD) %>%
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>%
  tidylog::select(-Variety)

gryld17 <- gryld17 %>%
  left_join(phenoVI_17, by = c("entity_id" = "Plot_ID"))

#### Heritability Function ####
h2_function <- function(dat, kinshipMat, ...) {
  effectvars <- names(dat) %in% c(
    "entity_id", "Variety", "rep", "year", "rep1", "rep2"
  )
  traits <- colnames(dat[, !effectvars])
  traits

  a <- c("VarAdditive", "VarEnvironment", "h2", "logLik")

  h2_all <- tibble::tibble(a)

  for (i in traits) {
    print(paste("Working on trait", i))
    h2_dat <- dat %>%
      tidylog::select("entity_id", "Variety", "rep1", "rep2", paste(i))
    h2 <- marker_h2(
      data.vector = h2_dat[, 5],
      geno.vector = h2_dat[, 2],
      covariates = h2_dat[, 3:4],
      K = kinshipMat,
      max.iter = 150
    )

    h2_f <- c(h2$va, h2$ve, h2$h2, h2$loglik)
    h2_all$newTrait <- h2_f
    colnames(h2_all)[colnames(h2_all) == "newTrait"] <- paste(i)
  }
  return(h2_all)
}

h2_17 <- h2_function(dat = gryld17, kinshipMat = kinMat)

h2_17 <- setDT(as.data.frame(t(h2_17)), keep.rownames = T)
names(h2_17) <- c("trait_id", "VarAdditive", "VarEnvironment", "h2", "logLik")
h2_17 <- h2_17[-1, ]
write.table(h2_17,
  file = "./Phenotype_Database/NarrowSenseHeritability_17.txt",
  row.names = F, col.names = T, quote = F, sep = "\t"
)
h2_17 <- fread("./Phenotype_Database/NarrowSenseHeritability_17.txt")

## 2018
gryld18 <- gryld18 %>%
  semi_join(snpMatrix, by = c("Variety" = "rn"))

snpMatrix18 <- snpMatrix %>%
  semi_join(gryld18, by = c("rn" = "Variety"))

phenoVI_18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
phenoVI_18 <- phenoVI_18 %>%
  tidylog::select(-GRYLD) %>%
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>%
  tidylog::select(-Variety)

gryld18 <- gryld18 %>%
  left_join(phenoVI_18, by = c("entity_id")) %>%
  tidylog::select(-year.x, -year.y)
colnames(gryld18)

h2_18 <- h2_function(dat = gryld18, kinshipMat = kinMat)

h2_18 <- setDT(as.data.frame(t(h2_18)), keep.rownames = T)
names(h2_18) <- c("trait_id", "VarAdditive", "VarEnvironment", "h2", "logLik")
h2_18 <- h2_18[-1, ]
write.table(h2_18,
  file = "./Phenotype_Database/NarrowSenseHeritability_18.txt",
  row.names = F, col.names = T, quote = F, sep = "\t"
)
h2_18 <- fread("./Phenotype_Database/NarrowSenseHeritability_18.txt")

## 2019

gryld19 <- gryld19 %>%
  semi_join(snpMatrix, by = c("Variety" = "rn"))

snpMatrix19 <- snpMatrix %>%
  semi_join(gryld19, by = c("rn" = "Variety"))
colnames(gryld19)

phenoVI_19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoVI_19 <- phenoVI_19 %>%
  tidylog::select(-GRYLD, -year) %>%
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>%
  tidylog::select(-Variety) 

gryld19 <- gryld19 %>%
  left_join(phenoVI_19, by = c("entity_id"))
colnames(gryld19)

h2_19 <- h2_function(dat = gryld19, kinshipMat = kinMat)

h2_19 <- setDT(as.data.frame(t(h2_19)), keep.rownames = T)
names(h2_19) <- c("trait_id", "VarAdditive", "VarEnvironment", "h2", "logLik")
h2_19 <- h2_19[-1, ]
write.table(h2_19,
  file = "./Phenotype_Database/NarrowSenseHeritability_19.txt",
  row.names = F, col.names = T, quote = F, sep = "\t"
)
h2_19 <- fread("./Phenotype_Database/NarrowSenseHeritability_19.txt")

##### Variance and covariance for correlation ####

coVariance_17 <- gryld17 %>%
  tidylog::select(-entity_id, -Variety, -rep, -year, -rep1, -rep2)
coVariance_17 <- cov(as.matrix(coVariance_17))

coVariance_18 <- gryld18 %>%
  tidylog::select(-entity_id, -Variety, -rep, -rep1, -rep2)
coVariance_18 <- cov(as.matrix(coVariance_18))

phenoVI_19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")

pheno19<- phenoVI_19 %>% 
  tidylog::select(-year) %>% 
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>% 
  tidylog::select(-entity_id,-Variety) %>% 
  glimpse()

coVariance_19<- cov(as.matrix(pheno19))

