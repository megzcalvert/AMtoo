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
    legend.key.size = unit(2, "lines"),
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
  arrange(Variety) %>%
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
  filter(year == "17") %>%
  arrange(Variety)

gryld18 <- pheno %>%
  filter(year == "18") %>%
  arrange(Variety)

gryld19 <- pheno %>%
  filter(year == "19") %>%
  arrange(Variety)

snpChip <- fread(
  file = "./Genotype_Database/SelectedImputedBeagleNumeric.txt",
  header = TRUE, check.names = F, sep = "\t"
)

snpMatrix <- t(snpChip[, c(-1, -2, -3)])

snpMatrix <- setDT(as.data.frame(snpMatrix), keep.rownames = T)

## 2017
gryld17 <- gryld17 %>%
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>%
  arrange(Variety)

snpMatrix17 <- snpMatrix %>%
  semi_join(gryld17, by = c("rn" = "Variety")) %>%
  arrange(rn) %>%
  column_to_rownames(var = "rn")

kinMat <- A.mat(as.matrix(snpMatrix17))

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
  arrange(Variety) %>%
  tidylog::select(-Variety)

gryld17 <- gryld17 %>%
  left_join(phenoVI_17, by = c("entity_id" = "Plot_ID")) %>%
  arrange(Variety)

#### Heritability Function ####
h2_function <- function(dat, snpMat, ...) {
  effectvars <- c("entity_id", "Variety", "rep", "year", "rep1", "rep2")
  traits <- dat %>%
    tidylog::select(-any_of(effectvars))
  traits <- colnames(traits)
  traits

  a <- c("VarAdditive", "VarEnvironment", "h2", "logLik")

  h2_all <- tibble::tibble(a)

  for (i in traits) {
    print(paste("Working on trait", i))
    h2_dat <- dat %>%
      tidylog::select("entity_id", "Variety", "rep1", "rep2", paste(i))
    names(h2_dat)[names(h2_dat) == paste(i)] <- "trait"

    snpMatrixUse <- snpMat %>%
      semi_join(h2_dat, by = c("rn" = "Variety")) %>%
      arrange(rn) %>%
      column_to_rownames(var = "rn")

    kinMat <- A.mat(as.matrix(snpMatrixUse))

    h2 <- marker_h2(
      data.vector = h2_dat$trait,
      geno.vector = h2_dat$Variety,
      covariates = h2_dat[, 3:4],
      K = kinMat,
      max.iter = 150
    )

    h2_f <- c(h2$va, h2$ve, h2$h2, h2$loglik)
    h2_all$newTrait <- h2_f
    colnames(h2_all)[colnames(h2_all) == "newTrait"] <- paste(i)
  }
  return(h2_all)
}

h2_17 <- h2_function(dat = gryld17, snpMat = snpMatrix)

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
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>%
  arrange(Variety)

snpMatrix18 <- snpMatrix %>%
  semi_join(gryld18, by = c("rn" = "Variety")) %>%
  arrange(rn)

phenoVI_18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
phenoVI_18 <- phenoVI_18 %>%
  tidylog::select(-GRYLD) %>%
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>%
  arrange(Variety) %>%
  tidylog::select(-Variety)

gryld18 <- gryld18 %>%
  left_join(phenoVI_18, by = c("entity_id")) %>%
  arrange(Variety) %>%
  tidylog::select(-year.x, -year.y)
colnames(gryld18)

h2_18 <- h2_function(dat = gryld18, snpMat = snpMatrix)

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
  semi_join(snpMatrix, by = c("Variety" = "rn")) %>%
  arrange(Variety)

snpMatrix19 <- snpMatrix %>%
  semi_join(gryld19, by = c("rn" = "Variety")) %>%
  arrange(rn)
colnames(gryld19)

phenoVI_19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoVI_19 <- phenoVI_19 %>%
  tidylog::select(-GRYLD, -year) %>%
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>%
  arrange(Variety) %>%
  tidylog::select(-Variety)

gryld19 <- gryld19 %>%
  left_join(phenoVI_19, by = c("entity_id")) %>%
  arrange(Variety)

colnames(gryld19)

h2_19 <- h2_function(dat = gryld19, snpMat = snpMatrix)

h2_19 <- setDT(as.data.frame(t(h2_19)), keep.rownames = T)
names(h2_19) <- c("trait_id", "VarAdditive", "VarEnvironment", "h2", "logLik")
h2_19 <- h2_19[-1, ]
write.table(h2_19,
  file = "./Phenotype_Database/NarrowSenseHeritability_19.txt",
  row.names = F, col.names = T, quote = F, sep = "\t"
)
h2_19 <- fread("./Phenotype_Database/NarrowSenseHeritability_19.txt")

##### Variance and covariance for correlation ####

pheno17 <- gryld17 %>%
  arrange(Variety) %>%
  tidylog::select(-entity_id, -Variety, -rep, -year, -rep1, -rep2)
coVariance_17 <- cov(as.matrix(pheno17))

pheno18 <- gryld18 %>%
  arrange(Variety) %>%
  tidylog::select(-entity_id, -Variety, -rep, -rep1, -rep2)
coVariance_18 <- cov(as.matrix(pheno18))

phenoVI_19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")

pheno19 <- phenoVI_19 %>%
  tidylog::select(-year) %>%
  unite("trait_id", ID:Date, sep = "_") %>%
  pivot_wider(names_from = trait_id, values_from = value) %>%
  arrange(Variety) %>%
  tidylog::select(-entity_id, -Variety) %>%
  glimpse()

coVariance_19 <- cov(as.matrix(pheno19))

#### function for expected gain from selection of correlated traits ####

expectedGainSelection <- function(heritability, phenotypic,
                                  traitInterest, intensity) {
  traits <- phenotypic %>%
    tidylog::select(-any_of(traitInterest))
  traits <- colnames(traits)

  heritability <- heritability %>%
    tidylog::select(trait_id, h2)

  correlatedTrait <- traits

  h2_interest <- heritability %>%
    filter(trait_id == traitInterest) %>%
    tidylog::select(h2)

  h2_interest <- h2_interest %>%
    pull(h2)

  gainCorrelated <- tibble::tibble(correlatedTrait, h2_interest)
  gainCorrelated <- add_column(gainCorrelated,
    h2_correlated = as.numeric(NA),
    geneticCorrelation = as.numeric(NA),
    expectedGain = as.numeric(NA)
  )
  for (i in traits) {
    print(paste("Working on trait", i))
    her <- heritability %>%
      filter(trait_id == i) %>%
      tidylog::select(h2)
    her <- her %>%
      pull(h2)

    gainCorrelated[
      match(i, gainCorrelated$correlatedTrait), "h2_correlated"
    ] <- her

    trait_interest <- phenotypic %>%
      pull(all_of(traitInterest))
    trait_correlated <- phenotypic %>%
      pull(all_of(i))

    geneticCorrelation <- cov(trait_interest, trait_correlated) /
      (sqrt(var(trait_interest) * var(trait_correlated)))

    print(paste("genetic correlation ", geneticCorrelation))

    gainCorrelated[
      match(i, gainCorrelated$correlatedTrait), "geneticCorrelation"
    ] <- geneticCorrelation

    expectedGain <- intensity * sqrt(h2_interest) * sqrt(her) *
      geneticCorrelation * (sd(trait_correlated) / sqrt(length(phenotypic)))
    gainCorrelated[
      match(i, gainCorrelated$correlatedTrait), "expectedGain"
    ] <- expectedGain
  }
  return(gainCorrelated)
}

gain17_2 <- expectedGainSelection(
  heritability = h2_17, phenotypic = pheno17,
  traitInterest = "GRYLD", intensity = 0.2
)

gain17_4 <- expectedGainSelection(
  heritability = h2_17, phenotypic = pheno17,
  traitInterest = "GRYLD", intensity = 0.4
)

gain17_6 <- expectedGainSelection(
  heritability = h2_17, phenotypic = pheno17,
  traitInterest = "GRYLD", intensity = 0.6
)

gain18_2 <- expectedGainSelection(
  heritability = h2_18, phenotypic = pheno18,
  traitInterest = "GRYLD", intensity = 0.2
)

gain18_4 <- expectedGainSelection(
  heritability = h2_18, phenotypic = pheno18,
  traitInterest = "GRYLD", intensity = 0.4
)

gain18_6 <- expectedGainSelection(
  heritability = h2_18, phenotypic = pheno18,
  traitInterest = "GRYLD", intensity = 0.6
)

gain19_2 <- expectedGainSelection(
  heritability = h2_19, phenotypic = pheno19,
  traitInterest = "GRYLD", intensity = 0.2
)

gain19_4 <- expectedGainSelection(
  heritability = h2_19, phenotypic = pheno19,
  traitInterest = "GRYLD", intensity = 0.4
)

gain19_6 <- expectedGainSelection(
  heritability = h2_19, phenotypic = pheno19,
  traitInterest = "GRYLD", intensity = 0.6
)

write.table(gain17_2, "./Phenotype_Database/ExpectedGainSelection_2_17.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain17_4, "./Phenotype_Database/ExpectedGainSelection_4_17.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain17_6, "./Phenotype_Database/ExpectedGainSelection_6_17.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain18_2, "./Phenotype_Database/ExpectedGainSelection_2_18.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain18_4, "./Phenotype_Database/ExpectedGainSelection_4_18.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain18_6, "./Phenotype_Database/ExpectedGainSelection_6_18.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain19_2, "./Phenotype_Database/ExpectedGainSelection_2_19.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain19_4, "./Phenotype_Database/ExpectedGainSelection_4_19.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(gain19_6, "./Phenotype_Database/ExpectedGainSelection_6_19.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

gain17_2<- fread("./Phenotype_Database/ExpectedGainSelection0.2_17.txt")
gain17_4<- fread("./Phenotype_Database/ExpectedGainSelection0.4_17.txt")
gain17_6<- fread("./Phenotype_Database/ExpectedGainSelection0.6_17.txt")
gain18_2<- fread("./Phenotype_Database/ExpectedGainSelection0.2_18.txt")
gain18_4<- fread("./Phenotype_Database/ExpectedGainSelection0.4_18.txt")
gain18_6<- fread("./Phenotype_Database/ExpectedGainSelection0.6_18.txt")
gain19_2<- fread("./Phenotype_Database/ExpectedGainSelection0.2_19.txt")
gain19_4<- fread("./Phenotype_Database/ExpectedGainSelection0.4_19.txt")
gain19_6<- fread("./Phenotype_Database/ExpectedGainSelection0.6_19.txt")

gain17 <- gain17_2 %>%
  mutate(intensity = as.factor(0.2)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

gain17 <- gain17_4 %>%
  mutate(intensity = as.factor(0.4)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d")) %>%
  bind_rows(gain17)

gain17 <- gain17_6 %>%
  mutate(intensity = as.factor(0.6)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d")) %>%
  bind_rows(gain17)

gain18 <- gain18_2 %>%
  mutate(intensity = as.factor(0.2)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  filter(trait_id != "height") %>% 
  filter(phenotype_date > "2018-02-01") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

gain18 <- gain18_4 %>%
  mutate(intensity = as.factor(0.4)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  filter(trait_id != "height") %>% 
  filter(phenotype_date > "2018-02-01") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d")) %>%
  bind_rows(gain18)

gain18 <- gain18_6 %>%
  mutate(intensity = as.factor(0.6)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  filter(trait_id != "height") %>% 
  filter(phenotype_date > "2018-02-01") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d")) %>%
  bind_rows(gain18)

gain19 <- gain19_2 %>%
  mutate(intensity = as.factor(0.2)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  filter(trait_id != "height") %>% 
  filter(phenotype_date > "2019-02-01") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d"))

gain19 <- gain19_4 %>%
  mutate(intensity = as.factor(0.4)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d")) %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  filter(trait_id != "height") %>% 
  filter(phenotype_date > "2019-02-01") %>% 
  bind_rows(gain19)

gain19 <- gain19_6 %>%
  mutate(intensity = as.factor(0.6)) %>%
  separate(correlatedTrait, c("trait_id", "phenotype_date"), sep = "_") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  filter(trait_id != "height") %>% 
  filter(phenotype_date > "2019-02-01") %>% 
  mutate(phenotype_date = as.Date(phenotype_date, format = "%Y-%m-%d")) %>%
  bind_rows(gain19)

plot17 <- gain17 %>%
  ggplot(aes(x = phenotype_date, y = expectedGain, colour = intensity)) +
  geom_point() +
  facet_wrap(~trait_id, scales = "free_y") +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "Expected gain in GRYLD from correlated phenotype 2016-2017")
plot17

plot18 <- gain18 %>%
  ggplot(aes(x = phenotype_date, y = expectedGain, colour = intensity)) +
  geom_point() +
  facet_wrap(~trait_id, scales = "free_y") +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "Expected gain in GRYLD from correlated phenotype 2017-2018")
plot18

plot19 <- gain19 %>%
  ggplot(aes(x = phenotype_date, y = expectedGain, colour = intensity)) +
  geom_point() +
  facet_wrap(~trait_id, scales = "free_y") +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "Expected gain in GRYLD from correlated phenotype 2018-2019")
plot19
