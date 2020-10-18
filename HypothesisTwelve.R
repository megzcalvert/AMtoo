rm(list = objects())
ls()

library(janitor)
library(tidyverse)
library(tidylog)
library(readr)
library(data.table)
library(rrBLUP)
library(broom)

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
      size = rel(1.5)
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
  clean_names() # janitor package, gives a consistent formatting

# Get rid of the ambiguous and the het codes
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

# Change to numeric format needed
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

# Number of markers per chromosome
chrSum <- plyr::count(snpChip, vars = "chrom")

# transpose and remove positions
snpMatrix <- t(snpChip[, c(-1, -2, -3)])
snpMatrix %>% glimpse()
snpMatrix[1:10, 1:10]

snpMatrix <- setDT(as.data.frame(snpMatrix), keep.rownames = T)
snpMatrix[1:10, 1:10]

##### Phenotypes ####

pheno17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoLong <- fread("./Phenotype_Database/Pheno_Long171819.txt")

dat17 <- fread("./Phenotype_Database/Hyp10BLUEs_17.txt")
dat18 <- fread("./Phenotype_Database/Hyp10BLUEs_18.txt")
dat19 <- fread("./Phenotype_Database/Hyp10BLUEs_19.txt")

###############################################################
####                    rrBlup trial                      ####

snpMatrix[1:5, 1:5]

## Keep only lines which have geno and pheno
# 2017
dat17 <- dat17 %>%
  semi_join(snpMatrix, by = "rn") %>%
  arrange(var = "rn")

snpMatrix17 <- snpMatrix %>%
  semi_join(dat17, by = "rn") %>%
  arrange(var = "rn") %>%
  column_to_rownames(var = "rn")

snpMatrix17 <- as.matrix(snpMatrix17)

snpMatrix17 <- snpMatrix17[match(
  dat17$rn,
  rownames(snpMatrix17)
), ] # align markers

# 2018
dat18 <- dat18 %>%
  semi_join(snpMatrix, by = "rn") %>%
  arrange(var = "rn")

snpMatrix18 <- snpMatrix %>%
  semi_join(dat18, by = "rn") %>%
  arrange(var = "rn") %>%
  column_to_rownames(var = "rn")

snpMatrix18 <- as.matrix(snpMatrix18)
snpMatrix18 <- snpMatrix18[match(
  dat18$rn,
  rownames(snpMatrix18)
), ] # align markers

# 2019
dat19 <- dat19 %>%
  semi_join(snpMatrix, by = "rn") %>%
  arrange(var = "rn")

snpMatrix19 <- snpMatrix %>%
  semi_join(dat19, by = "rn") %>%
  arrange(var = "rn") %>%
  column_to_rownames(var = "rn")

snpMatrix19 <- as.matrix(snpMatrix19)
snpMatrix19 <- snpMatrix19[match(
  dat19$rn,
  rownames(snpMatrix19)
), ] # align markers

##### Defining the training and test populations ####
# define the training and test populations
# training-80% validation-20%
dat17 <- dat17 %>%
  column_to_rownames(var = "rn")
str(dat17)

relationshipMat <- A.mat(snpMatrix17)

all(rownames(dat17) == rownames(relationshipMat))

y <- dat17[, "GRYLD"]
str(y)

y.NA <- dat17[, "GRYLD"]
mask <- sample(length(y), size = (0.2 * nrow(dat17)))
y.NA[mask] <- NA

# rrBLUP
g.blup <- mixed.solve(y.NA, K = relationshipMat)
#
pdf(file = "./Figures/GP_gryld17_phenoVSpred_rrBlup.pdf")
plot(g.blup$u, y,
  xlab = "Phenotype",
  ylab = "Pred. Gen. Value", cex = .8, bty = "L"
)
points(x = y[mask], y = g.blup$u[mask], col = 2, cex = .8, pch = 19)
legend("topleft",
  legend = c("training", "testing"), bty = "n",
  pch = c(1, 19), col = c("black", "red")
)
#
ymasked <- y[mask]
blupped <- g.blup$u[mask]
#
yunmasked <- y[-mask]
unblupped <- g.blup$u[-mask]
#
YLD_accuracy_gryldA <- cor(ymasked, blupped)
cor(x = yunmasked, y = unblupped)
YLD_accuracy_gryldA

## VI
y <- dat17[, "NDVI_20170512"]
str(y)

y.NA <- dat17[, "NDVI_20170512"]
mask <- sample(length(y), size = (0.2 * nrow(dat17)))
y.NA[mask] <- NA

# rrBLUP
g.blup <- mixed.solve(y.NA, K = relationshipMat)
#
pdf(file = "./Figures/GP_ndvi20170512_phenoVSpred_rrBlup.pdf")
plot(g.blup$u, y,
  xlab = "Phenotype",
  ylab = "Pred. Gen. Value", cex = .8, bty = "L"
)
points(x = y[mask], y = g.blup$u[mask], col = 2, cex = .8, pch = 19)
legend("topleft",
  legend = c("training", "testing"), bty = "n",
  pch = c(1, 19), col = c("black", "red")
)
#
ymasked <- y[mask]
blupped <- g.blup$u[mask]
#
yunmasked <- y[-mask]
unblupped <- g.blup$u[-mask]
#
accuracy_ndvi0512 <- cor(ymasked, blupped)
cor(x = yunmasked, y = unblupped)
accuracy_ndvi0512

y <- dat17[, "RedEdge_20170512"]
str(y)

y.NA <- dat17[, "NDVI_20170512"]
mask <- sample(length(y), size = (0.2 * nrow(dat17)))
y.NA[mask] <- NA

# rrBLUP
g.blup <- mixed.solve(y.NA, K = relationshipMat)
#
pdf(file = "./Figures/GP_re20170512_phenoVSpred_rrBlup.pdf")
plot(g.blup$u, y,
  xlab = "Phenotype",
  ylab = "Pred. Gen. Value", cex = .8, bty = "L"
)
points(x = y[mask], y = g.blup$u[mask], col = 2, cex = .8, pch = 19)
legend("topleft",
  legend = c("training", "testing"), bty = "n",
  pch = c(1, 19), col = c("black", "red")
)
#
ymasked <- y[mask]
blupped <- g.blup$u[mask]
#
yunmasked <- y[-mask]
unblupped <- g.blup$u[-mask]
#
accuracy_re0512 <- cor(ymasked, blupped)
cor(x = yunmasked, y = unblupped)
accuracy_re0512

##### Cross-Validation 2017 ####

genomicSelcectionCV <- function(dat, columnsToRemove, cycles, snpMatrix, ...) {
  traits <- dat %>%
    dplyr::select(!any_of(columnsToRemove)) %>%
    colnames()

  accuracy <- matrix(nrow = cycles, ncol = length(traits))
  colnames(accuracy) <- traits

  for (r in 1:cycles) {
    print(paste("Rep cycle: ", r))

    mask <- sample(nrow(dat), size = (0.2 * nrow(dat)))

    relationshipMat <- A.mat(snpMatrix)

    for (i in traits) {
      print(paste(i))

      y <- dat %>%
        pull(i)

      y.NA <- dat %>%
        pull(i)
      y.NA[mask] <- NA

      g.blup <- mixed.solve(y.NA, K = relationshipMat)
      ymasked <- y[mask]
      blupped <- g.blup$u[mask]

      yunmasked <- y[-mask]
      unblupped <- g.blup$u[-mask]

      accuracy[r, paste(i)] <- (cor(ymasked,
        y = blupped,
        use = "complete"
      ))
    }
  }
  return(accuracy)
}

accuracy17 <- genomicSelcectionCV(
  dat = dat17, columnsToRemove = c("rn"),
  cycles = 100,
  snpMatrix = snpMatrix17
)

write.csv(accuracy17,
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy17.txt",
  quote = F, row.names = F
)
colMeans(accuracy17)

### covariate

relationshipMat <- A.mat(snpMatrix17)

all(rownames(dat17) == rownames(relationshipMat))

y <- dat17[, "GRYLD"]
str(y)

y.NA <- dat17[, "GRYLD"]
mask <- sample(length(y), size = (0.2 * nrow(dat17)))
y.NA[mask] <- NA

# rrBLUP
g.blup <- mixed.solve(y.NA, K = relationshipMat, X = dat17[, "NDVI_20170512"])
#
pdf(file = "./Figures/GP_gryld17covarNdvi20170512_phenoVSpred_rrBlup.pdf")
plot(g.blup$u, y,
  xlab = "Phenotype",
  ylab = "Pred. Gen. Value", cex = .8, bty = "L"
)
points(x = y[mask], y = g.blup$u[mask], col = 2, cex = .8, pch = 19)
legend("topleft",
  legend = c("training", "testing"), bty = "n",
  pch = c(1, 19), col = c("black", "red")
)
#
ymasked <- y[mask]
blupped <- g.blup$u[mask]
#
yunmasked <- y[-mask]
unblupped <- g.blup$u[-mask]
#
YLD_accuracy_gryldA <- cor(ymasked, blupped)
cor(x = yunmasked, y = unblupped)
YLD_accuracy_gryldA

covariateCV <- function(dat, snpMat, mainTrait, cycles, ...) {
  
  traits <- dat %>%
    select(-all_of(mainTrait))
  traits <- colnames(traits)

  acc <- matrix(nrow = cycles, ncol = length(traits) + 1)
  colnames(acc) <- c("GRYLD_alone", traits)

  for (i in 1:cycles) {
    print(paste("Cycles: ", i))
    
    y <- dat %>% 
      pull(mainTrait)
    y.NA <- dat %>% 
      pull(mainTrait)
    mask <- sample(length(y), size = (0.2 * nrow(dat)))
    y.NA[mask] <- NA
    
    relationshipMat<- A.mat(snpMat)

    trt_answer <- mixed.solve(y.NA,
      K = relationshipMat,
      SE = TRUE
    )

    pred <- trt_answer$u
    trt_accuracy <- cor(y[mask], pred[mask])
    print(trt_accuracy)
    acc[i, 1] <- trt_accuracy

    for (t in traits) {
      print(paste(t))

      covar <- dat[, t]

      trt_covar1_answer <- mixed.solve(y.NA,
        K = relationshipMat,
        X = covar,
        SE = TRUE
      )

      pred_covar1 <- trt_covar1_answer$u
      trt_covar1_accuracy <- cor(y[mask], pred_covar1[mask])
      print(trt_covar1_accuracy)
      acc[i, t] <- trt_covar1_accuracy
    }
  }
  return(acc)
}

covariate17 <- covariateCV(
  dat = dat17, snpMat = snpMatrix17,
  mainTrait = "GRYLD", cycles = 100
)

write.csv(covariate17,
  "./R/rrBlup/HypothesisTwelve/GS_80_100_accuracy17_covariate.txt",
  quote = F, row.names = F
)
colMeans(covariate17)

# 2018
accuracy18 <- genomicSelcectionCV(
  dat = dat18, columnsToRemove = c("rn"),
  cycles = 100,
  snpMatrix = snpMatrix18
)

write.csv(accuracy18,
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy18.txt",
  quote = F, row.names = F
)
colMeans(accuracy18)

dat18<- dat18 %>% 
  column_to_rownames(var = "rn")

covariate18 <- covariateCV(
  dat = dat18, snpMat = snpMatrix18,
  mainTrait = "GRYLD", cycles = 100
)

write.csv(covariate18,
  "./R/rrBlup/HypothesisTwelve/GS_80_100_accuracy18_covariate.txt",
  quote = F, row.names = F
)
colMeans(covariate18)

# 2019
accuracy19 <- genomicSelcectionCV(
  dat = dat19, columnsToRemove = c("rn"),
  cycles = 100,
  snpMatrix = snpMatrix19
)

write.csv(accuracy19,
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy19.txt",
  quote = F, row.names = F
)
colMeans(accuracy19)

dat19<- dat19 %>% 
  column_to_rownames(var = "rn")

covariate19 <- covariateCV(
  dat = dat19, snpMat = snpMatrix19,
  mainTrait = "GRYLD", cycles = 100
)

write.csv(covariate19,
  "./R/rrBlup/HypothesisTwelve/GS_80_100_accuracy19_covariate.txt",
  quote = F, row.names = F
)
colMeans(covariate19)

###############################################################################
#### Examining results ####
getwd()

gs17_100 <- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy17.txt",
  header = T
)
gs17_covar_100 <- fread(
  "./R/rrBlup/HypothesisTwelve/GS_80_100_accuracy17_covariate.txt",
  header = T
)

gs18_100 <- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy18.txt",
  header = T
)
gs18_covar_100 <- fread(
  "./R/rrBlup/HypothesisTwelve/GS_80_100_accuracy18_covariate.txt",
  header = T
)

gs19_100 <- fread(
  "./R/rrBlup/HypothesisTwelve/GenomicSelection_80_100_accuracy19.txt",
  header = T
)
gs19_covar_100 <- fread(
  "./R/rrBlup/HypothesisTwelve/GS_80_100_accuracy19_covariate.txt",
  header = T
)

gs17_100<- gs17_100 %>% 
  mutate(Cycle = row_number()) %>% 
  select(Cycle, GRYLD, everything()) %>%  
  pivot_longer(cols = RedEdge_20170609:GNDVI_20170331,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20170613"),
    phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

gs17_covar_100<- gs17_covar_100 %>% 
  mutate(Cycle = row_number()) %>% 
  select(Cycle, GRYLD_alone, everything()) %>%  
  pivot_longer(cols = RedEdge_20170602:GNDVI_20170331,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20170613"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

gs18_100<- gs18_100 %>% 
  mutate(Cycle = row_number()) %>% 
  select(Cycle, GRYLD, everything()) %>%  
  pivot_longer(cols = RE_20180613:GNDVI_20171120,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20180615"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

gs18_covar_100<- gs18_covar_100 %>% 
  mutate(Cycle = row_number()) %>% 
  select(Cycle, GRYLD_alone, everything()) %>%  
  pivot_longer(cols = RE_20180613:GNDVI_20171120,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20180615"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

colnames(gs19_100)

gs19_100<- gs19_100 %>% 
  mutate(Cycle = row_number()) %>% 
  select(Cycle, GRYLD, everything()) %>%  
  pivot_longer(cols = RE_20190624:GNDVI_20190103,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20190627"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

gs19_covar_100<- gs19_covar_100 %>% 
  mutate(Cycle = row_number()) %>% 
  select(Cycle, GRYLD_alone, everything()) %>%  
  pivot_longer(cols = RE_20190624:GNDVI_20190103,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20190627"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"))

## Figures
gs17_fig<- ggplot(data = gs17_100,
                  mapping = aes(x = phenotype_date,
                                y = Correlation,
                                group = phenotype_date)) +
  geom_rect(ymin = mean(gs17_100$GRYLD) - sd(gs17_100$GRYLD),
            ymax = mean(gs17_100$GRYLD) + sd(gs17_100$GRYLD),
            xmin = as.Date("2017-03-31"),
            xmax = as.Date("2017-06-09"),
            colour = "#bdbdbd",
            fill = "#bdbdbd") +
  geom_hline(yintercept = mean(gs17_100$GRYLD),
             colour = "blue",
             linetype = 2) +
  geom_jitter(alpha = 0.05) +
  geom_violin(fill = NA,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "2016/2017 Season",
       subtitle = "No covariate")
gs17_fig

gs17_covar_fig<- ggplot(data = gs17_covar_100,
                  mapping = aes(x = phenotype_date,
                                y = Correlation,
                                group = phenotype_date)) +
  geom_rect(ymin = mean(gs17_covar_100$GRYLD_alone) - sd(gs17_covar_100$GRYLD_alone),
            ymax = mean(gs17_covar_100$GRYLD_alone) + sd(gs17_covar_100$GRYLD_alone),
            xmin = as.Date("2017-03-31"),
            xmax = as.Date("2017-06-09"),
            colour = "#bdbdbd",
            fill = "#bdbdbd") +
  geom_hline(yintercept = mean(gs17_covar_100$GRYLD_alone),
             colour = "blue",
             linetype = 2) +
  geom_jitter(alpha = 0.05) +
  geom_violin(fill = NA,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "2016/2017 Season",
       subtitle = "GP for GRYLD, VI as covariate")
gs17_covar_fig

gs18_fig<- ggplot(data = gs18_100,
                  mapping = aes(x = phenotype_date,
                                y = Correlation,
                                group = phenotype_date)) +
  geom_rect(ymin = mean(gs18_100$GRYLD) - sd(gs18_100$GRYLD),
            ymax = mean(gs18_100$GRYLD) + sd(gs18_100$GRYLD),
            xmin = as.Date("2017-11-20"),
            xmax = as.Date("2018-06-13"),
            colour = "#bdbdbd",
            fill = "#bdbdbd") +
  geom_hline(yintercept = mean(gs18_100$GRYLD),
             colour = "blue",
             linetype = 2) +
  geom_jitter(alpha = 0.05) +
  geom_violin(fill = NA,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "2017/2018 Season",
       subtitle = "No covariate")
gs18_fig

gs18_covar_fig<- ggplot(data = gs18_covar_100,
                        mapping = aes(x = phenotype_date,
                                      y = Correlation,
                                      group = phenotype_date)) +
  geom_rect(ymin = mean(gs18_covar_100$GRYLD_alone) - sd(gs18_covar_100$GRYLD_alone),
            ymax = mean(gs18_covar_100$GRYLD_alone) + sd(gs18_covar_100$GRYLD_alone),
            xmin = as.Date("2017-11-20"),
            xmax = as.Date("2018-06-13"),
            colour = "#bdbdbd",
            fill = "#bdbdbd") +
  geom_hline(yintercept = mean(gs18_covar_100$GRYLD_alone),
             colour = "blue",
             linetype = 2) +
  geom_jitter(alpha = 0.05) +
  geom_violin(fill = NA,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "2017/2018 Season",
       subtitle = "GP for GRYLD, VI as covariate")
gs18_covar_fig

gs19_fig<- ggplot(data = gs19_100,
                  mapping = aes(x = phenotype_date,
                                y = Correlation,
                                group = phenotype_date)) +
  geom_rect(ymin = mean(gs19_100$GRYLD) - sd(gs19_100$GRYLD),
            ymax = mean(gs19_100$GRYLD) + sd(gs19_100$GRYLD),
            xmin = as.Date("2019-01-03"),
            xmax = as.Date("2019-06-24"),
            colour = "#bdbdbd",
            fill = "#bdbdbd") +
  geom_hline(yintercept = mean(gs19_100$GRYLD),
             colour = "blue",
             linetype = 2) +
  geom_jitter(alpha = 0.05) +
  geom_violin(fill = NA,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "2018/2019 Season",
       subtitle = "No covariate")
gs19_fig

gs19_covar_fig<- ggplot(data = gs19_covar_100,
                        mapping = aes(x = phenotype_date,
                                      y = Correlation,
                                      group = phenotype_date)) +
  geom_rect(ymin = mean(gs19_covar_100$GRYLD_alone) - sd(gs19_covar_100$GRYLD_alone),
            ymax = mean(gs19_covar_100$GRYLD_alone) + sd(gs19_covar_100$GRYLD_alone),
            xmin = as.Date("2019-01-03"),
            xmax = as.Date("2019-06-24"),
            colour = "#bdbdbd",
            fill = "#bdbdbd") +
  geom_hline(yintercept = mean(gs19_covar_100$GRYLD_alone),
             colour = "blue",
             linetype = 2) +
  geom_jitter(alpha = 0.05) +
  geom_violin(fill = NA,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~trait_id, scales = "free") +
  coord_cartesian(ylim = c(-1,1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(title = "2018/2019 Season",
       subtitle = "GP for GRYLD, VI as covariate")
gs19_covar_fig
