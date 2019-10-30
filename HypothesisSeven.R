rm(list = objects())
ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(asreml)
library(lme4)

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
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 1
    ),
    strip.text = element_text(
      colour = "black",
      size = rel(2)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")
set.seed(1642)

#### Load data ####
pheno17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
phenoLong <- fread("./Phenotype_Database/Pheno_Long171819.txt")

glimpse(pheno17)
glimpse(pheno18)
glimpse(pheno19)
glimpse(phenoLong)

phenoLong <- phenoLong %>%
  dplyr::rename(Plot_ID = entity_id) %>%
  tidylog::select(Plot_ID, Variety, block, rep, range, column)

pheno17$Date <- as.Date(pheno17$Date, format = "%Y-%m-%d")
pheno17$Date <- format(pheno17$Date, "%Y%m%d")
pheno18$Date <- as.Date(pheno18$Date, format = "%Y-%m-%d")
pheno18$Date <- format(pheno18$Date, "%Y%m%d")
pheno19$Date <- as.Date(pheno19$Date, format = "%Y-%m-%d")
pheno19$Date <- format(pheno19$Date, "%Y%m%d")

pheno17 <- pheno17 %>%
  unite("ID", c("ID", "Date")) %>%
  spread(key = ID, value = value) %>%
  tidylog::select(
    Plot_ID, Variety, GRYLD,
    GNDVI_20170331:RedEdge_20170609
  ) %>%
  tidylog::inner_join(phenoLong) %>%
  distinct()

pheno18 <- pheno18 %>%
  dplyr::rename(Plot_ID = entity_id) %>%
  unite("ID", c("ID", "Date")) %>%
  spread(key = ID, value = value) %>%
  tidylog::select(
    Plot_ID, Variety, GRYLD,
    GNDVI_20171120:RE_20180613
  ) %>%
  tidylog::inner_join(phenoLong) %>%
  distinct()

pheno19 <- pheno19 %>%
  dplyr::rename(Plot_ID = entity_id) %>%
  unite("ID", c("ID", "Date")) %>%
  spread(key = ID, value = value) %>%
  tidylog::select(
    Plot_ID, Variety, GRYLD,
    GNDVI_20190103:RE_20190624
  ) %>%
  tidylog::inner_join(phenoLong) %>%
  distinct()

###############################################################################
####                  ASREML to calculate heritability                  ####

pheno17$Plot_ID <- as.factor(pheno17$Plot_ID)
pheno17$Variety <- as.factor(pheno17$Variety)
pheno17$block <- as.factor(pheno17$block)
pheno17$rep <- as.factor(pheno17$rep)
pheno17$range <- as.factor(pheno17$range)
pheno17$column <- as.factor(pheno17$column)

pheno18$Plot_ID <- as.factor(pheno18$Plot_ID)
pheno18$Variety <- as.factor(pheno18$Variety)
pheno18$block <- as.factor(pheno18$block)
pheno18$rep <- as.factor(pheno18$rep)
pheno18$range <- as.factor(pheno18$range)
pheno18$column <- as.factor(pheno18$column)

pheno19$Plot_ID <- as.factor(pheno19$Plot_ID)
pheno19$Variety <- as.factor(pheno19$Variety)
pheno19$block <- as.factor(pheno19$block)
pheno19$rep <- as.factor(pheno19$rep)
pheno19$range <- as.factor(pheno19$range)
pheno19$column <- as.factor(pheno19$column)

# 2017 trial

# Without the auto-reggression correlation

asreml.license.status()

t17 <- asreml(
  fixed = GRYLD ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno17
)
plot(t17)
coef(t17)$random
fitted(t17)
summary(t17)
resid(t17)

h <- as.data.frame(summary(t17)$varcomp)
h

h2 <- as.data.frame(h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
h2

## 2017 all
effectvars <- names(pheno17) %in% c(
  "block", "rep", "Variety", "year", "column",
  "range", "Plot_ID"
)
traits <- colnames(pheno17[, !effectvars])
H2_2017 <- data.frame(traits)
H2_2017$Heritability <- NA
fieldInfo <- pheno17 %>%
  tidylog::select(Variety, rep, block, column, range)
ntraits <- 1:nrow(H2_2017)

for (i in ntraits) {
  print(paste("Working on trait", H2_2017[i, 1]))
  j <- H2_2017[i, 1]

  data <- cbind(fieldInfo, pheno17[, paste(j)])
  names(data) <- c("Variety", "rep", "block", "column", "range", "Trait")

  t17 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + rep + rep:block,
    data = data
  )
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock17_", H2_2017[i, 1], ".pdf"))
  plot(t17)
  dev.off()
  h <- as.data.frame(summary(t17)$varcomp)
  print(paste("Creating Data Frame", j))
  print(h)

  h2 <- (h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
  h2
  H2_2017[i, 2] <- h2
}

dev.off()

# lme4 comparison 2017
calcH2r <- function(dat, fill = NA, ...) {
  r <- length(table(dat$rep))

  effectvars <- names(dat) %in% c(
    "block", "rep", "Variety", "year", "column",
    "range", "Plot_ID"
  )

  t <- colnames(dat[, !effectvars])

  for (i in t) {
    print(paste("Working on trait", i))
    h <- as.data.frame(VarCorr(
      lmer(paste0(i, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)
    ))
    H2 <- h[1, 4] / (h[1, 4] + (h[4, 4] / r))
    print(H2)
  }
}

calcH2r(pheno17)

## 2018
dev.off()

t17 <- asreml(
  fixed = GRYLD ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno18
)
plot(t17)

h <- as.data.frame(summary(t17)$varcomp)
h

h2 <- as.data.frame(h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
h2

## 2018 all
effectvars <- names(pheno18) %in% c(
  "block", "rep", "Variety", "year", "column",
  "range", "Plot_ID"
)
traits <- colnames(pheno18[, !effectvars])
H2_2018 <- data.frame(traits)
H2_2018$Heritability <- NA
fieldInfo <- pheno18 %>%
  tidylog::select(Variety, rep, block, column, range)
ntraits <- 1:nrow(H2_2018)

for (i in ntraits) {
  print(paste("Working on trait", H2_2018[i, 1]))
  j <- H2_2018[i, 1]
  print(paste("Creating Data Frame", j))
  data <- cbind(fieldInfo, pheno18[, paste(j)])
  names(data) <- c("Variety", "rep", "block", "column", "range", "Trait")

  t17 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + rep + rep:block,
    data = data
  )
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock18_", H2_2018[i, 1], ".pdf"))
  plot(t17)
  dev.off()
  h <- as.data.frame(summary(t17)$varcomp)
  print(h)

  h2 <- (h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
  h2
  H2_2018[i, 2] <- h2
}

# lme4 comparison 2018
calcH2r <- function(dat, fill = NA, ...) {
  r <- length(table(dat$rep))

  effectvars <- names(dat) %in% c(
    "block", "rep", "Variety", "year", "column",
    "range", "Plot_ID"
  )

  t <- colnames(dat[, !effectvars])

  for (i in t) {
    print(paste("Working on trait", i))
    h <- as.data.frame(VarCorr(
      lmer(paste0(i, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)
    ))
    H2 <- h[1, 4] / (h[1, 4] + (h[4, 4] / r))
    print(H2)
  }
}

calcH2r(pheno18)

## 2019
dev.off()

t19 <- asreml(
  fixed = GRYLD ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno19
)
plot(t19)

h <- as.data.frame(summary(t19)$varcomp)
h

h2 <- as.data.frame(h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
h2

## 2019 all
effectvars <- names(pheno19) %in% c(
  "block", "rep", "Variety", "year", "column",
  "range", "Plot_ID"
)
traits <- colnames(pheno19[, !effectvars])
H2_2019 <- data.frame(traits)
H2_2019$Heritability <- NA
fieldInfo <- pheno19 %>%
  tidylog::select(Variety, rep, block, column, range)
ntraits <- 1:nrow(H2_2019)

for (i in ntraits) {
  print(paste("Working on trait", H2_2019[i, 1]))
  j <- H2_2019[i, 1]
  print(paste("Creating Data Frame", j))
  data <- cbind(fieldInfo, pheno19[, paste(j)])
  names(data) <- c("Variety", "rep", "block", "column", "range", "Trait")

  t19 <- asreml(
    fixed = Trait ~ 1,
    random = ~ Variety + rep + rep:block,
    data = data
  )
  pdf(paste0("./Figures/AsremlPlots/ASREML_repBlock19_", H2_2019[i, 1], ".pdf"))
  plot(t19)
  dev.off()
  h <- as.data.frame(summary(t19)$varcomp)
  print(h)

  h2 <- (h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
  h2
  H2_2019[i, 2] <- h2
}

# lme4 comparison 2019
calcH2r <- function(dat, fill = NA, ...) {
  r <- length(table(dat$rep))

  effectvars <- names(dat) %in% c(
    "block", "rep", "Variety", "year", "column",
    "range", "Plot_ID"
  )

  t <- colnames(dat[, !effectvars])

  for (i in t) {
    print(paste("Working on trait", i))
    h <- as.data.frame(VarCorr(
      lmer(paste0(i, "~ (1|Variety) + (1|rep) + (1|rep:block)"), data = dat)
    ))
    H2 <- h[1, 4] / (h[1, 4] + (h[4, 4] / r))
    print(H2)
  }
}

calcH2r(pheno19)

H2_2017 <- H2_2017 %>%
  separate(traits, c("Trait", "Date"), sep = "_")
H2_2017$Date <- as.Date(H2_2017$Date, format = "%Y%m%d")

H2_2017 %>%
  tidylog::filter(Trait != "GRYLD") %>%
  ggplot(aes(x = Date, y = Heritability)) +
  geom_point() +
  facet_wrap(~Trait, scales = "free") +
  geom_hline(yintercept = 0.8007513, linetype = 2, colour = "blue") +
  labs(title = "Broad-sense heritability of VI over Date 2016/2017") +
  coord_cartesian(ylim = c(0, 1))

write.table(H2_2017, "./Phenotype_Database/Hyp7_heritability17.txt",
  quote = F,
  sep = "\t", col.names = T, row.names = F
)

H2_2018 <- H2_2018 %>%
  separate(traits, c("Trait", "Date"), sep = "_")
H2_2018$Date <- as.Date(H2_2018$Date, format = "%Y%m%d")

H2_2018 %>%
  tidylog::filter(Trait != "GRYLD") %>%
  tidylog::filter(Trait != "height") %>%
  tidylog::filter(Date > as.Date("20180101", format = "%Y%m%d")) %>%
  ggplot(aes(x = Date, y = Heritability)) +
  geom_point() +
  geom_hline(yintercept = 5.502152e-01, linetype = 2, colour = "blue") +
  facet_wrap(~Trait, scales = "free") +
  scale_x_date(date_breaks = "10 days", date_labels = "%b%d") +
  labs(title = "Broad-sense heritability of VI over Date 2017/2018") +
  coord_cartesian(ylim = c(0, 1))

write.table(H2_2018, "./Phenotype_Database/Hyp7_heritability18.txt",
  quote = F,
  sep = "\t", col.names = T, row.names = F
)

H2_2019 <- H2_2019 %>%
  separate(traits, c("Trait", "Date"), sep = "_")
H2_2019$Date <- as.Date(H2_2019$Date, format = "%Y%m%d")

H2_2019 %>%
  tidylog::filter(Trait != "GRYLD") %>%
  tidylog::filter(Trait != "height") %>%
  tidylog::filter(Date > as.Date("20190301", format = "%Y%m%d")) %>%
  ggplot(aes(x = Date, y = Heritability)) +
  geom_point() +
  geom_hline(yintercept = 0.7247998, linetype = 2, colour = "blue") +
  facet_wrap(~Trait, scales = "free") +
  scale_x_date(date_breaks = "10 days", date_labels = "%b%d") +
  labs(title = "Broad-sense heritability of VI over Date 2018/2019") +
  coord_cartesian(ylim = c(0, 1))

write.table(H2_2019, "./Phenotype_Database/Hyp7_heritability19.txt",
  quote = F,
  sep = "\t", col.names = T, row.names = F
)

###### Examine weird residuals #####

t18 <- asreml(
  fixed = GNDVI_20171215 ~ 1,
  random = ~ Variety + rep + rep:block,
  data = pheno18
)
plot(t18)

residuals <- setDT(as.data.frame(t18[["residuals"]]), keep.rownames = T)
names(residuals) <- c("plot", "residual")
fits <- as.data.frame(fitted.asreml(t18))
names(fits) <- "fitted"
plotInfo <- pheno18 %>%
  tidylog::select(block, rep, range, column) %>%
  bind_cols(residuals) %>%
  bind_cols(fits)

plotInfo %>%
  ggplot(aes(x = fitted, y = residual, colour = rep)) +
  geom_point() +
  labs(
    x = "fitted", y = "Residual",
    title = "Residual plots for GNDVI_20171215"
  )

coef(t18)$random

h2 <- as.data.frame(h[3, 1] / (h[3, 1] + (h[4, 1] / 2)))
h2
