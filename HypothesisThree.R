rm(list = objects())
ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
require(lubridate)
library(car)
library(tidylog)
library(broom)
library(readxl)
library(lme4)
library(Hmisc)
library(psych)

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
      size = rel(1)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")
set.seed(1642)
# useful infos **reproducible research**
sessionInfo()

#### loading data ####
pheno <- fread("./Phenotype_Database/Pheno_Long171819.txt")

phenoGryld <- pheno %>%
  filter(!str_detect(entity_id, "19RKY7")) %>%
  filter(trait_id == "GRYLD") %>%
  filter(phenotype_value > 0) %>% # checking for strange
  filter(phenotype_value < 10) %>% # checking for strange
  tidylog::select(entity_id, phenotype_value, Variety, year) %>%
  dplyr::rename(GRYLD = phenotype_value)

phenoGryld17 <- phenoGryld %>%
  filter(year == "17")

phenoGryld18 <- phenoGryld %>%
  filter(year == "18")

phenoGryld19 <- phenoGryld %>%
  filter(year == "19")

##### 2017 HTP VI data load ####
path <- "./Phenotype_Database/2017_Ashland_AM3_traits_UAS/2017_ASH_AM_vis.xlsx"
htp17 <- path %>%
  excel_sheets() %>%
  purrr::set_names() %>%
  map(read_excel, path = path)

htp17long <- map2_df(htp17, names(htp17), ~ mutate(.x, ID = .y)) %>%
  group_by(ID) %>%
  gather(key = Date, value = value, `20170331`:`20170609`) %>%
  left_join(phenoGryld17, by = c("Plot_ID" = "entity_id")) %>%
  tidylog::select(-year)

htp17long$Date <- as.Date(htp17long$Date, "%Y%m%d")

str(htp17long)

htp17Wide <- htp17long %>%
  unite("trait_date", c("ID", "Date")) %>%
  spread(trait_date, value)

##### 2018 HTP VI data load ####

htpPheno <- c("GNDVI", "GRVI", "height", "NDRE", "NDVI", "Nir", "RE")

htpFileLoad <- function(htp, f, path, ...) {
  for (i in htp) {
    fileNames <- list.files(
      path = path,
      full.names = T,
      pattern = paste0("_", i)
    )

    traitNames <- basename(fileNames) %>%
      str_remove_all(c(".csv"))
    load.file <- function(filename) {
      d <- fread(
        file = filename, header = TRUE, check.names = F,
        data.table = F
      )
      d
    }

    data <- lapply(fileNames, load.file)
    names(data) <- traitNames
    data <- plyr::ldply(data, data.frame, .id = "Phenotype")
    print(colnames(data))
    data <- data[, 1:3]
    names(data) <- c("Phenotype", "entity_id", "phenotype_value")
    t <- tidyr::spread(
      data = data, key = "Phenotype",
      value = "phenotype_value"
    )
    head(t)
    f <- left_join(f, t)
  }
  return(f)
}

htp18Wide <- htpFileLoad(htpPheno, phenoGryld18,
  path =
    "./Phenotype_Database/2018_Ashland_AM3_traits_UAS"
)

htp18Long <- htp18Wide %>%
  gather(key = Trait, value = value, `20171120_GNDVI`:`20180613_RE`) %>%
  separate(Trait, c("Date", "ID"), sep = "_")
htp18Long$Date <- as.Date(htp18Long$Date, "%Y%m%d")

str(htp18Long)

##### 2019 HTP data load ####
htpPheno <- c("GNDVI", "NDRE", "NDVI", "Nir", "RE")
htp19Wide <- htpFileLoad(
  htp = htpPheno, f = phenoGryld19,
  path =
    "./Phenotype_Database/2019_Ashland_AM3_traits_UAS"
)

htp19Long <- htp19Wide %>%
  gather(key = Trait, value = value, `20190103_GNDVI`:`20190624_RE`) %>%
  separate(Trait, c("Date", "ID"), sep = "_")
htp19Long$Date <- as.Date(htp19Long$Date, "%Y%m%d")

#### Scatterplot of VI vs GRYLD ####

correlationPlots <- function(dat, htp, path, xaxis, value, year, ...) {
  plotList <- list()
  for (i in htp) {
    d <- dat %>%
      filter(ID == paste(i))
    thisPlot <- ggplot(
      data = d,
      mapping = aes_string(
        x = xaxis,
        y = value
      )
    ) +
      geom_point(size = 1) +
      geom_smooth(method = "lm") +
      facet_wrap(~Date, scales = "free") +
      labs(
        x = xaxis, y = paste(i),
        title = paste0(i, " vs ", xaxis, " ", year)
      )

    thisPlot
    plotList[[i]] <- thisPlot
    ggsave(paste0(i, "_", year, ".png"), thisPlot,
      path = path, width = 24,
      height = 20,
      units = "cm"
    )
    print(i)
  }
  return(plotList)
}

htpPheno <- unique(htp17long$ID)

corplots17 <- correlationPlots(
  dat = htp17long, htp = htpPheno,
  path = "./Figures/Correlations",
  xaxis = "GRYLD", value = "value", year = "2017"
)

htpPheno <- unique(htp18Long$ID)

corplots18 <- correlationPlots(
  dat = htp18Long, htp = htpPheno,
  path = "./Figures/Correlations",
  xaxis = "GRYLD", value = "value", year = "2018"
)

htpPheno <- unique(htp19Long$ID)

corplots19 <- correlationPlots(
  dat = htp19Long, htp = htpPheno,
  path = "./Figures/Correlations",
  xaxis = "GRYLD", value = "value", year = "2019"
)

#### Data correlation with significance ####

phenoMatrix17 <- htp17Wide %>%
  tidylog::select(-"Plot_ID", -"Variety")

corMat17 <- corr.test(as.matrix(phenoMatrix17),
  method = "pearson",
  adjust = "holm"
)

phenoMatrix18 <- htp18Wide %>%
  tidylog::select(-entity_id, -Variety, -year)

corMat18 <- corr.test(as.matrix(phenoMatrix18),
  method = "pearson",
  adjust = "holm"
)

phenoMatrix19 <- htp19Wide %>%
  tidylog::select(-entity_id, -Variety, -year)

corMat19 <- corr.test(as.matrix(phenoMatrix19),
  method = "pearson",
  adjust = "holm"
)

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
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

flatCor17 <- flattenCorrMatrix(corMat17$r, corMat17$p)
viGryld17cor <- flatCor17 %>%
  filter(row == "GRYLD")
Gryld17corCI <- corMat17$ci

flatCor18 <- flattenCorrMatrix(corMat18$r, corMat18$p)
viGryld18cor <- flatCor18 %>%
  filter(row == "GRYLD")
Gryld18corCI <- corMat18$ci

flatCor19 <- flattenCorrMatrix(corMat19$r, corMat19$p)
viGryld19cor <- flatCor19 %>%
  filter(row == "GRYLD")
Gryld19corCI <- corMat19$ci

#### Proportion of variance in GRYLD explained by VI by year ####

htp17long <- as_tibble(htp17long)
# check analysis
fit <- lm(htp17Wide$GRYLD ~ htp17Wide$`GNDVI_2017-03-31`, data = htp17Wide)
summary(fit)

reg17GNDVI <- htp17long %>%
  tidylog::filter(ID == "GNDVI") %>%
  tidylog::select(-Variety, -Plot_ID, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "GNDVI") %>%
  unite("VI_date", c("VI", "Date"))

reg17GRVI <- htp17long %>%
  tidylog::filter(ID == "GRVI") %>%
  tidylog::select(-Variety, -Plot_ID, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "GRVI") %>%
  unite("VI_date", c("VI", "Date"))

reg17NDVI <- htp17long %>%
  tidylog::filter(ID == "NDVI") %>%
  tidylog::select(-Variety, -Plot_ID, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NDVI") %>%
  unite("VI_date", c("VI", "Date"))

reg17NDRE <- htp17long %>%
  tidylog::filter(ID == "NDRE") %>%
  tidylog::select(-Variety, -Plot_ID, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NDRE") %>%
  unite("VI_date", c("VI", "Date"))

reg17NIR <- htp17long %>%
  tidylog::filter(ID == "NIR") %>%
  tidylog::select(-Variety, -Plot_ID, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NIR") %>%
  unite("VI_date", c("VI", "Date"))

reg17RE <- htp17long %>%
  tidylog::filter(ID == "RedEdge") %>%
  tidylog::select(-Variety, -Plot_ID, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "RE") %>%
  unite("VI_date", c("VI", "Date"))

htp18Long <- as_tibble(htp18Long)
# check analysis
fit <- lm(htp18Wide$GRYLD ~ htp18Wide$`20171120_GNDVI`)
summary(fit)

reg18GNDVI <- htp18Long %>%
  tidylog::filter(ID == "GNDVI") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "GNDVI") %>%
  unite("VI_date", c("VI", "Date"))

reg18GRVI <- htp18Long %>%
  tidylog::filter(ID == "GRVI") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "GRVI") %>%
  unite("VI_date", c("VI", "Date"))

reg18NDVI <- htp18Long %>%
  tidylog::filter(ID == "NDVI") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NDVI") %>%
  unite("VI_date", c("VI", "Date"))

reg18NDRE <- htp18Long %>%
  tidylog::filter(ID == "NDRE") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NDRE") %>%
  unite("VI_date", c("VI", "Date"))

reg18NIR <- htp18Long %>%
  tidylog::filter(ID == "Nir") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NIR") %>%
  unite("VI_date", c("VI", "Date"))

reg18RE <- htp18Long %>%
  tidylog::filter(ID == "RE") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "RE") %>%
  unite("VI_date", c("VI", "Date"))

htp19Long <- as_tibble(htp19Long)
# check analysis
fit <- lm(htp19Wide$GRYLD ~ htp19Wide$`20190103_GNDVI`)
summary(fit)

reg19GNDVI <- htp19Long %>%
  tidylog::filter(ID == "GNDVI") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "GNDVI") %>%
  unite("VI_date", c("VI", "Date"))

reg19GRVI <- htp19Long %>%
  tidylog::filter(ID == "GRVI") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "GRVI") %>%
  unite("VI_date", c("VI", "Date"))

reg19NDVI <- htp19Long %>%
  tidylog::filter(ID == "NDVI") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NDVI") %>%
  unite("VI_date", c("VI", "Date"))

reg19NDRE <- htp19Long %>%
  tidylog::filter(ID == "NDRE") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NDRE") %>%
  unite("VI_date", c("VI", "Date"))

reg19NIR <- htp19Long %>%
  tidylog::filter(ID == "Nir") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "NIR") %>%
  unite("VI_date", c("VI", "Date"))

reg19RE <- htp19Long %>%
  tidylog::filter(ID == "RE") %>%
  tidylog::select(-Variety, -entity_id, -ID) %>%
  nest(data = c(value, GRYLD)) %>%
  mutate(
    test = map(data, ~ lm(.x$GRYLD ~ .x$value)),
    glanced = map(test, glance),
    tidied = map(test, tidy)
  ) %>%
  unnest(c(glanced)) %>%
  rename(
    F_statistic = statistic,
    Ftest_pvalue = p.value
  ) %>%
  unnest(c(tidied)) %>%
  select(-data, -test) %>%
  filter(term != "(Intercept)") %>%
  add_column(VI = "RE") %>%
  unite("VI_date", c("VI", "Date"))

VarGryldVI17 <- bind_rows(
  reg17GNDVI, reg17GRVI, reg17NDRE, reg17NDVI, reg17NIR,
  reg17RE
) %>%
  separate(VI_date, c("VI", "Date"), sep = "_") %>%
  glimpse()

VarGryldVI18 <- bind_rows(
  reg18GNDVI, reg18GRVI, reg18NDRE, reg18NDVI,
  reg18NIR, reg18RE
) %>%
  separate(VI_date, c("VI", "Date"), sep = "_") %>%
  glimpse()

VarGryldVI19 <- bind_rows(
  reg19GNDVI, reg19NDRE, reg19NDVI,
  reg19NIR, reg19RE
) %>%
  separate(VI_date, c("VI", "Date"), sep = "_") %>%
  glimpse()

VarGryldVI17$Date <- as.Date(VarGryldVI17$Date)

VarGryldVI18$Date <- as.Date(VarGryldVI18$Date)

VarGryldVI19$Date <- as.Date(VarGryldVI19$Date)

VarGryldVI17 %>%
  ggplot(aes(x = Date, y = adj.r.squared, colour = p.value)) +
  geom_point() +
  facet_wrap(~VI, scales = "free") +
  scale_color_gradient(low = "#e41a1c", high = "#000000") +
  theme(plot.subtitle = element_text(size = rel(1.75))) +
  labs(
    title = bquote(R^2 ~ " for linear regression models explaining GRYLD"),
    subtitle = "2016/2017 season"
  )

VarGryldVI18 %>%
  ggplot(aes(x = Date, y = adj.r.squared, colour = p.value)) +
  geom_point() +
  facet_wrap(~VI, scales = "free") +
  scale_color_gradient(low = "#e41a1c", high = "#000000") +
  theme(plot.subtitle = element_text(size = rel(1.75))) +
  labs(
    title = bquote(R^2 ~ " for linear regression models explaining GRYLD"),
    subtitle = "2017/2018 season"
  )

VarGryldVI19 %>%
  ggplot(aes(x = Date, y = adj.r.squared, colour = p.value)) +
  geom_point() +
  facet_wrap(~VI, scales = "free") +
  scale_color_gradient(low = "#e41a1c", high = "#000000") +
  theme(plot.subtitle = element_text(size = rel(1.75))) +
  labs(
    title = bquote(R^2 ~ " for linear regression models explaining GRYLD"),
    subtitle = "2018/2019 season"
  )

write.table(VarGryldVI17,
  "./Phenotype_Database/linearRegression_VIbyDate17.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(VarGryldVI18,
  "./Phenotype_Database/linearRegression_VIbyDate18.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(VarGryldVI19,
  "./Phenotype_Database/linearRegression_VIbyDate19.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(phenoMatrix17, "./Phenotype_Database/phenoMatrix17.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(phenoMatrix18, "./Phenotype_Database/phenoMatrix18.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(phenoMatrix19, "./Phenotype_Database/phenoMatrix19.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(htp18Long, "./Phenotype_Database/pheno18_htpLong.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(htp17long, "./Phenotype_Database/pheno17_htpLong.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(htp19Long, "./Phenotype_Database/pheno19_htpLong.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
