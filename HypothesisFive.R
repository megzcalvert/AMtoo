rm(list = objects())
ls()

library(readr)
library(data.table)
library(tidyverse)
library(janitor)
library(tidylog)
library(broom)
library(Hmisc)
library(psych)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")

pheno17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")
hddt17 <- fread("./Phenotype_Database/HDDT2017.txt")
hddt17$hddt17 <- as.Date(hddt17$hddt17)
hddt18 <- fread("./Phenotype_Database/HDDT2018.txt")
hddt18$hddt18 <- as.Date(hddt18$hddt18)
hddt19 <- fread("./Phenotype_Database/HDDT2019.txt")
hddt19$hddt19 <- as.Date(hddt19$hddt19)

colnames(pheno17)

# Check
trial <- pheno17 %>%
  filter(ID == "NDVI") %>%
  filter(Date == "2017-03-31")
cor.test(trial$value, trial$GRYLD)

nested17 <- pheno17 %>%
  tidylog::select(-Plot_ID, -Variety) %>%
  group_by(Date, ID) %>%
  nest() %>%
  mutate(
    correlation = map(data, ~ cor.test(.x$GRYLD, .x$value)),
    tidyCor = map(correlation, glance)
  ) %>%
  unnest(tidyCor) %>%
  tidylog::select(-data, -correlation)

nested18 <- pheno18 %>%
  filter(ID != "height") %>%
  tidylog::select(-entity_id, -Variety, -year) %>%
  group_by(Date, ID) %>%
  nest() %>%
  mutate(
    correlation = map(data, ~ cor.test(.x$GRYLD, .x$value)),
    tidyCor = map(correlation, glance)
  ) %>%
  unnest(tidyCor) %>%
  tidylog::select(-data, -correlation)

nested19 <- pheno19 %>%
  filter(ID != "height") %>%
  tidylog::select(-entity_id, -Variety, -year) %>%
  group_by(Date, ID) %>%
  nest() %>%
  mutate(
    correlation = map(data, ~ cor.test(.x$GRYLD, .x$value)),
    tidyCor = map(correlation, glance)
  ) %>%
  unnest(tidyCor) %>%
  tidylog::select(-data, -correlation)

nested17$Date <- as.Date(nested17$Date)
nested18$Date <- as.Date(nested18$Date)
nested19$Date <- as.Date(nested19$Date)

nested17 %>%
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1) +
  geom_hline(
    yintercept = 0, linetype = 2,
    colour = "darkgrey"
  ) +
  geom_vline(aes(xintercept = hddt17), hddt17,
    colour = "#74c476",
    alpha = 0.75
  ) +
  theme_bw() +
  # scale_colour_manual(values = c("#e41a1c","#377eb8","#4daf4a",
  #                               "#984ea3","#ff7f00","#525252")) +
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%d%b"
  ) +
  labs(title = "Correlation with CI 2017") +
  ylab("Pearson correlation co-efficient") +
  theme(
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 16),
    title = element_text(size = 20),
    # legend.position = "bottom",
    legend.text = element_text(size = 14)
  )

nested18 %>%
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1) +
  geom_hline(
    yintercept = 0, linetype = 2,
    colour = "darkgrey"
  ) +
  geom_vline(aes(xintercept = hddt18), hddt18,
    colour = "#74c476",
    alpha = 0.75
  ) +
  theme_bw() +
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%m/%d"
  ) +
  labs(title = "Correlation with CI 2018") +
  ylab("Pearson correlation co-efficient") +
  theme(
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 16),
    title = element_text(size = 20),
    # legend.position = "bottom",
    legend.text = element_text(size = 14)
  )

nested19 %>%
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1) +
  geom_hline(
    yintercept = 0, linetype = 2,
    colour = "darkgrey"
  ) +
  geom_vline(aes(xintercept = hddt19), hddt19,
    colour = "#74c476",
    alpha = 0.75
  ) +
  theme_bw() +
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%m/%d"
  ) +
  labs(title = "Correlation with CI 2019") +
  ylab("Pearson correlation co-efficient") +
  theme(
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 16),
    title = element_text(size = 20),
    # legend.position = "bottom",
    legend.text = element_text(size = 14)
  )
