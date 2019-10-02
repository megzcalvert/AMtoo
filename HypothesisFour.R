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
      size = rel(1)
    ),
    complete = F
  )

theme_set(custom_theme)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")
set.seed(1642)

#### Loading data ####
pheno17 <- fread("./Phenotype_Database/pheno17_htpLong.txt")
pheno18 <- fread("./Phenotype_Database/pheno18_htpLong.txt")
pheno19 <- fread("./Phenotype_Database/pheno19_htpLong.txt")

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
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%d%b"
  ) +
  labs(title = "Correlation with CI 2017") +
  ylab("Pearson correlation co-efficient") 

nested18 %>%
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1) +
  geom_hline(
    yintercept = 0, linetype = 2,
    colour = "darkgrey"
  ) +
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%d%b"
  ) +
  labs(title = "Correlation with CI 2018") +
  ylab("Pearson correlation co-efficient") 

nested19 %>%
  ggplot(aes(x = Date, y = estimate, color = ID)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = 1) +
  geom_hline(
    yintercept = 0, linetype = 2,
    colour = "darkgrey"
  ) +
  scale_x_date(
    date_breaks = "1 week",
    date_labels = "%d%b"
  ) +
  labs(title = "Correlation with CI 2019") +
  ylab("Pearson correlation co-efficient") 

write.table(nested17, "./Phenotype_Database/Correlation_VI_2017.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(nested18, "./Phenotype_Database/Correlation_VI_2018.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
write.table(nested19, "./Phenotype_Database/Correlation_VI_2019.txt",
  sep = "\t", quote = F, row.names = F, col.names = T
)
