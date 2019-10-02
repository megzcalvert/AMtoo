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
pheno17 <- pheno17 %>%
  tidylog::select(-GRYLD) %>%
  mutate(Date = as.Date(Date))

colnames(pheno18)
pheno18 <- pheno18 %>%
  tidylog::select(-GRYLD, -year) %>%
  filter(ID != "height") %>%
  mutate(Date = as.Date(Date))

colnames(pheno19)
pheno19 <- pheno19 %>%
  tidylog::select(-GRYLD, -year) %>%
  filter(ID != "height") %>%
  mutate(Date = as.Date(Date))

pheno17 %>%
  ggplot(aes(x = Date, y = value, colour = Plot_ID)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  labs(title = "VI by line over time 2016/2017") +
  guides(colour = FALSE)

pheno18 %>%
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  labs(title = "VI by line over time 2017/2018") +
  guides(colour = FALSE)

pheno18 %>%
  filter(Date > "2018-04-1") %>%
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  labs(title = "VI by line over time 2017/2018") +
  guides(colour = FALSE)

pheno19 %>%
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  labs(title = "VI by line over time 2018/2019") +
  guides(colour = FALSE)

pheno19 %>%
  filter(Date > "2019-04-1") %>%
  ggplot(aes(x = Date, y = value, colour = entity_id)) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  labs(title = "VI by line over time 2018/2019") +
  guides(colour = FALSE)

# Anova of linear reg
pheno17$Date <- as.factor(pheno17$Date)
nested17 <- pheno17 %>%
  tidylog::select(-Plot_ID) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested17, "./Phenotype_Database/ANOVA_VIbyDateVariety17.txt",
  quote = F, row.names = F, col.names = T, sep = "\t"
)

pheno18$Date <- as.factor(pheno18$Date)
nested18All <- pheno18 %>%
  tidylog::select(-entity_id) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested18All, "./Phenotype_Database/ANOVA_VIbyDateVariety18.txt",
  quote = F, row.names = F, col.names = T, sep = "\t"
)

nested18AfterV <- pheno18 %>%
  tidylog::select(-entity_id) %>%
  filter(Date != "2017-11-20") %>%
  filter(Date != "2017-11-27") %>%
  filter(Date != "2017-12-05") %>%
  filter(Date != "2017-12-15") %>%
  filter(Date != "2017-12-18") %>%
  filter(ID != "height") %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested18AfterV,
  "./Phenotype_Database/ANOVA_VIbyDateVariety18_afterVern.txt",
  quote = F, row.names = F, col.names = T, sep = "\t"
)

pheno19$Date <- as.factor(pheno19$Date)
nested19All <- pheno19 %>%
  tidylog::select(-entity_id) %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested19All, "./Phenotype_Database/ANOVA_VIbyDateVariety19.txt",
  quote = F, row.names = F, col.names = T, sep = "\t"
)

nested19AfterV <- pheno19 %>%
  tidylog::select(-entity_id) %>%
  filter(Date != "2019-01-03") %>%
  group_by(ID) %>%
  do(tidy(anova(lm(value ~ Date + Variety + Date:Variety, data = .))))

write.table(nested19AfterV,
  "./Phenotype_Database/ANOVA_VIbyDateVariety19_afterVern.txt",
  quote = F, row.names = F, col.names = T, sep = "\t"
)

nested17 %>%
  tidylog::filter(ID != "RedEdge") %>%
  tidylog::filter(ID != "NIR") %>%
  ggplot(aes(x = p.value, fill = ID)) +
  geom_histogram()
