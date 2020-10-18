rm(list = objects())
ls()

library(readr)
library(ggbiplot)
library(data.table)
library(tidyverse)
library(janitor)
library(GGally)
require(lubridate)
library(car)
library(tidylog)
library(broom)
library(readxl)
library(gvlma)

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

#### Load data ####

hddt17 <- fread("./Phenotype_Database/HDDT2017.txt")
hddt17$hddt17 <- as.Date(hddt17$hddt17)
hddt18 <- fread("./Phenotype_Database/HDDT2018.txt")
hddt18$hddt18 <- as.Date(hddt18$hddt18)
hddt19 <- fread("./Phenotype_Database/HDDT2019.txt")
hddt19$hddt19 <- as.Date(hddt19$hddt19)

htpPheno <- c("GNDVI", "GRVI", "height", "NDRE", "NDVI", "Nir", "RE")

#####################################################################
#### 2017 HTP vegetaion indices over time and relation to hddt ####

path <- "./Phenotype_Database/2017_Ashland_AM3_traits_UAS/2017_ASH_AM_vis.xlsx"
htp17 <- path %>%
  excel_sheets() %>%
  purrr::set_names() %>%
  map(read_excel, path = path)

htp17long <- map2_df(htp17, names(htp17), ~ mutate(.x, ID = .y)) %>%
  group_by(ID) %>%
  gather(key = Date, value = value, `20170331`:`20170609`) %>%
  mutate(Date = as.Date(Date, "%Y%m%d")) %>%
  glimpse()

htp17long %>%
  ggplot(aes(x = Date, y = value, colour = Plot_ID)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  scale_x_date(breaks = "1 week", date_labels = "%b %d") +
  geom_vline(
    alpha = 0.5, colour = "#1d91c0",
    aes(xintercept = hddt17), hddt17
  ) +
  labs(title = "2016/2017 season VI with HDDT") +
  guides(colour = FALSE)

#####################################################################
#### 2018 HTP vegetaion indices over time and relation to hddt ####

htpFileLoad <- function(htp, f, ...) {
  for (i in htp) {
    fileNames <- list.files(
      path = "./Phenotype_Database/2018_Ashland_AM3_traits_UAS",
      full.names = T,
      pattern = paste0("_", i)
    )

    traitNames <- basename(fileNames) %>%
      str_remove_all(c(".csv"))
    load.file <- function(filename) {
      d <- fread(file = filename, header = TRUE, check.names = F, data.table = F)
      d
    }

    data <- lapply(fileNames, load.file)
    names(data) <- traitNames
    data <- plyr::ldply(data, data.frame, .id = "Phenotype")
    print(colnames(data))
    t <- data %>%
      pivot_wider(
        id_cols = Plot_ID,
        names_from = Phenotype,
        values_from = paste0(colnames(data)[3])
      )

    head(t)
    f <- left_join(f, t, by = c("plots18" = "Plot_ID"))
  }
  return(f)
}

htp18 <- htpFileLoad(htpPheno, hddt18)

htp18Long <- htp18 %>%
  gather(key = Trait, value = value, `20171120_GNDVI`:`20180613_RE`) %>%
  separate(Trait, c("Date", "ID"), sep = "_") %>%
  mutate(Date = as.Date(Date, "%Y%m%d"))

htp18Long %>%
  tidylog::filter(ID != "height") %>%
  ggplot(aes(x = Date, y = value, colour = plots18)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  scale_x_date(breaks = "1 week", date_labels = "%b %d") +
  geom_vline(
    alpha = 0.5, colour = "#1d91c0",
    aes(xintercept = hddt18), hddt18
  ) +
  guides(colour = FALSE)

htp18Long %>%
  filter(Date >= "2018-04-01") %>%
  tidylog::filter(ID != "height") %>%
  ggplot(aes(x = Date, y = value, colour = plots18)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  scale_x_date(breaks = "1 week", date_labels = "%b %d") +
  geom_vline(
    alpha = 0.5, colour = "#1d91c0",
    aes(xintercept = hddt18), hddt18
  ) +
  labs(title = "2017/2018 season VI with HDDT") +
  guides(colour = FALSE)

##### 2019 HTP data load ####

htpPheno <- c("GNDVI", "NDRE", "NDVI", "Nir", "RE")

htpFileLoad <- function(htp, f, ...) {
  for (i in htp) {
    fileNames <- list.files(
      path = "./Phenotype_Database/2019_Ashland_AM3_traits_UAS",
      full.names = T,
      pattern = paste0("_", i)
    )

    traitNames <- basename(fileNames) %>%
      str_remove_all(c(".csv"))
    load.file <- function(filename) {
      d <- fread(file = filename, header = TRUE, check.names = F, data.table = F)
      d
    }

    data <- lapply(fileNames, load.file)
    names(data) <- traitNames
    data <- plyr::ldply(data, data.frame, .id = "Phenotype")
    print(colnames(data))
    t <- data %>%
      pivot_wider(
        id_cols = Plot_ID,
        names_from = Phenotype,
        values_from = paste0(colnames(data)[3])
      )
    head(t)
    f <- left_join(f, t, by = c("plots19" = "Plot_ID"))
  }
  return(f)
}

htp19Wide <- htpFileLoad(htpPheno, hddt19,
  path =
    "./Phenotype_Database/2019_Ashland_AM3_traits_UAS"
)

htp19Long <- htp19Wide %>%
  gather(key = Trait, value = value, `20190103_GNDVI`:`20190624_RE`) %>%
  separate(Trait, c("Date", "ID"), sep = "_")
htp19Long$Date <- as.Date(htp19Long$Date, "%Y%m%d")

htp19Long %>%
  filter(Date >= "2019-04-01") %>%
  tidylog::filter(ID != "height") %>%
  ggplot(aes(x = Date, y = value, colour = plots19)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ID, scales = "free", ncol = 2) +
  scale_x_date(breaks = "1 week", date_labels = "%b %d") +
  geom_vline(
    alpha = 0.5, colour = "#1d91c0",
    aes(xintercept = hddt19), hddt19
  ) +
  labs(title = "2018/2019 season VI with HDDT") +
  guides(colour = FALSE)

