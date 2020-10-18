rm(list = objects())
ls()

library(tidyverse)
library(tidylog)
library(data.table)
library(sysfonts)
library(ggpubr)
library(cowplot)
library(janitor)
library(extrafont)
library(sysfonts)
library(GGally)
library(ggrepel)

loadfonts()
fonts()

font_paths()
font_families()
fontFil <- font_files()
fontTab <- fonttable()

pdfFonts()
pstFont <- postscriptFonts()
font_add(
  family = "Times New Roman",
  regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
)
font_add(
  family = "Arial",
  regular = "/System/Library/Fonts/Supplemental/Arial.ttf"
)
font_families()

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
      size = rel(2)
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
    legend.title = element_text(size = rel(2)),
    panel.grid.major = element_line(
      colour = "#949494",
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
      size = rel(1.75),
      hjust = 0,
      vjust = 1
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 1
    ),
    strip.text = element_text(
      colour = "black",
      size = rel(1.75)
    ),
    text = element_text(family = "Arial"),
    complete = F
  )

theme_set(custom_theme)

getwd()
setwd("~/Dropbox/Research_Poland_Lab/AM Panel/")
set.seed(1642)
###############################################################################
#### Heading date ####
hddt17 <- fread("./Phenotype_Database/HDDT2017.txt")
hddt18 <- fread("./Phenotype_Database/HDDT2018.txt")
hddt19 <- fread("./Phenotype_Database/HDDT2019.txt")

hddt17 <- hddt17 %>%
  tidylog::mutate(hddt17 = as.Date(hddt17))
hddt18 <- hddt18 %>%
  tidylog::mutate(hddt18 = as.Date(hddt18))
hddt19 <- hddt19 %>%
  tidylog::mutate(hddt19 = as.Date(hddt19))

#### Heritability ####

heritability17 <- fread("./Phenotype_Database/Hyp7_heritability17.txt")
heritability18 <- fread("./Phenotype_Database/Hyp7_heritability18.txt")
heritability19 <- fread("./Phenotype_Database/Hyp7_heritability19.txt")

herGryld17 <- heritability17 %>%
  filter(Trait == "GRYLD") %>% 
  tidylog::mutate(Year = "17")
heritability17 <- heritability17 %>%
  mutate(Date = as.Date(heritability17$Date, format = "%Y-%m-%d"),
         Year = "17",
         hddt_min = as.Date("2017-04-24"),
         hddt_max = as.Date("2017-05-15")) %>%
  filter(Trait != "GRYLD") %>%
  filter(Trait != "GRVI") %>%
  filter(Trait != "NIR") %>%
  filter(Trait != "RedEdge") %>%
  filter(Trait != "numHddt17")

her17 <- heritability17 %>%
  ggplot(aes(
    x = Date,
    y = Heritability
  )) +
  geom_rect(aes(
    xmin = as.Date("2017-04-24"),
    xmax = as.Date("2017-05-15"),
    ymax = 1,
    ymin = 0
  ),
  colour = "#bdbdbd",
  fill = "#bdbdbd"
  ) +
  geom_point() +
  geom_hline(aes(yintercept = herGryld17$Heritability),
    linetype = 2
  ) +
  facet_wrap(~Trait, scales = "free", ncol = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(
    subtitle = "2016/2017 Season"
  )
her17

ggsave("./Figures/Figure3_heritability17.png",
  height = 10,
  width = 32.5, units = "cm", dpi = 350
)

herGryld18 <- heritability18 %>%
  filter(Trait == "GRYLD") %>% 
  tidylog::mutate(Year = "18")
heritability18 <- heritability18 %>%
  mutate(Date = as.Date(heritability18$Date, format = "%Y-%m-%d"),
         Year = "18",
         hddt_max = as.Date("2018-05-11"),
         hddt_min = as.Date("2018-05-25")) %>%
  filter(Date > "2018-02-01") %>%
  filter(Trait != "GRYLD") %>%
  filter(Trait != "GRVI") %>%
  filter(Trait != "Nir") %>%
  filter(Trait != "RE") %>%
  filter(Trait != "height") %>%
  filter(Trait != "numHddt18")

her18 <- heritability18 %>%
  ggplot(aes(
    x = Date,
    y = Heritability
  )) +
  geom_rect(aes(
    xmin = as.Date("2018-05-11"),
    xmax = as.Date("2018-05-25"),
    ymax = 1,
    ymin = 0
  ),
  colour = "#bdbdbd",
  fill = "#bdbdbd"
  ) +
  geom_point() +
  geom_hline(aes(yintercept = herGryld18$Heritability),
    linetype = 2
  ) +
  facet_wrap(~Trait, scales = "free", ncol = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(
    subtitle = "2017/2018 Season"
  )
her18

ggsave("./Figures/Figure3_heritability18.png",
  height = 10,
  width = 32.5, units = "cm", dpi = 350
)

herGryld19 <- heritability19 %>%
  filter(Trait == "GRYLD") %>% 
  tidylog::mutate(Year = "19")
heritability19 <- heritability19 %>%
  mutate(Date = as.Date(heritability19$Date, format = "%Y-%m-%d"),
         Year = "19",
         hddt_min = as.Date("2019-05-08"),
         hddt_max = as.Date("2019-05-27")) %>%
  filter(Date > "2019-02-01") %>%
  filter(Trait != "GRYLD") %>%
  filter(Trait != "Nir") %>%
  filter(Trait != "RE") %>%
  filter(Trait != "height") %>%
  filter(Trait != "numHddt19")

her19 <- heritability19 %>%
  ggplot(aes(
    x = Date,
    y = Heritability
  )) +
  geom_rect(aes(
    xmin = as.Date("2019-05-08"),
    xmax = as.Date("2019-05-27"),
    ymax = 1,
    ymin = 0
  ),
  colour = "#bdbdbd",
  fill = "#bdbdbd"
  ) +
  geom_point() +
  geom_hline(aes(yintercept = herGryld19$Heritability),
    linetype = 2
  ) +
  facet_wrap(~Trait, scales = "free", ncol = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  ) +
  labs(
    subtitle = "2018/2019 Season"
  )
her19

ggsave("./Figures/Figure3_heritability19.png",
  height = 10,
  width = 32.5, units = "cm", dpi = 350
)

ggarrange(her17, her18, her19, nrow = 3, ncol = 1)

ggsave("./AMPanel_Manuscript/Figures/Figure3_heritability.jpg",
  height = 27.5,
  width = 35, units = "cm", dpi = 350
)

her_all<- bind_rows(heritability17, heritability18, heritability19)
her_all %>% 
  ggplot(aes(
    x = Date,
    y = Heritability
  )) +
  geom_rect(aes(
    xmin = hddt_min,
    xmax = hddt_max,
    ymax = 1,
    ymin = 0
  ),
  colour = "#99d8c9",
  fill = "#e5f5f9"
  ) +
  geom_point(size = 3) +
  geom_hline(aes(yintercept = herGryld19$Heritability),
             linetype = 2
  ) +
  facet_wrap(Year~Trait, scales = "free", ncol = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_date(
    date_labels = "%m/%d",
    date_breaks = "2 weeks"
  )

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/AMpanel/Figures/Figure3_heritability.png",
       height = 27.5,
       width = 42.5, units = "cm", dpi = 320
)


#### Correlations and regressions ####
cor17 <- fread("./Phenotype_Database/Correlation_VI_2017.txt")

cor18 <- fread("./Phenotype_Database/Correlation_VI_2018.txt")

cor19 <- fread("./Phenotype_Database/Correlation_VI_2019.txt")

cor17 <- cor17 %>%
  mutate(Date = as.Date(cor17$Date, format = "%Y-%m-%d"),
         Year = "17",
         hddt_min = as.Date("2017-04-24"),
         hddt_max = as.Date("2017-05-15")) %>%
  filter(ID != "GRVI") %>%
  filter(ID != "NIR") %>%
  filter(ID != "RE") %>%
  filter(ID != "RedEdge")

corLin17 <- ggplot(
  data = cor17,
  aes(x = Date, y = estimate)
) +
  geom_rect(aes(
    xmin = as.Date("2017-04-24"),
    xmax = as.Date("2017-05-15"),
    ymax = 1,
    ymin = -1
  ),
  colour = "#bdbdbd",
  fill = "#bdbdbd"
  ) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~ID, ncol = 3, scales = "free") +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d") +
  labs(
    y = "Correlation",
    subtitle = "2016/2017 Season"
  )
corLin17

ggsave("./Figures/Figure2_heritability17.png",
  height = 10,
  width = 32.5, units = "cm", dpi = 350
)

cor18 <- cor18 %>%
  mutate(Date = as.Date(cor18$Date, format = "%Y-%m-%d"),
         Year = "18",
         hddt_max = as.Date("2018-05-11"),
         hddt_min = as.Date("2018-05-25")) %>%
  filter(Date > "2018-02-02") %>%
  filter(ID != "GRVI") %>%
  filter(ID != "NIR") %>%
  filter(ID != "Nir") %>%
  filter(ID != "RE") %>%
  filter(ID != "RedEdge")

corLin18 <- ggplot(
  data = cor18,
  aes(x = Date, y = estimate)
) +
  geom_rect(aes(
    xmin = as.Date("2018-05-11"),
    xmax = as.Date("2018-05-25"),
    ymax = 1,
    ymin = -1
  ),
  colour = "#bdbdbd",
  fill = "#bdbdbd"
  ) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~ID, ncol = 3, scales = "free") +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d") +
  labs(
    y = "Correlation",
    subtitle = "2017/2018 Season"
  )
corLin18

ggsave("./Figures/Figure2_heritability18.png",
  height = 10,
  width = 32.5, units = "cm", dpi = 350
)

cor19 <- cor19 %>%
  mutate(Date = as.Date(cor19$Date, format = "%Y-%m-%d"),
         Year = "19",
         hddt_min = as.Date("2019-05-08"),
         hddt_max = as.Date("2019-05-27")) %>%
  filter(Date > "2019-02-02") %>%
  filter(ID != "GRVI") %>%
  filter(ID != "NIR") %>%
  filter(ID != "Nir") %>%
  filter(ID != "RE") %>%
  filter(ID != "RedEdge")

corLin19 <- ggplot(
  data = cor19,
  aes(x = Date, y = estimate)
) +
  geom_rect(aes(
    xmin = as.Date("2019-05-08"),
    xmax = as.Date("2019-05-27"),
    ymax = 1,
    ymin = -1
  ),
  colour = "#bdbdbd",
  fill = "#bdbdbd"
  ) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~ID, ncol = 3, scales = "free") +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d") +
  labs(
    y = "Correlation",
    subtitle = "2018/2019 Season"
  )
corLin19

ggsave("./Figures/Figure2_heritability19.png",
  height = 10,
  width = 32.5, units = "cm", dpi = 350
)

ggarrange(corLin17, corLin18, corLin19,
  nrow = 3, ncol = 1, legend = "right",
  common.legend = TRUE
)

ggsave("./AMPanel_Manuscript/Figures/Figure2_Correlation.jpg",
  height = 27.5,
  width = 35, units = "cm", dpi = 350
)

cor_all<- bind_rows(cor17, cor18, cor19)

ggplot(
  data = cor_all,
  aes(x = Date, y = estimate)
) +
  geom_rect(aes(
    xmin = hddt_min,
    xmax = hddt_max,
    ymax = 1,
    ymin = -1
  ),
  colour = "#99d8c9",
  fill = "#e5f5f9"
  ) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) +
  facet_wrap(Year~ID, ncol = 3, scales = "free") +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d") +
  labs(y = "Correlation")

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/AMpanel/Figures/Figure2.png",
       height = 27.5,
       width = 42.5, units = "cm", dpi = 320
)

#### GRYLD vs HDDT ####
pheno <- fread("./Phenotype_Database/Pheno_Long171819.txt")

pheno17 <- pheno %>%
  filter(year == "17") %>%
  filter(trait_id == "GRYLD") %>%
  tidylog::mutate(phenotype_date = as.Date(phenotype_date)) %>%
  left_join(hddt17, by = c("entity_id" = "plots17")) %>% 
  rename(hddt = hddt17)

pheno18 <- pheno %>%
  filter(year == "18") %>%
  filter(trait_id == "GRYLD") %>%
  tidylog::mutate(phenotype_date = as.Date(phenotype_date)) %>%
  left_join(hddt18, by = c("entity_id" = "plots18")) %>% 
  rename(hddt = hddt18)

pheno19 <- pheno %>%
  filter(year == "19") %>%
  filter(trait_id == "GRYLD") %>%
  tidylog::mutate(phenotype_date = as.Date(phenotype_date)) %>%
  left_join(hddt19, by = c("entity_id" = "plots19")) %>% 
  rename(hddt = hddt19)

headGr17 <- ggplot(
  data = pheno17,
  mapping = aes(
    x = hddt,
    y = phenotype_value
  )
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  scale_x_date(date_breaks = "3 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 10)) +
  annotate(geom = "text", x = as.Date("2017-04-24"),
           y = 9.6,
           label = "r = -0.399",
           size = 6,
           hjust = 0) +
  labs(
    subtitle = "2016/2017 Season",
    x = "Date",
    y = "Grain Yield (t/ha)"
  )

headGr17

headGr18 <- ggplot(
  data = pheno18,
  mapping = aes(
    x = hddt,
    y = phenotype_value
  )
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  scale_x_date(date_breaks = "3 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 10)) +
  annotate(geom = "text", x = as.Date("2018-05-10"),
           y = 9.6,
           label = "r = -0.053",
           size = 6,
           hjust = 0) +
  labs(
    subtitle = "2017/2018 Season",
    x = "Date",
    y = "Grain Yield (t/ha)"
  )

headGr18

headGr19 <- ggplot(
  data = pheno19,
  mapping = aes(
    x = hddt,
    y = phenotype_value
  )
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  scale_x_date(date_breaks = "3 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 10)) +
  annotate(geom = "text", x = as.Date("2019-05-08"),
           y = 9.6,
           label = "r = -0.491",
           size = 6,
           hjust = 0) +
  labs(
    subtitle = "2018/2019 Season",
    x = "Date",
    y = "Grain Yield (t/ha)"
  )

headGr19

ggarrange(headGr17,headGr18,headGr19, ncol = 1)
ggsave("./AMPanel_Manuscript/Figures/Figure1_hddtVSgryld.jpg",
       height = 35,
       width = 20, units = "cm", dpi = 350
)

head_all<- bind_rows(pheno17,pheno18,pheno19)

ggplot(
  data = head_all,
  mapping = aes(
    x = hddt,
    y = phenotype_value
  )
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  scale_x_date(date_breaks = "3 days", date_labels = "%m/%d") +
  coord_cartesian(ylim = c(0, 10)) +
  facet_wrap(~year, ncol = 2,
             scales = "free_x") +
  stat_cor(method = "pearson",
           label.y = 7.75, size = 6) +
  labs(x = "Date",
       y = "Grain yield (t/ha)")

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/AMpanel/Figures/Figure1.png",
       height = 35,
       width = 20, units = "cm", dpi = 350
)

#### PCA for pop struc ####
phenoLines<- pheno %>% 
  select(Variety) %>% 
  distinct()

snpChip <- fread(
  file = "./Genotype_Database/SelectedImputedBeagleNumeric_outliers.txt",
  header = TRUE, check.names = F, sep = "\t"
)

chrSum <- plyr::count(snpChip, vars = "chrom")
snpChip <- snpChip %>%
  tidylog::filter(chrom != "UN") %>%
  tidylog::mutate(akron = as.numeric(akron))

snpMatrix <- t(snpChip[, c(-1, -2, -3)])

snpMatrix<- setDT(as.data.frame(snpMatrix),keep.rownames = TRUE) %>% 
  semi_join(phenoLines, by = c("rn" = "Variety")) %>% 
  column_to_rownames(var = "rn") 
snpMatrix<- as.matrix(snpMatrix)

pcaMethods::checkData(snpMatrix) # Check PCA assumptions

pcaAM <- pcaMethods::pca(snpMatrix, nPcs = 5) # SVD PCA

sumPCA <- as.data.frame(summary(pcaAM))

Scores <- as.data.frame(pcaMethods::scores(pcaAM))
Scores <- setDT(Scores, keep.rownames = TRUE)

Scores <- Scores %>%
  tidylog::mutate(founders = if_else(
    rn %in% c(
      "everest", "kanmark", "joe", "wb4458", "zenda",
      "bobdole", "overley", "monument", "jagger", "tam401",
      "tam107", "vona", "guymon", "prairie_red",
      "duke", "warrior", "deliver", "ok_rising", "jagalene"
    ), "yes", "no"
  ))

pc1 <- ggplot(data = Scores, aes(x = PC1,
              y = PC2)) +
  geom_point(position = "jitter", size = 3) +
  gghighlight::gghighlight(founders == "yes",
                           label_key = rn,
                           label_params = list(box.padding = 1,
                                               size = 7.5)
  ) +
  theme(aspect.ratio = 1:1) +
  labs(
    x = expression(paste("PC1 ","R"^{2}," = 5.2%")),
    y = expression(paste("PC2 ","R"^{2}," = 3.8%"))
  )
pc1

pc2 <- ggplot(data = Scores, aes(x = PC1,
                                 y = PC3)) +
  geom_point(position = "jitter", size = 3) +
  gghighlight::gghighlight(founders == "yes",
                           label_key = rn,
                           label_params = list(box.padding = 1,
                                               size = 7.5)
  ) +
  theme(aspect.ratio = 1:1) +
  labs(
    x = expression(paste("PC1 ","R"^{2}," = 5.2%")),
    y = expression(paste("PC3 ","R"^{2}," = 3.5%"))
  )
pc2

pc3 <- ggplot(data = Scores, aes(x = PC2,
                                 y = PC3)) +
  geom_point(position = "jitter") +
  gghighlight::gghighlight(founders == "yes",
                           label_key = rn,
                           label_params = list(box.padding = 1,
                                               size = 7.5)
  ) +
  theme(aspect.ratio = 1:1) +
  labs(
    x = expression(paste("PC2 ","R"^{2}," = 3.8%")),
    y = expression(paste("PC3 ","R"^{2}," = 3.5%"))
  )
pc3

ggarrange(pc1,pc2, ncol = 1)
ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/AMpanel/Figures/Figure4.png",
       height = 35,
       width = 20, units = "cm", dpi = 320
)

#### Genomic Prediction ####
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
  mutate(Cycle = row_number(),
         year = "17",
         Covariate = "No") %>% 
  select(Cycle, GRYLD, year, Covariate, everything()) %>%  
  pivot_longer(cols = RedEdge_20170609:GNDVI_20170331,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20170613"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"),
         date_min = as.Date("2017-03-31"),
         date_max = as.Date("2017-06-13"))

gs17_100_sum<- gs17_100 %>% 
  tidylog::group_by(trait_id,phenotype_date) %>% 
  summarise(Average = mean(Correlation), 
            stdev = sd(Correlation))

gs17_100<- gs17_100 %>% 
  inner_join(gs17_100_sum, by = c("trait_id","phenotype_date"))

gs17_covar_100<- gs17_covar_100 %>% 
  mutate(Cycle = row_number(),
         year = "17",
         Covariate = "Yes") %>% 
  select(Cycle, GRYLD_alone, year, Covariate, everything()) %>%  
  pivot_longer(cols = RedEdge_20170602:GNDVI_20170331,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20170613"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"),
         date_min = as.Date("2017-03-31"),
         date_max = as.Date("2017-06-13"))

gs17_covar_100_sum<- gs17_covar_100 %>% 
  tidylog::group_by(trait_id,phenotype_date) %>% 
  summarise(Average = mean(Correlation), 
            stdev = sd(Correlation)) 

gs17_covar_100<- gs17_covar_100 %>% 
  inner_join(gs17_covar_100_sum, by = c("trait_id","phenotype_date"))

gs18_100<- gs18_100 %>% 
  mutate(Cycle = row_number(),
         year = "18",
         Covariate = "No") %>% 
  select(Cycle, GRYLD, year, Covariate, everything()) %>%  
  pivot_longer(cols = RE_20180613:GNDVI_20171120,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  filter(phenotype_date > "20180101") %>% 
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20180615"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"),
         date_min = as.Date("2018-04-04"),
         date_max = as.Date("2018-06-15"))

gs18_100_sum<- gs18_100 %>% 
  tidylog::group_by(trait_id,phenotype_date) %>% 
  summarise(Average = mean(Correlation), 
            stdev = sd(Correlation)) 

gs18_100<- gs18_100 %>% 
  inner_join(gs18_100_sum, by = c("trait_id","phenotype_date"))

gs18_covar_100<- gs18_covar_100 %>% 
  mutate(Cycle = row_number(),
         year = "18",
         Covariate = "Yes") %>% 
  select(Cycle, GRYLD_alone, year, Covariate, everything()) %>%  
  pivot_longer(cols = RE_20180613:GNDVI_20171120,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  filter(phenotype_date > "20180101") %>% 
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20180615"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"),
         date_min = as.Date("2018-04-04"),
         date_max = as.Date("2018-06-15"))

gs18_covar_100_sum<- gs18_covar_100 %>% 
  tidylog::group_by(trait_id,phenotype_date) %>% 
  summarise(Average = mean(Correlation), 
            stdev = sd(Correlation)) 

gs18_covar_100<- gs18_covar_100 %>% 
  inner_join(gs18_covar_100_sum, by = c("trait_id","phenotype_date"))

gs19_100<- gs19_100 %>% 
  mutate(Cycle = row_number(),
         year = "19",
         Covariate = "No") %>% 
  select(Cycle, GRYLD, year, Covariate, everything()) %>%
  pivot_longer(cols = RE_20190624:GNDVI_20190103,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  filter(phenotype_date > "20190201") %>% 
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20190627"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"),
         date_min = as.Date("2019-04-12"),
         date_max = as.Date("2019-06-27"))

gs19_100_sum<- gs19_100 %>% 
  tidylog::group_by(trait_id,phenotype_date) %>% 
  summarise(Average = mean(Correlation), 
            stdev = sd(Correlation)) 

gs19_100<- gs19_100 %>% 
  inner_join(gs19_100_sum, by = c("trait_id","phenotype_date"))

gs19_covar_100<- gs19_covar_100 %>% 
  mutate(Cycle = row_number(),
         year = "19",
         Covariate = "Yes") %>% 
  select(Cycle, GRYLD_alone, year, Covariate, everything()) %>%  
  pivot_longer(cols = RE_20190624:GNDVI_20190103,
               names_to = "trait_id",
               values_to = "Correlation") %>% 
  separate(trait_id, into = c("trait_id","phenotype_date"),
           sep = "_") %>% 
  filter(phenotype_date > "20190201") %>%
  filter(trait_id != "GRVI") %>%
  filter(trait_id != "NIR") %>%
  filter(trait_id != "Nir") %>%
  filter(trait_id != "RE") %>%
  filter(trait_id != "RedEdge") %>% 
  mutate(phenotype_date = replace_na(phenotype_date, "20190627"),
         phenotype_date = as.Date(phenotype_date, format = "%Y%m%d"),
         date_min = as.Date("2019-04-12"),
         date_max = as.Date("2019-06-27"))

gs19_covar_100_sum<- gs19_covar_100 %>% 
  tidylog::group_by(trait_id,phenotype_date) %>% 
  summarise(Average = mean(Correlation), 
            stdev = sd(Correlation)) 

gs19_covar_100<- gs19_covar_100 %>% 
  inner_join(gs19_covar_100_sum, by = c("trait_id","phenotype_date"))

## Figures

gs_all_100<- bind_rows(gs17_100,gs18_100,gs19_100) 

write.table(gs_all_100, "./Phenotype_Database/gs_all.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
gs_all_100<- fread("./Phenotype_Database/gs_all.txt")

gs_covar_all_100<- bind_rows(gs17_covar_100, gs18_covar_100, gs19_covar_100)
write.table(gs_covar_all_100, "./Phenotype_Database/gs_covar_all.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
gs_covar_all_100<- fread("./Phenotype_Database/gs_covar_all.txt")

gs_all_100 %>% 
  ggplot(aes(x = phenotype_date, y = Correlation)) +
  geom_rect(aes(ymin = mean(GRYLD) - sd(GRYLD),
                ymax = mean(GRYLD) + sd(GRYLD),
                xmin = date_min,
                xmax = date_max),
            colour = "#99d8c9",
            fill = "#e5f5f9") +
  geom_jitter(size = 0.5, alpha = 0.25, colour = '#bdbdbd') +
  geom_point(aes(x = phenotype_date, y = Average),
             size = 2) +
  geom_errorbar(aes(x = phenotype_date, 
                    ymin = Average - stdev,
                    ymax = Average + stdev)) +
  geom_hline(aes(yintercept = mean(GRYLD)), linetype = 2) +
  facet_wrap(year~trait_id, scales = "free_x") +
  coord_cartesian(ylim = c(-0.5,1)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d") +
  labs(x = "Date")

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/AMpanel/Figures/Figure5.png",
       height = 30,
       width = 42.5, units = "cm", dpi = 320
)

gs_covar_all_100 %>% 
  ggplot(aes(x = phenotype_date, y = Correlation)) +
  geom_rect(aes(ymin = mean(GRYLD_alone) - sd(GRYLD_alone),
                ymax = mean(GRYLD_alone) + sd(GRYLD_alone),
                xmin = date_min,
                xmax = date_max),
            colour = "#99d8c9",
            fill = "#e5f5f9") +
  geom_jitter(size = 0.5, alpha = 0.25, colour = '#bdbdbd') +
  geom_point(aes(x = phenotype_date, y = Average),
             size = 2) +
  geom_errorbar(aes(x = phenotype_date, 
                    ymin = Average - stdev,
                    ymax = Average + stdev)) +
  geom_hline(aes(yintercept = mean(GRYLD_alone)), linetype = 2) +
  facet_wrap(year~trait_id, scales = "free_x") +
  coord_cartesian(ylim = c(-0.5,1)) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d") +
  labs(x = "Date")

ggsave("~/OneDrive - Kansas State University/Dissertation_Calvert/AMpanel/Figures/Figure6.png",
       height = 30,
       width = 42.5, units = "cm", dpi = 320
)


#### Association Mapping ####
######### Reading in files as a list of data frames

load.file <- function(filename) {
  d <- fread(file = filename, header = TRUE, check.names = F, data.table = F)
  d
}

## rrBLUP results

fileNames <- list.files(
  path = "./R/rrBlup/HypothesisEleven/PC4_2017",
  full.names = T,
  pattern = "_2017_4PC.txt$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all("_2017_4PC.txt")

rrB4_17 <- lapply(fileNames, load.file)

names(rrB4_17) <- traitNames

fileNames <- list.files(
  path = "./R/rrBlup/HypothesisEleven/PC4_2018",
  full.names = T,
  pattern = "_2018_4PC.txt$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all("_2018_4PC.txt")

rrB4_18 <- lapply(fileNames, load.file)

names(rrB4_18) <- traitNames

fileNames <- list.files(
  path = "./R/rrBlup/HypothesisEleven/PC4_2019",
  full.names = T,
  pattern = "_2019_4PC.txt$"
)
traitNames <- basename(fileNames) %>%
  str_remove_all("_2019_4PC.txt")

rrB4_19 <- lapply(fileNames, load.file)

names(rrB4_19) <- traitNames

snpPos <- rrB4_17$GNDVI_20170331[, 1:3]

rrB4_17F <- snpPos %>%
  bind_cols(map_dfr(rrB4_17, 4)) %>%
  mutate(
    Position = as.numeric(Position)
  ) %>%
  arrange(Chromosome, Position)

rrB4_18F <- snpPos %>%
  bind_cols(map_dfr(rrB4_18, 4)) %>%
  mutate(
    Position = as.numeric(Position)
  ) %>%
  arrange(Chromosome, Position)

rrB4_19F <- snpPos %>%
  bind_cols(map_dfr(rrB4_19, 4)) %>%
  mutate(
    Position = as.numeric(Position)
  ) %>%
  arrange(Chromosome, Position)

##### Comparison plots between GWAS methods and PCs ####

don_rrB4_17 <- rrB4_17F %>%
  # Compute chromosome size
  group_by(Chromosome) %>%
  summarise(chr_len = max(Position),
            .groups = "drop_last") %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_17F, ., by = c("Chromosome" = "Chromosome")) %>%
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate(BPcum = Position + tot)

axisdf <- don_rrB4_17 %>%
  group_by(Chromosome) %>%
  dplyr::summarize(
    center = (max(BPcum) + min(BPcum)) / 2,
    maximum = max(BPcum),
    .groups = "drop_last"
  )

don_rrB4_18 <- rrB4_18F %>%
  # Compute chromosome size
  group_by(Chromosome) %>%
  summarise(chr_len = max(Position),
            .groups = "drop_last") %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_18F, ., by = c("Chromosome" = "Chromosome")) %>%
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate(BPcum = Position + tot)

don_rrB4_19 <- rrB4_19F %>%
  # Compute chromosome size
  group_by(Chromosome) %>%
  summarise(chr_len = max(Position),
            .groups = "drop_last") %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  tidylog::select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(rrB4_19F, ., by = c("Chromosome" = "Chromosome")) %>%
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate(BPcum = Position + tot)

ggplot(don_rrB4_18, aes(x = BPcum, colour = as.factor(Chromosome))) +
  geom_point(aes(y = GRYLD),
             alpha = 0.5, size = 1
  ) +
  scale_color_manual(values = rep(c("#006837", "#41ae76", "#a1d99b"), 22)) +
  geom_hline(yintercept = -log10(0.05 / nrow(don_rrB4_18)), linetype = 2) +
  geom_vline(xintercept = 6979839046, linetype = 3, size = 1) +
  # custom X axis:
  scale_x_continuous(
    label = c(
      "1A", "1B", "1D",
      "2A", "2B", "2D",
      "3A", "3B", "3D",
      "4A", "4B", "4D",
      "5A", "5B", "5D",
      "6A", "6B", "6D",
      "7A", "7B", "7D"
    ),
    breaks = axisdf$center
  ) +
  scale_y_continuous(expand = c(0, 0.05)) + 
  theme(legend.position = "none") +
  labs(
    title = "GWAS results GRYLD 2018",
    subtitle = "Bonferroni Threshold alpha = 0.05",
    x = "Chromosome",
    y = "-log10(P)"
  ) +
  annotate(geom = "text", x = 6100000000, y = 6, label = "Rht-1B", size = 7.5)

ggsave("./AMPanel_Manuscript/Figures/GRYLD_association_2018.jpg",
       height = 30,
       width = 45, units = "cm", dpi = 350
)

blups18<- fread("./Phenotype_Database/ASREMLBlup_2018.txt")

gryld2018_snp<- don_rrB4_18 %>% 
  select(rs_number, chrom, pos, GRYLD) %>% 
  filter(GRYLD == max(GRYLD))
gryld2018_snp<- as.data.frame(snpMatrix) %>% 
  select(V8064)
gryld2018_snp<- setDT(as.data.frame(gryld2018_snp),keep.rownames = TRUE)
colnames(gryld2018_snp)

gryld2018<- blups18 %>% 
  select(Taxa, GRYLD) %>%
  left_join(gryld2018_snp, by = c("Taxa" = "rn"))

p2<- gryld2018 %>% 
  ggplot(aes(x = as.factor(V8064), y = GRYLD)) +
  geom_boxplot() + 
  labs(x = "SNP 8064",
       title = "GWAS GRYLD 2017/2018 Season")
p2

ggpubr::ggexport(p2, 
                 filename = "./AMPanel_Manuscript/Figures/snpEffects_gryld_2018.png",
                 width = 750,
                 height = 500)

