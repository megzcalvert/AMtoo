rm(list = objects())
ls()

library(readr)
library(data.table)
library(tidyverse)
library(RMySQL)
library(janitor)
library(GGally)
require(lubridate)
library(car)
library(tidylog)
library(broom)
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
# useful infos **reproducible research**
sessionInfo()

#### Load data from database ####
## Connect to database
wheatgenetics <- dbConnect(MySQL(),
  user = rstudioapi::askForPassword("Database user"),
  dbname = "wheatgenetics",
  host = "beocat.cis.ksu.edu",
  password =
    rstudioapi::askForPassword("Database password"),
  port = 6306
)

# SQL Query to get all AM Panel phenotype data

pheno_query <- "SELECT phenotype.entity_id,
phenotype.trait_id,
phenotype.phenotype_value,
phenotype.phenotype_date,
phenotype.phenotype_person,
plot.plot_name AS 'Variety',
plot.block,
plot.rep,
plot.range,
plot.column
FROM wheatgenetics.phenotype LEFT JOIN wheatgenetics.plot 
ON plot.plot_id = phenotype.entity_id
WHERE wheatgenetics.phenotype.entity_id LIKE '18ASH30%' OR
wheatgenetics.phenotype.entity_id LIKE '17ASH1%' OR
wheatgenetics.phenotype.entity_id LIKE '19RKY0%';"

# run the query to get plot information
pheno <- dbGetQuery(wheatgenetics, pheno_query)

# save original data
saveRDS(pheno, "./Phenotype_Database/Pheno171819.RDS")

# disconnect from database
dbDisconnect(wheatgenetics)

#### Data summaries of phenotypic data ####

## clean workspace
rm(list = objects())
ls()

getwd()

# Read in original data
pheno_long <- readRDS("./Phenotype_Database/Pheno171819.RDS")

glimpse(pheno_long)

# What lines are in the data and how often
iniNames <- tabyl(pheno_long$Variety)
names(iniNames)[1:2] <- c("Variety", "Count")

## Long Format Summaries with modified names
pheno_long <- pheno_long %>%
  mutate(
    Variety = tolower(Variety),
    Variety = str_replace_all(Variety, " ", "_"),
    Variety = str_replace_all(Variety, "-", "_"),
    Variety = str_replace_all(Variety, "'", "")
  ) %>%
  glimpse()

# Function to add a column based on a portion of text in another column
ff <- function(x, patterns, replacements = patterns, fill = NA, ...) {
  stopifnot(length(patterns) == length(replacements))

  ans <- rep_len(as.character(fill), length(x))
  empty <- seq_along(x)

  for (i in seq_along(patterns)) {
    greps <- grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] <- replacements[[i]]
    empty <- empty[!greps]
  }

  return(ans)
}

# Adding a year column based on the entity_id
pheno_long$year <- ff(pheno_long$entity_id,
  c("17ASH", "18ASH", "19RKY"),
  c("17", "18", "19"),
  "NA",
  ignore.case = TRUE
)

# Summary of lines with modified names
modNames <- tabyl(pheno_long$Variety)
names(modNames)[1:2] <- c("Variety", "Count")

# Summary of variables
dataMaid::summarize(pheno_long[, c("Variety", "trait_id", "rep", "year")])
knitr::kable(tabyl(pheno_long$trait_id),
  caption = "Number of observations per trait"
)
knitr::kable(tabyl(pheno_long$rep),
  caption = "Number of observations of reps"
)
knitr::kable(tabyl(pheno_long$year),
  caption = "Number of observations per year"
)

### Removing lines that do not fit the assumptions
pheno_long <- pheno_long %>%
  filter(Variety != "blank") %>%
  filter(!is.na(Variety)) %>%
  filter(rep != 0) %>%
  filter(!is.na(rep)) %>%
  filter(trait_id != "NOTES")

finalNames <- tabyl(pheno_long$Variety)

glimpse(pheno_long)
tabyl(pheno_long$trait_id)
tabyl(pheno_long$rep)
tabyl(pheno_long$year)

#### Conversion to Wide format and other data cleaning ####

# Removing confounding phenotypes  - Jesse's Code
######################################################

pheno_awns <- pheno_long %>%
  filter(trait_id == "AWNS")

## set the date as a Date type
## remove some rows that are missing the phenotype date
pheno_pcthead <- pheno_long %>%
  filter(trait_id == "PCTHEAD") %>%
  mutate(phenotype_date = as.Date(phenotype_date)) %>%
  filter(!is.na(phenotype_date))

## set the date as a Date type
pheno_long <- pheno_long %>%
  mutate(phenotype_date = as.Date(phenotype_date)) %>%
  filter(trait_id != "AWNS") %>%
  filter(trait_id != "PCTHEAD")

## sort the rows by plot and date
pheno_pcthead <- pheno_pcthead[order(
  pheno_pcthead$entity_id,
  pheno_pcthead$phenotype_date
), ]

runLogReg <- function(dates, pctHEAD) {
  res <- data.frame(phi1 = NA, phi2 = NA, phi3 = NA, hday = NA, hdate = NA)

  if (sum(!is.na(pctHEAD)) == 0) {
    return(list(res = res, dates = dates, pctHEAD = pctHEAD, pred = NA))
  } ## skip if no phenotype data
  if (max(pctHEAD, na.rm = TRUE) < 60) {
    return(list(res = res, dates = dates, pctHEAD = pctHEAD, pred = NA))
  } ## skip plots that never reach 50% heading

  days <- yday(dates) ## convert to day of year

  ## add some dummy variables for 0% and 100% at 10, 20, 30 days before/after the phenotyping range ##
  days <- c(days, min(days) - c(10, 20, 30), max(days) + c(10, 20, 30))
  pctH <- c(pctHEAD, c(0, 0, 0, 100, 100, 100))

  ## find initial starting values for phi2 and phi3,  fix phi1 at 100
  phi1 <- 100
  phi2 <- coef(lm(logit(pctH / phi1) ~ days))[1]
  phi3 <- coef(lm(logit(pctH / phi1) ~ days))[2]

  heading_model <- try(nls(pctH ~ 100 / (1 + exp(-(phi2 + phi3 * days))), start = list(phi2 = phi2, phi3 = phi3), trace = FALSE), silent = TRUE)

  if (class(heading_model) == "try-error") {
    return(list(res = res, dates = dates, pctHEAD = pctHEAD, pred = NA))
  } ## skip if doesn't converge

  # update model coefficients for predictions
  phi2 <- coef(heading_model)[1]
  phi3 <- coef(heading_model)[2]

  ## get predicted value
  pred <- getPred(phi1, phi2, phi3)

  hday <- pred$x[which.min(abs(pred$y - 50))] ## get heading day of year
  hdate <- as.Date(paste(year(dates[1]), "-01-01", sep = "")) + hday ## convert to date

  res <- data.frame(phi1, phi2, phi3, hday, hdate) ## set values

  return(list(res = res, dates = dates, pctHEAD = pctHEAD, pred = pred))
}

getPred <- function(phi1, phi2, phi3) {
  x <- c((80 * 10):(150 * 10)) / 10 # construct a range of x values for all days of year
  y <- phi1 / (1 + exp(-(phi2 + phi3 * x))) # predicted y values
  pred <- data.frame(x, y) # create the prediction data frame

  return(pred)
}

#################### TEST NON-LINEAR FIT #############################

## heading date calculation results
## get list of plots
plots <- unique(pheno_pcthead$entity_id)
## work with only 17ASH plots
plots17 <- plots[grepl("17ASH", plots)]
## data frame to store results
head.dates <- data.frame(
  plot_id = plots17,
  day = NA,
  phi1 = NA,
  phi2 = NA,
  phi3 = NA
)

## first example
test.plot <- "17ASH10002"

sample <- pheno_pcthead[pheno_pcthead$entity_id == test.plot, ]
sample$phenotype_value <- as.numeric(sample$phenotype_value)

test <- runLogReg(
  dates = sample$phenotype_date,
  pctHEAD = sample$phenotype_value
)
hddt17 <- as.data.frame(test$res)

## data frame to store results

for (i in plots17) {
  ## subsample the data for the plot for analysis
  print(paste("Working on:", i))
  ## run for visual scores
  sample <- pheno_pcthead[pheno_pcthead$entity_id == i, ]
  sample$phenotype_value <- as.numeric(sample$phenotype_value)
  vis.logistic <- runLogReg(
    dates = sample$phenotype_date,
    pctHEAD = sample$phenotype_value
  )
  n <- as.data.frame(vis.logistic$res)
  hddt17 <- rbind(hddt17, n)
}

hddt17 <- hddt17[2:nrow(hddt17), ]
hddt17 <- cbind(plots17, hddt17)

plotDates17 <- hddt17[, c(1, 5:6)]
names(plotDates17) <- c("plots17", "numHddt17", "hddt17")

## work with only 18ASH plots
plots18 <- plots[grepl("18ASH", plots)]
## data frame to store results
head.dates <- data.frame(
  plot_id = plots18,
  day = NA,
  phi1 = NA,
  phi2 = NA,
  phi3 = NA
)

## first example
test.plot <- "18ASH30001"

sample <- pheno_pcthead[pheno_pcthead$entity_id == test.plot, ]
sample$phenotype_value <- as.numeric(sample$phenotype_value)

test <- runLogReg(
  dates = sample$phenotype_date,
  pctHEAD = sample$phenotype_value
)
hddt18 <- as.data.frame(test$res)

for (i in plots18) {
  ## subsample the data for the plot for analysis
  print(paste("Working on:", i))
  ## run for visual scores
  sample <- pheno_pcthead[pheno_pcthead$entity_id == i, ]
  sample$phenotype_value <- as.numeric(sample$phenotype_value)
  vis.logistic <- runLogReg(
    dates = sample$phenotype_date,
    pctHEAD = sample$phenotype_value
  )
  n <- as.data.frame(vis.logistic$res)
  hddt18 <- rbind(hddt18, n)
}

hddt18 <- hddt18[2:nrow(hddt18), ]
hddt18 <- cbind(plots18, hddt18)

plotDates18 <- hddt18[, c(1, 5:6)]
names(plotDates18) <- c("plots18", "numHddt18", "hddt18")

## work with only 19RKY plots
## get list of plots
plots <- unique(pheno_pcthead$entity_id)
## work with only 17ASH plots
plots19 <- plots[grepl("19RKY", plots)]
## data frame to store results
head.dates <- data.frame(
  plot_id = plots19,
  day = NA,
  phi1 = NA,
  phi2 = NA,
  phi3 = NA
)

## first example
test.plot <- "19RKY00728"

sample <- pheno_pcthead[pheno_pcthead$entity_id == test.plot, ]
sample$phenotype_value <- as.numeric(sample$phenotype_value)

test <- runLogReg(
  dates = sample$phenotype_date,
  pctHEAD = sample$phenotype_value
)
hddt19 <- as.data.frame(test$res)

## data frame to store results

for (i in plots19) {
  ## subsample the data for the plot for analysis
  print(paste("Working on:", i))
  ## run for visual scores
  sample <- pheno_pcthead[pheno_pcthead$entity_id == i, ]
  sample$phenotype_value <- as.numeric(sample$phenotype_value)
  vis.logistic <- runLogReg(
    dates = sample$phenotype_date,
    pctHEAD = sample$phenotype_value
  )
  n <- as.data.frame(vis.logistic$res)
  hddt19 <- rbind(hddt19, n)
}

hddt19 <- hddt19[2:nrow(hddt19), ]
hddt19 <- cbind(plots19, hddt19)

plotDates19 <- hddt19[, c(1, 5:6)]
names(plotDates19) <- c("plots19", "numHddt19", "hddt19")

#####################################################################

pheno_long$phenotype_value <- as.numeric(
  as.character(pheno_long$phenotype_value)
)

write.table(plotDates17, "./Phenotype_Database/HDDT2017.txt",
  col.names = T, row.names = F, sep = "\t", quote = F
)
write.table(plotDates18, "./Phenotype_Database/HDDT2018.txt",
  col.names = T, row.names = F, sep = "\t", quote = F
)
write.table(plotDates19, "./Phenotype_Database/HDDT2019.txt",
            col.names = T, row.names = F, sep = "\t", quote = F
)
write.table(pheno_long,
  file = "./Phenotype_Database/Pheno_Long171819.txt",
  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)
write.table(pheno_awns,
  file = "./Phenotype_Database/Pheno_LongAwns.txt",
  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

rm(
  hddt17, hddt18, head.dates, iniNames, modNames, n, pheno_long, pheno_pcthead,
  sample, test, vis.logistic, i, plots, plots17, plots18, plots19, getPred, 
  runLogReg, test.plot
)

pheno_long <- fread("./Phenotype_Database/Pheno_Long171819.txt", header = T)
plotDates17 <- fread("./Phenotype_Database/HDDT2017.txt")
plotDates18 <- fread("./Phenotype_Database/HDDT2018.txt")
plotDates19 <- fread("./Phenotype_Database/HDDT2019.txt")

# Selecting pheno of choice

unique(pheno_long$trait_id)

pheno17 <- pheno_long %>%
  filter(year == "17") %>%
  filter(trait_id == "GRYLD") %>%
  dplyr::arrange(entity_id) %>%
  left_join(plotDates17,
    by = c("entity_id" = "plots17")
  ) %>%
  tidylog::select(entity_id, phenotype_value, hddt17, numHddt17) %>%
  mutate(hddt17 = as.Date(hddt17))

pheno18 <- pheno_long %>%
  filter(year == "18") %>%
  filter(trait_id == "GRYLD") %>%
  dplyr::arrange(entity_id) %>%
  left_join(plotDates18,
    by = c("entity_id" = "plots18")
  ) %>%
  tidylog::select(entity_id, phenotype_value, hddt18, numHddt18) %>%
  mutate(hddt18 = as.Date(hddt18))

pheno19 <- pheno_long %>%
  filter(year == "19") %>%
  filter(trait_id == "GRYLD") %>%
  dplyr::arrange(entity_id) %>%
  left_join(plotDates19,
            by = c("entity_id" = "plots19")
  ) %>%
  tidylog::select(entity_id, phenotype_value, hddt19, numHddt19) %>%
  mutate(hddt19 = as.Date(hddt19))

str(pheno17)

# LOESS - Locally Weighted Scatterplot Smoothing
# popular for regression through time lines

ggplot(pheno17, aes(x = hddt17, y = phenotype_value)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.25) +
  # geom_smooth(method = "loess",color = "#7fbf7b",alpha = 0.25) +
  scale_x_date(date_breaks = "3 day", date_labels = "%b %d") +
  labs(
    y = "GRYLD",
    x = "HDDT",
    title = "GRYLD vs HDDT 2016/2017"
  ) +
  annotate("text",
    x = as.Date("2017-05-13"), y = 6.5,
    label = "italic(r) == -0.468",
    parse = TRUE, size = 8
  )

ggplot(pheno18, aes(x = hddt18, y = phenotype_value)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.25) +
  # geom_smooth(method = "loess",color = "#7fbf7b",alpha = 0.25) +
  scale_x_date(date_breaks = "3 day", date_labels = "%b %d") +
  labs(
    y = "GRYLD",
    x = "HDDT",
    title = "GRYLD vs HDDT 2017/2018"
  ) +
  annotate("text",
    x = as.Date("2018-05-19"), y = 5,
    label = "italic(r) == 0.0008",
    parse = TRUE, size = 8
  )

ggplot(pheno19, aes(x = hddt19, y = phenotype_value)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.25) +
  # geom_smooth(method = "loess",color = "#7fbf7b",alpha = 0.25) +
  scale_x_date(date_breaks = "3 day", date_labels = "%b %d") +
  labs(
    y = "GRYLD",
    x = "HDDT",
    title = "GRYLD vs HDDT 2018/2019"
  ) +
  annotate("text",
           x = as.Date("2019-05-23"), y = 9.5,
           label = "italic(r) == -0.510",
           parse = TRUE, size = 8
  )

tidy(cor.test(pheno17$numHddt17, pheno17$phenotype_value))
tidy(cor.test(pheno18$numHddt18, pheno18$phenotype_value))
tidy(cor.test(pheno19$numHddt19, pheno19$phenotype_value))

linReg2017 <- lm(phenotype_value ~ numHddt17, data = pheno17)
summary(linReg2017)
(linReg2017tidy <- tidy(linReg2017))
(linReg2017augmented <- augment(linReg2017))
(linReg2017glance <- glance(linReg2017))

outlierTest(linReg2017)
plot(linReg2017)

linReg2018 <- lm(phenotype_value ~ numHddt18, data = pheno18)
summary(linReg2018)
(linReg2018tidy <- tidy(linReg2018))
linReg2018augmented <- augment(linReg2018)
linReg2018glance <- glance(linReg2018)

outlierTest(linReg2018)
plot(linReg2018)

linReg2019 <- lm(phenotype_value ~ numHddt19, data = pheno19)
summary(linReg2019)
(linReg2019tidy <- tidy(linReg2019))
linReg2019augmented <- augment(linReg2019)
linReg2019glance <- glance(linReg2019)

outlierTest(linReg2019)
plot(linReg2019)

# Influential Observations
# added variable plots
avPlot(linReg2017, variable = "numHddt17")
# Cook's D plot
# identify D values > 4/(n-k-1)
cutoff <- 4 / ((nrow(pheno17) - length(linReg2017$coefficients) - 2))
plot(linReg2017, which = 4, cook.levels = cutoff)
# Influence Plot
influencePlot(linReg2017,
  id.method = "identify",
  main = "Influence Plot",
  sub = "Circle size is proportial to Cook's Distance"
)

avPlot(linReg2018, variable = "numHddt18")
# Cook's D plot
# identify D values > 4/(n-k-1)
cutoff <- 4 / ((nrow(pheno18) - length(linReg2018$coefficients) - 2))
plot(linReg2018, which = 4, cook.levels = cutoff)
# Influence Plot
influencePlot(linReg2018,
  id.method = "identify",
  main = "Influence Plot",
  sub = "Circle size is proportial to Cook's Distance"
)

avPlot(linReg2019, variable = "numHddt19")
# Cook's D plot
# identify D values > 4/(n-k-1)
cutoff <- 4 / ((nrow(pheno19) - length(linReg2019$coefficients) - 2))
plot(linReg2019, which = 4, cook.levels = cutoff)
# Influence Plot
influencePlot(linReg2019,
              id.method = "identify",
              main = "Influence Plot",
              sub = "Circle size is proportial to Cook's Distance"
)

# Non-normality

sresid <- studres(linReg2017)
hist(sresid,
  freq = FALSE,
  main = "Distribution of Studentized Residuals"
)
xfit <- seq(min(sresid), max(sresid), length = 40)
yfit <- dnorm(xfit)
lines(xfit, yfit)

sresid <- studres(linReg2018)
hist(sresid,
  freq = FALSE,
  main = "Distribution of Studentized Residuals"
)
xfit <- seq(min(sresid), max(sresid), length = 40)
yfit <- dnorm(xfit)
lines(xfit, yfit)

sresid <- studres(linReg2019)
hist(sresid,
     freq = FALSE,
     main = "Distribution of Studentized Residuals"
)
xfit <- seq(min(sresid), max(sresid), length = 40)
yfit <- dnorm(xfit)
lines(xfit, yfit)

# Evaluate homoscedasticity
# non-constant error variance test
ncvTest(linReg2017)
# plot studentized residuals vs. fitted values
spreadLevelPlot(linReg2017)

ncvTest(linReg2018)
# plot studentized residuals vs. fitted values
spreadLevelPlot(linReg2018)

ncvTest(linReg2019)
# plot studentized residuals vs. fitted values
spreadLevelPlot(linReg2019)

# Evaluate Nonlinearity
# component + residual plot
crPlots(linReg2017)
crPlots(linReg2018)
crPlots(linReg2019)

# Test for Autocorrelated Errors
durbinWatsonTest(linReg2017)
durbinWatsonTest(linReg2018)
durbinWatsonTest(linReg2019)


gvmodel <- gvlma(linReg2017)
summary(gvmodel)

gvmodel <- gvlma(linReg2018)
summary(gvmodel)

gvmodel <- gvlma(linReg2019)
summary(gvmodel)
