library(MultiPhen)
library(readr)
library(tidyverse)
library(readr)
library(data.table)
library(ggplot2)
library(qqman)

rm(list = objects()); ls()

all.opts<- mPhen.options("all", descr = T)

misc.optsDesc<- mPhen.options("misc", descr = T)
pheno.optsDesc = mPhen.options("pheno.input", descr = T)

pheno.opts = mPhen.options("pheno.input")
misc.opts<- mPhen.options("misc")
options("mPhen.log10p"=FALSE)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

print(names(pheno))
print(pheno$limit)

opts.desc<- mPhen.options("regression",descr=TRUE)
phenoDat<- pheno$pheno

opts<- mPhen.options("regression")
phenoDat<- pheno$pheno
phenoOb<- mPhen.preparePheno(pheno, opts = opts)
numPhenos<- length(phenoOb$phenN)

geno.optsDesc<- mPhen.options("geno.input", descr = T)
geno.opts<- mPhen.options("geno.input")
geno.opts$mPhen.format = "GT"
geno.opts$mPhen.batch = 1e+10
geno.opts$mPhen.numGenoPCs = 5

geno<- mPhen.readGenotypes(
  "~/Dropbox/Research_Poland_Lab/AM Panel/plink/amPanel/AMsnpChipImputed.gt.vcf.gz",
  opts = geno.opts,
  indiv = rownames(pheno$pheno))
genoOb<- geno$genoData
dim(genoOb)

print(opts)

resDir = "~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate"
towrite = list(long.txt = FALSE,   wide.txt = TRUE)
toplot = list(
  .manh = TRUE,
  .qq = TRUE,
  .heatm = TRUE,
  .fprint = FALSE
)

plotopts = mPhen.options("plot")
plotopts.optsDesc = mPhen.options("plot", descr = TRUE)
plotopts$mPhen.noPhenosPerPlot = 1
plotopts$mPhen.Colv = F
plotopts$mPhen.onlyShowJoint = TRUE


###### 2018 Models

opts$mPhen.variable.selection=FALSE
opts$mPhen.JointModel=TRUE
opts$mPhen.inverseRegress=TRUE
opts$mPhen.adjustSinglePv=FALSE
opts$mPhen.calcHWE = FALSE

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs Back Selection
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,3,4,8,10,17,19,20,22,23,24,
                                                  27,29,30,31,32,34,36,41,45,51,
                                                  54,55,57,58,59,60,62,63,64,71,
                                                  72,77,78,80,83,86,87,89,90,93,
                                                  94,95,97,98)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluBckSelec5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluBckSelec5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluBckSelec5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)
str(singTest)
trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluBckSelec5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluBckSelec5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)


pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs Back Selection with control for Multicollinearity
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,3,4,8,10,19,23,24,27,32,34,
                                                  41,45,51,62,77,78,80,83,86,87,
                                                  89,90,93,94,95,97,98)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluBckSelecVIF5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluBckSelecVIF5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluBckSelecVIF5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)
str(singTest)
trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluBckSelecVIF5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluBckSelecVIF5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)


pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs Back Selection with control for Multicollinearity and significance
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,3,4,8,19,41,51,62,77,89)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluBckSelecVIFSig5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluBckSelecVIFSig5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluBckSelecVIFSig5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluBckSelecVIFSig5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluBckSelecVIFSig5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs GRWT correlated network
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,19,34,60,61,75,76)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluCorr5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluCorr5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluCorr5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]


qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluCorr5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluCorr5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs GRWT correlated*heritability network
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,8,4,5,92,93,46,84,94,47,75,48,
                                                  49,22,36,63,77,87,23,37,51,64,78,88,98)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluCorrHer5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluCorrHer5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluCorrHer5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluCorrHer5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluCorrHer5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs LASSO variable Selection
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,52,55,27,71,19,57,63,26,59,24,20,41,17,62,54)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluLasSelec5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluLasSelec5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluLasSelec5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluLasSelec5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluLasSelec5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs LASSO Selected Variables CV
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,3,8,19,20,38,62,76,77,89)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluLasSelecCV5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluLasSelecCV5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluLasSelecCV5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluLasSelecCV5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluLasSelecCV5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

## BLUEs LASSO Selected Minimum Lambda
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(1,3,4,5,8,12,19,32,35,38,40,45,46,47,48,50,51,62,64,65,76,77,78,83,85,88,89)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2018bluLasSelecML5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2018bluLasSelecML5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2018bluLasSelecML5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2018bluLasSelecML5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2018 BluLasSelecML5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno18_NaNa.txt",
  opts = pheno.opts)

##### 2017


pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)
  
geno.optsDesc<- mPhen.options("geno.input", descr = T)
geno.opts<- mPhen.options("geno.input")
geno.opts$mPhen.format = "GT"
geno.opts$mPhen.batch = 1e+10
geno.opts$mPhen.numGenoPCs = 5

geno<- mPhen.readGenotypes(
  "~/Dropbox/Research_Poland_Lab/AM Panel/plink/amPanel/AMsnpChipImputed.gt.vcf.gz",
  opts = geno.opts,
  indiv = rownames(pheno$pheno))
genoOb<- geno$genoData
dim(genoOb)

## Blu Back Selection
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,1,3,7,9,10,11,12,14,19,21,27,28,32,36,37,39,40,45,46,47,48,51,53,54,55,58,59,60,61,62,63,64,65,66,72)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017BluBckSelec5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017BluBckSelec5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017BluBckSelec5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017BluBckSelec5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluBckSelec5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## Blu Back Selection with MultiCollinearity
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,1,3,7,9,10,14,21,32,36,39,45,46,48,51,55,58,59,62,63,66,72)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017BluBckSelecVIF5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017BluBckSelecVIF5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017BluBckSelecVIF5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017BluBckSelecVIF5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluBckSelecVIF5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## BLUE Back Selection with MultiCollinearity
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,3,7,9,10,21,32,45,55,59,63,66,72)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017BluBckSelecVIFSig5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017BluBckSelecVIFSig5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017BluBckSelecVIFSig5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017BluBckSelecVIFSig5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluBckSelecVIFSig5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## BLUEs GRWT correlated network
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,10,11,12,13,18,19,20,21,26,27,28,29,35,52,53)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017bluCorr5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017bluCorr5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017bluCorr5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]


qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017bluCorr5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluCorr5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## BLUEs GRWT correlated*heritability network
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,29,20,28,21,19,27,13,12,11,26,18,10,35,36,34)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017bluCorrHer5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017bluCorrHer5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017bluCorrHer5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017bluCorrHer5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluCorrHer5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## BLUEs LASSO variable Selection
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,12,27,58,10,11,28,36,14,59)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017bluLasSelec5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017bluLasSelec5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017bluLasSelec5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017bluLasSelec5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluLasSelec5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## BLUEs LASSO Selected Variables CV
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,3,7,16,22,26,27,29,30,65,73)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017bluLasSelecCV5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017bluLasSelecCV5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017bluLasSelecCV5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017bluLasSelecCV5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluLasSelecCV5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)

pheno<- mPhen.readPhenoFiles(
  "~/Dropbox/Research_Poland_Lab/AM Panel/Phenotype_Database/cleanPheno17_NaNa.txt",
  opts = pheno.opts)

## BLUEs LASSO Selected Minimum Lambda
print(names(pheno))
print(pheno$limit)
pheno$limit$phenotypes = pheno$limit$phenotypes[c(2,1,3,4,7,9,10,16,21,22,24,27,29,30,33,43,44,47,54,57,59,60,63,65,69,71,72,73)]
print(pheno$limit)

phenoOb<- mPhen.preparePheno(pheno, opts = opts)

results2017bluLasSelecML5 = mPhen.assoc(genoOb, phenoOb,  opts)

output = mPhen.writeOutput(
  results2017bluLasSelecML5,
  output = resDir,
  geno = genoOb,
  towrite = towrite,
  toplot = toplot,
  opts = plotopts
)

dev.off()

singTest<- read.table("~/Dropbox/Research_Poland_Lab/AM Panel/R/MultiPhen/Multivariate/results2017bluLasSelecML5.wide.txt",
                      header = T, sep = "")
singTest<- separate(singTest, "rsid", c("CHR","BP"), sep = "_")
singTest$CHR <- as.numeric(singTest$CHR)
singTest$BP <- as.numeric(singTest$BP)

trials<- singTest[which(singTest$label == "pval"), c(5,3:4,6:ncol(singTest))]

qqman.plot <- function(x, ...) {
  md <- names(x) %in% c("SNP", "BP", "CHR")
  traits <- names(x[!md])
  info<- x[,md]
  
  for (i in traits) {
    P <- x[[i]]
    dat<- cbind(info, P)
    mypath <- file.path("~","Dropbox","Research_Poland_Lab","AM Panel","R",
                        "MultiPhen","Figures",paste("2017bluLasSelecML5", i, ".pdf", sep = ""))
    pdf(file = mypath, width = 20, height = 10)
    qqman::manhattan(dat,
                     main = paste("2017 BluLasSelecML5 ",i),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0,20),
                     col = c("blue","grey40","black"),
                     chrlabs = c("1","1B","1D",
                                 "2A","2B","2D",
                                 "3A","3B","3D",
                                 "4A","4B","4D",
                                 "5A","5B","5D",
                                 "6A","6B","6D",
                                 "7A","7B","7D"),
                     genomewideline = -log10(0.05 / nrow(dat)),
                     #suggestiveline = -log10(1 - (1 - 0.05)^(1/nrow(dat))),
                     logp = T)
    dev.off()
    
  }
}

qqman.plot(trials)


