##########################
### LUAD vs LUSC: KLIC ###
##########################

rm(list=ls())
# setwd("/home/ac2051/rds/hpc-work/klic-code/luad-lusc/")
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/klic-code/luad-lusc/")

# library(devtools)
# install_github("acabassi/coca")
# currently using version 1.0.0-alpha: DOI 10.5281/zenodo.2562342
library(coca)
# install_github("acabassi/klic") 
# currently using version 1.0.0-alpha: DOI 10.5281/zenodo.2562339
library(klic)

library(mclust)

load("reduced-datasets.RData")
load("annotations.RData")

klic <- klic(data = lu_reduced, M = 3, individualMaxK = 10, individualClAlgorithm = "hclust",
             globalMaxK = 10, scale = TRUE, annotations = groups, 
             ccClMethods = "hc", ccDistHCs = c("euclidean", "binary", "euclidean"))
