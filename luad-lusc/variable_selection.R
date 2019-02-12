########################################
### LUAD vs LUSC: variable selection ###
########################################

rm(list=ls())
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/klic-code/luad-lusc/")

load("merged-data.RData")

library(protoclust)

K <- 100

dist_RNAseq <- dist(t(lu$RNA_seq), method = "euclidean")
hc_RNAseq <- protoclust(dist_RNAseq)
cut_RNAseq <- protocut(hc_RNAseq, k = K)

dist_mutation <- dist(lu$mutation, method = "binary")
hc_mutation <- protoclust(dist_mutation)
cut_mutation <- protocut(hc_mutation, k = K)

dist_miRNAseq <- dist(lu$miRNAseq, method = "euclidean")
hc_miRNAseq <- protoclust(dist_miRNAseq)
cut_miRNAseq <- protocut(hc_miRNAseq, k = K)

save(cut_RNAseq, cut_mutation, cut_miRNAseq, file = "cut_protoclusts.R")