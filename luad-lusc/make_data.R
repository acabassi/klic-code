rm(list=ls())

setwd("~/Dropbox/PhD Projects/MultipleKernels/Code/LUAD-LUSC/make-data/")

load("TCGA_LUSC.RData")
load("TCGA_LUAD.RData")

library(TCGA2STAT)

### Merge the LUAD OMICs-data into one R object
# methylation doesn't have any macth

luad1 <- OMICSBind(dat1 = luad_RNAseq$dat,   dat2 = luad_mutation$dat)
luad2 <- OMICSBind(dat1 = luad1$X, dat2 = luad_miRNAseq$dat)
luad3 <- OMICSBind(dat1 = luad1$Y, dat2 = luad2$Y)

colnames(luad2$X) == colnames(luad3$X)
colnames(luad2$X) == colnames(luad3$Y)

luad <- list(RNA_seq = luad2$X, mutation = luad3$X, miRNAseq = luad3$Y)

### Merge the LUSC OMICs-data into one R object
lusc1 <- OMICSBind(dat1 = lusc_RNAseq$dat,   dat2 = lusc_mutation$dat)
lusc2 <- OMICSBind(dat1 = lusc1$X, dat2 = lusc_miRNAseq$dat)
lusc3 <- OMICSBind(dat1 = lusc1$Y, dat2 = lusc2$Y)

colnames(lusc2$X) == colnames(lusc3$X)
colnames(lusc2$X) == colnames(lusc3$Y)

lusc <- list(RNA_seq = lusc2$X, mutation = lusc3$X, miRNAseq = lusc3$Y)

save(luad, file = "TCGA_LUAD_match.RData")
save(lusc, file = "TCGA_LUSC_match.RData")
