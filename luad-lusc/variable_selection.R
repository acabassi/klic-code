########################################
### LUAD vs LUSC: variable selection ###
########################################

rm(list=ls())
# setwd("/home/ac2051/rds/hpc-work/klic-code/luad-lusc/")
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/klic-code/luad-lusc/")

load("merged-data.RData")

library(protoclust)
library(pheatmap)

K <- 100

# dist_RNAseq <- dist(t(lu$RNA_seq), method = "euclidean")
# hc_RNAseq <- protoclust(dist_RNAseq)
# cut_RNAseq <- protocut(hc_RNAseq, k = K)
# 
# dist_mutation <- dist(t(lu$mutation), method = "binary")
# hc_mutation <- protoclust(dist_mutation)
# cut_mutation <- protocut(hc_mutation, k = K)
# 
# dist_miRNAseq <- dist(t(lu$miRNAseq), method = "euclidean")
# hc_miRNAseq <- protoclust(dist_miRNAseq)
# cut_miRNAseq <- protocut(hc_miRNAseq, k = K)

# save(cut_RNAseq, cut_mutation, cut_miRNAseq, file = "cut_protoclusts.RData")

load("cut_protoclusts.RData")

groups <- data.frame(Type = lu$label)
rownames(groups) <- rownames(lu$RNA_seq)

lu_reduced <- list()
lu_reduced$RNAseq <- lu$RNA_seq[,as.numeric(cut_RNAseq$protos)]
lu_reduced$mutation <- lu$mutation[,as.numeric(cut_mutation$protos)]
lu_reduced$miRNAseq <- lu$miRNAseq[,as.numeric(cut_miRNAseq$protos)]

save(lu_reduced, file = "reduced-datasets.RData")
save(groups, file = "annotations.RData")

png('reduced_RNAseq.png')
pheatmap(lu_reduced$RNAseq, annotation_row = groups,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

png('reduced_mutations.png')
pheatmap(lu_reduced$mutation, annotation_row = groups,
         show_colnames = FALSE, show_rownames = FALSE)
dev.off()

png('reduced_miRNAseq.png')
pheatmap(lu_reduced$miRNAseq, annotation_row = groups,
         show_rownames = FALSE, show_colnames = FALSE)
dev.off()

