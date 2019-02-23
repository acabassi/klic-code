##########################
### LUAD vs LUSC: COCA ###
##########################

rm(list=ls())
setwd("/home/ac2051/rds/hpc-work/klic-code/luad-lusc/")
# setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/klic-code/luad-lusc/")

library(devtools)
install_github("acabassi/coca")
# currently using version 1.0.0-alpha: DOI 10.5281/zenodo.2562342
library(coca)
library(mclust)
library(pheatmap)
library(ggbiplot)

load("merged-data.RData")
load("annotations.RData")

### PCA 

pca_RNAseq <- prcomp(lu$RNA_seq, center = TRUE, scale = TRUE)
ggbiplot(pca_RNAseq, var.axes = F, groups = groups$Type)
ggsave("pca_RNAseq.png")

# Here I try to do PCA on the mutation but IT IS BINARY DATA
pca_mutation <- prcomp(lu$mutation, center = FALSE, scale = FALSE)
ggbiplot(pca_mutation, var.axes = F, groups = groups$Type)
ggsave("pca_mutation.png")

pca_miRNAseq <- prcomp(lu$miRNAseq, center = TRUE, scale = TRUE)
ggbiplot(pca_miRNAseq, var.axes = F, groups = groups$Type)
ggsave("pca_miRNAseq.png")

### Scale non-binary datasets

lu$RNAseq <- scale(lu$RNA_seq, scale = TRUE, center = TRUE)
lu$RNA_seq <- NULL
lu$miRNAseq <- scale(lu$miRNAseq, scale = TRUE, center = TRUE)

### Choice of the number of clusters

dist_RNAseq <- dist(lu$RNAseq, method = "euclidean")
hc_RNAseq <- hclust(dist_RNAseq, method = "average")

groups_RNAseq <- groups
for(i in 2:3){
    groups_RNAseq[[i]] <- as.factor(cutree(hc_RNAseq, k = i))
}

table(groups_RNAseq[[2]])
# Not very interesting

# Let's try out different settings

hc_RNAseq <- hclust(dist_RNAseq, method = "ward.D")

for(i in 2:10){
    groups_RNAseq[[i]] <- as.factor(cutree(hc_RNAseq, k = i))
}

table(groups_RNAseq[[2]])
# Maybe the ward.D linkage gives more sensible results

png('RNAseq.png')
pheatmap(lu$RNAseq[,1:500], annotation_row = groups_RNAseq,
         show_rownames = FALSE, show_colnames = FALSE, annotation_legend = TRUE, 
         cluster_rows = FALSE)
dev.off()

adjustedRandIndex(groups_RNAseq[[3]], groups_RNAseq$Type) # 0.1261197

# Mutations

dist_mutation <- dist(lu$mutation, method = "binary")
# dist_mutation <- proxy::dist(lu$mutation, 
                             # by_rows = TRUE, method = "Jaccard")
hc_mutation <- hclust(dist_mutation, method = "median")

groups_mutation <- groups
for(i in 2:20){
    groups_mutation[[i]] <- as.factor(cutree(hc_mutation, k = i))
}

png('mutations.png')
pheatmap(lu$mutation[,1:500], annotation_row = groups_mutation,
         show_rownames = FALSE, show_colnames = FALSE, annotation_legend = FALSE, 
         cluster_rows = FALSE, clustering_distance_rows = "binary")
dev.off()

adjustedRandIndex(groups_mutation[[3]], groups_mutation$Type) # -0.003949293

### miRNAseq

dist_miRNAseq <- dist(lu$miRNAseq, method = "manhattan")
hc_miRNAseq <- hclust(dist_miRNAseq, method = "ward.D")

groups_miRNAseq <- groups
for(i in 2:5){
    groups_miRNAseq[[i]] <- as.factor(cutree(hc_miRNAseq, k = i))
}

png('miRNAseq.png')
pheatmap(lu$miRNAseq[,1:500], annotation_row = groups_miRNAseq,
         show_rownames = FALSE, show_colnames = FALSE, annotation_legend = FALSE, 
         cluster_rows = FALSE, clustering_distance_rows = "euclidean")
dev.off()

adjustedRandIndex(groups_miRNAseq[[4]], groups_mutation$Type) # 0.1102866

### COCA 

moc <- buildMOC(lu_reduced, M = 3, maxK = 10, methods = "hclust", 
         distances = c("euclidean", "binary", "euclidean"),  savePNG = TRUE)

pheatmap(moc$moc, cluster_cols = F)

datasetIndicator <- as.integer(c(rep(1,moc$K[1]), rep(2, moc$K[2]), rep(3, moc$K[3])))
datasetNames <- c(rep("RNAseq", moc$K[1]), rep("mutation", moc$K[2]), rep("miRNAseq", moc$K[3]))

load("annotations.RData")

# png("matrix-of-clusters.png")
plotMOC(moc$moc, datasetIndicator = datasetIndicator, datasetNames = datasetNames,
        annotations = groups, clc = T, clr = T)
# dev.off()

coca <- coca(moc$moc, maxK = 5, ccClMethod = 'hc', ccDistHC = "binary")

adjustedRandIndex(coca$clusterLabels, groups$Type)

groups$COCA <- as.factor(coca$clusterLabels)

png("matrix-of-clusters.png")
plotMOC(moc$moc, datasetIndicator = datasetIndicator, datasetNames = datasetNames,
        annotations = groups, clc = F, clr = T)
dev.off()
