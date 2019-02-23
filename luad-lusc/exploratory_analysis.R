##########################################
### LUAD vs LUSC: Exploratory analysis ###
##########################################

rm(list=ls())
setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/klic-code/luad-lusc/")

load("TCGA_LUAD_match.RData")
load("TCGA_LUSC_match.RData")

library(pheatmap)
library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)

dim(luad$RNA_seq) # 20502 X 59
dim(luad$mutation) # 13932 X 59 - binary data 
dim(luad$miRNAseq) # 1046 X 59

dim(lusc$RNA_seq) # 20502 X 52
dim(lusc$mutation) #  13665 X 52
dim(lusc$miRNAseq) # 1046 X 52

sum(rownames(luad$RNA_seq)%in%rownames(lusc$RNA_seq)) # 20502 
sum(rownames(luad$mutation)%in%rownames(lusc$mutation)) # 11197 - Should I take only the common ones or not? 
sum(rownames(luad$miRNAseq)%in%rownames(lusc$miRNAseq)) # 1046

lu <- list(RNA_seq = rbind(t(luad$RNA_seq), t(lusc$RNA_seq)), 
           mutation = rbind(t(luad$mutation[rownames(luad$mutation)%in%rownames(lusc$mutation),]), 
                            t(lusc$mutation[rownames(lusc$mutation)%in%rownames(luad$mutation),])),
           miRNAseq = rbind(t(luad$miRNAseq), t(lusc$miRNAseq)),
           label = c(rep("LUAD", 59), rep("LUSC", 52)))

lu$RNA_seq <- lu$RNA_seq[,-which(colSums(lu$RNA_seq)==0)]
lu$mutation <- lu$mutation[,-which(colSums(lu$miRNAseq)==0)]
lu$miRNAseq <- lu$miRNAseq[,-which(colSums(lu$miRNAseq)==0)]

save(lu, file = "merged-data.RData")

### PCA
pca_RNAseq <- prcomp(lu$RNA_seq, center = TRUE, scale = TRUE)
ggbiplot(pca_RNAseq, var.axes = F, groups = lu$label)
ggsave("pca_RNAseq.png")

# Here I try to do PCA on the mutation but IT IS BINARY DATA
pca_mutation <- prcomp(lu$mutation, center = FALSE, scale = FALSE)
ggbiplot(pca_mutation, var.axes = F, groups = lu$label)
ggsave("pca_mutation.png")

pca_miRNAseq <- prcomp(lu$miRNAseq, center = TRUE, scale = TRUE)
ggbiplot(pca_miRNAseq, var.axes = F, groups = lu$label)
ggsave("pca_miRNAseq.png")
