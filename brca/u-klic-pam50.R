###################################################
### Clustering BRCA data with unsupervised KLIC ###
###################################################

rm(list = ls())

setwd("~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/klic-code/brca")

# library(devtools)
# install_github("acabassi/klic")
library(klic)
library(mclust)

### Load datasets ###
load("Breast50.RData")
X <- X50

# Transpose and scale datasets
for(i in 1:4){
  X[[i]] <- scale(t(X[[i]]))
}

### Load number of clusters ###
load("coca50_num_clusters.RData")

### KLIC

klicOutput <- klic(X, M = 4, individualK = num_clusters_each, globalK = k_final,
                   savePlots = FALSE, fileName = "test")

# Save output
save(klicOutput, file = "u-klic-50.RData")

###### Plots ######

### Load PAM50 subtype ###
survivalBRCA <- read.csv("Table1Nature.csv")
PAM50 <- survivalBRCA$PAM50.mRNA
names(PAM50) <- survivalBRCA$Complete.TCGA.ID
patientIDs <- rownames(X[[1]])
patientIDs <- substr(patientIDs,1,12)
patientIDs <- gsub(".", "-", patientIDs, fixed = TRUE)
PAM50 <- PAM50[patientIDs]
table(PAM50)
PAM50 <- as.integer(factor(PAM50))
table(PAM50)

### Plot data matrix with PAM50 subtypes ###
library(heatmaply)
PAM50 <- as.matrix(PAM50)
colnames(PAM50) <- "PAM50"
heatmaply(X[[2]], scale = 'column', trace = 'none', col = viridis(100),
          RowSideColors = PAM50, key = TRUE,
          dendrogram = 'row',
          label_names = c('Individual', 'Variable', 'Value'),
          xlab = "Variables", ylab = "Individuals",
          showticklabels = c(FALSE, FALSE), fontsize_row = 12,
          plot_method = "plotly", colorbar_xanchor = "left")

### Plot consensus matrices ###

for(i in 1:4){
  plotSimilarityMatrix(klicOutput$consensusMatrices[,,i], y = as.integer(factor(PAM50)), clusLabels = PAM50, file_name = paste("cm",i,"_PAM50_50.png", sep=""), savePNG = TRUE, semiSupervised = TRUE, myLegend = c("Basal-like",'HER2-enriched','Luminal A','Luminal B','Normal-like'))
  plotSimilarityMatrix(klicOutput$consensusMatrices[,,i], y = as.integer(factor(PAM50)), clusLabels = klicOutput$globalClusterLabels, file_name = paste("cm",i,"_clusters_50.png", sep=""), savePNG = TRUE, semiSupervised = TRUE, myLegend = c("Basal-like",'HER2-enriched','Luminal A','Luminal B','Normal-like'))
}

plotSimilarityMatrix(klicOutput$weightedKM, y = as.integer(factor(PAM50)), clusLabels = klicOutput$globalClusterLabels, file_name = "weightedCM_clusters_50.png", savePNG = TRUE, semiSupervised = TRUE, myLegend = c("Basal-like",'HER2-enriched','Luminal A','Luminal B','Normal-like'))
plotSimilarityMatrix(klicOutput$weightedKM, y = (PAM50==1)+1, clusLabels = klicOutput$globalClusterLabels, file_name = "weightedCM_binaryInd_50.png", savePNG = TRUE, semiSupervised = TRUE, myLegend = c("Basal-like",'Other'))
