#################################
### ARI of breast cancer data ###
#################################

rm(list=ls())

# Load clusterings found with KLIC
#load("ukic.RData")
load("coca.RData")
#load("skic2.RData")

# IDs of the tumour samples..
GE = read.csv("context1_GE.csv", header = TRUE)
namesExp = names(GE)[2:349]
namesExp = substr(namesExp,1,12)
namesExp = gsub(".", "-", namesExp, fixed = TRUE)

# UKICclusters <- as.matrix(UKICclusters)
# rownames(UKICclusters)<- namesExp

COCAclusters <- as.matrix(COCAclusters$clusterLabels)
rownames(COCAclusters) <- namesExp

# SKICclusters <- as.matrix(SKICclusters)
# rownames(SKICclusters) <- namesExp

survivalBRCA <- read.csv("Table1Nature.csv")
names(survivalBRCA)

# prepare the data to make sure they match the two clusters from above
PAM50 <- survivalBRCA$PAM50.mRNA
PAM50 <- as.matrix(PAM50)
rownames(PAM50) <- survivalBRCA$Complete.TCGA.ID
PAM50 <- PAM50[rownames(COCAclusters),]

PAM50bi <- rep(NA, 348)

library(R.matlab)
table(PAM50)
PAM50bi[PAM50 == 'Basal-like'] = 1
PAM50bi[PAM50 == 'HER2-enriched'] = 2
PAM50bi[PAM50 == 'Luminal A'] = 2
PAM50bi[PAM50 == 'Luminal B'] = 2
PAM50bi[PAM50 == 'Normal-like'] = 2
table(PAM50bi)

PAM50[PAM50 == 'Basal-like'] = 1
PAM50[PAM50 == 'HER2-enriched'] = 2
PAM50[PAM50 == 'Luminal A'] = 3
PAM50[PAM50 == 'Luminal B'] = 4
PAM50[PAM50 == 'Normal-like'] = 5
table(PAM50)

save(PAM50bi, file = "PAM50_binary.RData")
save(PAM50, file = "PAM50.RData")

library(mclust)
# adjustedRandIndex(PAM50bi, UKICclusters)
adjustedRandIndex(PAM50bi, COCAclusters)
# adjustedRandIndex(UKICclusters, COCAclusters)
# adjustedRandIndex(PAM50bi, SKICclusters)
# adjustedRandIndex(UKICclusters, SKICclusters)
