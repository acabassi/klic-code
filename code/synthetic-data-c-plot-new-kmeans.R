#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(MASS)
library(ggplot2)
library(reshape)
library(reshape2)

### Define ggplot2 theme ###
my_theme <-  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.major.x = element_blank() ,
    panel.grid.major.y = element_line(size=.1, color="black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
)

### Load results ###
n_experiments <- 100
all_weights_6clusters <- all_weights_3clusters <- 
  all_ari_one_6clusters <- all_ari_one_3clusters <- 
  all_cophene_6clusters <- all_cophene_3clusters <-
  array(NA, c(2, n_experiments))
all_ari_all_6clusters <- all_ari_all_3clusters <- rep(NA, n_experiments)


for(j in 1:n_experiments){
  load(paste0("../results/ari-c-", j,"-new-kmeans.RData"))
  
  # K = 6
  all_weights_6clusters[,j] <- weights_6clusters
  all_ari_all_6clusters[j] <- ari_all_6clusters
  all_ari_one_6clusters[,j] <- ari_one_6clusters
  all_cophene_6clusters[,j] <- cophenetic_6clusters
  
  # K = 3
  all_weights_3clusters[,j] <- weights_3clusters
  all_ari_all_3clusters[j] <- ari_all_3clusters
  all_ari_one_3clusters[,j] <- ari_one_3clusters
  all_cophene_3clusters[,j] <- cophenetic_3clusters
}

##################################### K = 3 ####################################

### Adjusted Rand index

ari <- rbind(all_ari_one_3clusters, all_ari_all_3clusters)
rownames(ari) <- c("6", "3", "3+6")
ari<-t(ari)

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")

ari.m$Datasets <- factor(ari.m$Datasets,
                       levels = c("3","6", "3+6"), ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme
ggsave("../figures/ari-c-new-3clusters-kmeans.jpg",
       device = "jpeg", width = 7, height = 8, units = "cm")

### Weights 

rownames(all_weights_3clusters) <- c("6", "3")
colnames(all_weights_3clusters) <- 1:n_experiments

weights_3clusters.m <- melt(t(all_weights_3clusters))
head(weights_3clusters.m) 
colnames(weights_3clusters.m) <- c("Experiments", "Dataset", "Weight")
weights_3clusters.m$Dataset <- factor(weights_3clusters.m$Dataset,
                         levels = c("3","6"), ordered = TRUE)


ggplot(data = weights_3clusters.m, aes(x=Dataset, y=Weight)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme
ggsave("../figures/weights-3clusters-new-kmeans.jpg", device = "jpeg",
       width = 3.5, height = 8, units = "cm")

### Cophenetic correlation coefficient

rownames(all_cophene_3clusters) <- c("6", "3")
colnames(all_cophene_3clusters) <- 1:n_experiments

cophene_3clusters.m <- melt(t(all_cophene_3clusters))
head(cophene_3clusters.m) 
colnames(cophene_3clusters.m) <- c("Experiments", "Dataset", "Correlation")
cophene_3clusters.m$Dataset <- factor(cophene_3clusters.m$Dataset,
                                      levels = c("3","6"), ordered = TRUE)


ggplot(data = cophene_3clusters.m, aes(x=Dataset, y=Correlation)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme
ggsave("../figures/cophene-3clusters-new-kmeans.jpg", device = "jpeg",
       width = 3.5, height = 8, units = "cm")

##################################### K = 6 ####################################

### Adjusted Rand index

ari <- rbind(all_ari_one_6clusters, all_ari_all_6clusters)
rownames(ari) <- c("6", "3", "3+6")
ari<-t(ari)

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")

ari.m$Datasets <- factor(ari.m$Datasets,
                         levels = c("3","6", "3+6"), ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme
ggsave("../figures/ari-c-6clusters-new-kmeans.jpg",
       device = "jpeg", width = 7, height = 8, units = "cm")

### Weights 

rownames(all_weights_6clusters) <- c("6", "3")
colnames(all_weights_6clusters) <- 1:n_experiments

weights_6clusters.m <- melt(t(all_weights_6clusters))
head(weights_6clusters.m) 
colnames(weights_6clusters.m) <- c("Experiments", "Dataset", "Weight")
weights_6clusters.m$Dataset <- factor(weights_6clusters.m$Dataset,
                                      levels = c("3","6"), ordered = TRUE)


ggplot(data = weights_6clusters.m, aes(x=Dataset, y=Weight)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme
ggsave("../figures/weights-6clusters-new-kmeans.jpg", device = "jpeg",
       width = 3.5, height = 8, units = "cm")

### Cophenetic correlation coefficient

rownames(all_cophene_6clusters) <- c("6", "3")
colnames(all_cophene_6clusters) <- 1:n_experiments

cophene_6clusters.m <- melt(t(all_cophene_6clusters))
head(cophene_6clusters.m) 
colnames(cophene_6clusters.m) <- c("Experiments", "Dataset", "Correlation")
cophene_6clusters.m$Dataset <- factor(cophene_6clusters.m$Dataset,
                                      levels = c("3","6"), ordered = TRUE)


ggplot(data = cophene_6clusters.m, aes(x=Dataset, y=Correlation)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme
ggsave("../figures/cophene-6clusters-new-kmeans.jpg", device = "jpeg",
       width = 3.5, height = 8, units = "cm")
