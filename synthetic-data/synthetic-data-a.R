#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(MASS)
library(ggplot2)
library(GGally)
library(reshape)
library(klic)
library(coca)
library(mclust)
library(reshape2)
library(ggplot2)

# Define ggplot2 theme 
my_theme <-  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.major.x = element_blank() ,
    panel.grid.major.y = element_line(size=.1, color="black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
)

### Data generation ### 
### A. Six clusters, same cluster separability in each dataset ###

n_experiments <- 100
n_separation_levels <- 1
n_datasets_same_rho <- 3
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

Sigma <- diag(n_variables)

N <- n_obs_per_cluster*n_clusters
P <- n_variables
data <- array(NA, c(N, P, n_datasets_same_rho, n_experiments))

separation_level <- 4 # Same for all datasets

set.seed(151)

for(i in 1:n_datasets_same_rho){
  for(j in 1:n_experiments){
    mu = rep(NA, N)
    for(k in 1:n_clusters){
      mu = rep((k*(separation_level-1)/2), n_variables)
      data[((k-1)*n_obs_per_cluster+1):(k*n_obs_per_cluster),,i,j] <- mvrnorm(n = n_obs_per_cluster, mu, Sigma)
    }
  }
}

uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

### Plot ###

dataset_example <- data[,,1,1]
example6_dataframe <- as.data.frame(dataset_example)
example6_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6_dataframe, # data.frame with variables
        title="Datasets with the same cluster separability",  # title of the plot
        mapping=ggplot2::aes(colour = cluster))
ggsave("synthetic-data-a.pdf", width = 7, height = 7)

save(data, file = "synthetic-data-a.RData")

### Clustering one dataset at a time ###

CM <- array(NA, c(n_obs_per_cluster*n_clusters, n_obs_per_cluster*n_clusters,
                  n_datasets_same_rho, n_experiments))

ari_one <- array(NA, c(n_datasets_same_rho, n_experiments))

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters

for(i in 1:n_datasets_same_rho){
  for(j in 1:n_experiments){
    # Use consensus clustering to find kernel matrix
    CM_temp <- consensusCluster(data[,,i,j], n_clusters)
    # Shift the eigenvalues of the kernel matrix so that it is positive semi-definite
    CM_temp <- spectrumShift(CM_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeans(CM_temp, parameters)$clustering
    # Compute ARI
    ari_one[i,j] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
  }
}

### All together ###

ari_all <- rep(NA, n_experiments)
weights <- array(NA, c(n_datasets_same_rho, n_experiments))

data_for_klic <- list()
for(j in 1:n_experiments){
  
  # Build list of datasets, input for KLIC
  for(i in 1:n_datasets_same_rho){
    data_for_klic[[i]] <- data[,,i,j]
  }
  
  # Run KLIC
  klicOutput <- klic(data_for_klic, n_datasets_same_rho,
                      individualK = rep(n_clusters, n_datasets_same_rho),
                      globalK = n_clusters)
  
  # Extract cluster labels and weights
  klic_labels <- klicOutput$globalClusterLabels
  weights[,j] <- colMeans(klicOutput$weights)
  
  # Compute ARI
  ari_all[j] <- adjustedRandIndex(klic_labels, cluster_labels) 
}

### COCA ###

ari_coca <- rep(NA, n_experiments)
moc <- array(NA, c(dim(data)[1],n_clusters*n_datasets_same_rho))

for(i in 1:n_experiments){
  
  # Over-write the matrix-of-clusters each time
  count <- 0
  
  # Find clusters in each dataset
  for(j in 1:n_datasets_same_rho){
    
    kmeans_cluster_labels <- kmeans(data[,,j,i], n_clusters)$cluster
    
    # Fill the matrix-of-clusters
    for(l in 1:n_clusters){
      count <- count + 1
      moc[which(kmeans_cluster_labels == l), count] <- 1
      moc[which(kmeans_cluster_labels != l), count] <- 0
    }
  }
  
  # Use COCA to find final clusters
  coca_cluster_labels <- coca(moc, n_clusters)$clusterLabels
  # Compute ARI 
  ari_coca[i] <- adjustedRandIndex(cluster_labels, coca_cluster_labels)
}

# Save results
save(ari_one, ari_all, weights, ari_coca, file = "ari-a.RData")

# Load results
load("ari-a.RData")

dim(ari_one)
dim(as.matrix(ari_all))
ari <- cbind(t(ari_one), as.matrix(ari_all))
colnames(ari) <- c("A", "B", "C", "A+B+C")

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")
ari.m$Datasets <- factor(ari.m$Datasets,
                              levels = c("A","B", "C", "A+B+C"), ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) + geom_boxplot() + ylim(0,1) + my_theme
ggsave("ari-a.pdf", width = 15, height = 10, units = "cm")

dim(weights)
rownames(weights) <- c("A", "B", "C")
weights.m <- melt(t(weights))
head(weights.m)
colnames(weights.m) <- c("Experiment", "Dataset", "Weight")
ggplot(data = weights.m, aes(x=Dataset, y=Weight)) + geom_boxplot() + ylim(0,1) + my_theme
ggsave("weights-a.pdf", width = 12, height = 10, units = "cm")
