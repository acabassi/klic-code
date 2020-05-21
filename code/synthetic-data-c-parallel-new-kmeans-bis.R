#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(mclust)
library(coca)
library(klic)

j <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### Data generation ### 
### C1. Three vs six clusters ###

n_types <- 2
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

N <- n_obs_per_cluster*n_clusters
P <- n_variables

Sigma <- diag(n_variables)

uno <- rep(1, n_obs_per_cluster)
cluster_labels_6clusters <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)
cluster_labels_3clusters <- c(uno, uno, uno*2, uno*2, uno*3, uno*3)

load("../data/synthetic-data-c-new.RData")

### Clustering one dataset at a time ###

# Initialise parameters for kernel k-means
parameters <- list(cluster_count = n_clusters)
ari_one_6clusters <- ari_one_3clusters <-
  cophenetic <-  rep(NA, n_types)
CM<- array(NA, c(N, N, n_types))

for(i in 1:n_types){
    # Use consensus clustering to find kernel matrix
    CM_temp <- consensusCluster(data[,,i,j], n_clusters/i, clMethod = "kmeans")
    # Shift the eigenvalues of the kernel matrix so that it is positive
    # semi-definite
    CM[,,i] <- spectrumShift(CM_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeans(CM[,,i], parameters)$clustering
    # Compute ARI
    ari_one_6clusters[i] <- adjustedRandIndex(
      cluster_labels_6clusters, kkmeans_labels)
    ari_one_3clusters[i] <- adjustedRandIndex(
      cluster_labels_3clusters, kkmeans_labels)
    # Compute cophenetic correlation coefficient
    cophenetic[i] <- copheneticCorrelation(CM[,,i])
}

##################################### K = 6 ####################################

# Run KLIC
parameters <- list(iteration_count = 1000, cluster_count = n_clusters)
klicOutput <- lmkkmeans(CM, parameters)

# Extract cluster labels and weights
klic_labels_6clusters <- klicOutput$clustering
weights_6clusters <- colMeans(klicOutput$Theta)

# Compute ARI
ari_all_6clusters <- adjustedRandIndex(
  klic_labels_6clusters, cluster_labels_6clusters) 

##################################### K = 3 ####################################

# Run KLIC
parameters <- list(iteration_count = 1000, cluster_count = n_clusters/2)
klicOutput <- lmkkmeans(CM, parameters)

# Extract cluster labels and weights
klic_labels_3clusters <- klicOutput$clustering
weights_3clusters <- colMeans(klicOutput$Theta)

# Compute ARI
ari_all_3clusters <- adjustedRandIndex(
  klic_labels_3clusters, cluster_labels_3clusters) 

### Save results ###
save(
  ari_one_6clusters, ari_all_6clusters, weights_6clusters, 
  ari_one_3clusters, ari_all_3clusters, weights_3clusters,
  cophenetic,
  file = paste0("../results/ari-c-", j,"-new-kmeans-bis.RData"))
