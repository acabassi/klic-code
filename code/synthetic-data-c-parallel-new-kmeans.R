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
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

load("../data/synthetic-data-c-new.RData")

##################################### K = 6 ####################################

### Clustering one dataset at a time ###

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters
ari_one_6clusters <- cophenetic_6clusters <- rep(NA, n_types)
CM_6clusters <- array(NA, c(N, N, n_types))

for(i in 1:n_types){
    # Use consensus clustering to find kernel matrix
    CM_temp <- consensusCluster(data[,,i,j], n_clusters, clMethod = "kmeans")
    # Shift the eigenvalues of the kernel matrix so that it is positive
    # semi-definite
    CM_6clusters[,,i] <- spectrumShift(CM_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeans(CM_6clusters[,,i], parameters)$clustering
    # Compute ARI
    ari_one_6clusters[i] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
    cophenetic_6clusters[i] <- copheneticCorrelation(CM_6clusters[,,i])
}

### Combining subsets of three datasets ###

# Build list of datasets, input for KLIC
data_for_klic <- list()
for(l in 1:n_types){
  data_for_klic[[l]] <- data[,,l,j]
}

# Run KLIC
klicOutput <- klic(data_for_klic, n_types,
                   individualK = rep(n_clusters, n_types),
                   globalK = n_clusters
                   )

# Extract cluster labels and weights
klic_labels <- klicOutput$globalClusterLabels
weights_6clusters <- colMeans(klicOutput$weights)

# Compute ARI
ari_all_6clusters <- adjustedRandIndex(klic_labels, cluster_labels) 

##################################### K = 3 ####################################

n_clusters <- n_clusters/2

### Clustering one dataset at a time ###

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters
ari_one_3clusters <- cophenetic_3clusters <- rep(NA, n_types)
CM_3clusters <- array(NA, c(N, N, n_types))

for(i in 1:n_types){
  # Use consensus clustering to find kernel matrix
  CM_temp <- consensusCluster(data[,,i,j], n_clusters, clMethod = "kmeans")
  # Shift the eigenvalues of the kernel matrix so that it is positive
  # semi-definite
  CM_3clusters[,,i] <- spectrumShift(CM_temp)
  # Use kernel k-means to find clusters
  kkmeans_labels <- kkmeans(CM_3clusters[,,i], parameters)$clustering
  # Compute ARI
  ari_one_3clusters[i] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
  cophenetic_3clusters[i] <- copheneticCorrelation(CM_3clusters[,,i])
}

### Combining the two datasets ###

# Build list of datasets, input for KLIC
data_for_klic <- list()
for(l in 1:n_types){
  data_for_klic[[l]] <- data[,,l,j]
}

# Run KLIC
klicOutput <- klic(data_for_klic, n_types,
                   individualK = rep(n_clusters, n_types),
                   globalK = n_clusters, ccClMethods = "kmeans")

# Extract cluster labels and weights
klic_labels <- klicOutput$globalClusterLabels
weights_3clusters <- colMeans(klicOutput$weights)

# Compute ARI
ari_all_3clusters <- adjustedRandIndex(klic_labels, cluster_labels) 

### Save results ###
save(
  ari_one_6clusters, ari_all_6clusters, weights_6clusters, cophenetic_6clusters,
  ari_one_3clusters, ari_all_3clusters, weights_3clusters, cophenetic_3clusters,
  file = paste0("../results/ari-c-", j,"-new-kmeans.RData"))
