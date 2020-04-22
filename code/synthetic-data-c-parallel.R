#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(mclust)
library(coca)
library(klic)

j <- 1 # as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### Data generation ### 
### C1. Three vs six clusters ###

n_types <- 3
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

N <- n_obs_per_cluster*n_clusters
P <- n_variables

Sigma <- diag(n_variables)

uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

load("../data/synthetic-data-c.RData")

### Clustering one dataset at a time ###

ari_one <- rep(NA, n_types)

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters

for(i in 1:n_types){
    # Use consensus clustering to find kernel matrix
    CM_temp <- consensusCluster(data[,,i,j], n_clusters)
    # Shift the eigenvalues of the kernel matrix so that it is positive
    # semi-definite
    CM_temp <- spectrumShift(CM_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeans(CM_temp, parameters)$clustering
    # Compute ARI
    ari_one[i] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
}

### Combining subsets of three datasets ###

n_subsets <- 2
n_datasets_per_subset <- 2

ari_all <- rep(NA, n_subsets)
weights <- array(NA, c(n_datasets_per_subset, n_subsets))

subsets <- rbind(c(1,2),c(2,3))

data_for_klic <- list()

for(i in 1:n_subsets){
  
  # Build list of datasets, input for KLIC
  count <- 0
  datasets_in_subset <- subsets[i,]
  for(l in datasets_in_subset){
    count <- count + 1
    data_for_klic[[count]] <- data[,,l,j]
  }
  
  # Run KLIC
  klicOutput <- klic(data_for_klic, n_datasets_per_subset,
                     individualK = rep(n_clusters, n_datasets_per_subset),
                     globalK = n_clusters)
  
  # Extract cluster labels and weights
  klic_labels <- klicOutput$globalClusterLabels
  weights[,i] <- colMeans(klicOutput$weights)
  
  # Compute ARI
  ari_all[i] <- adjustedRandIndex(klic_labels, cluster_labels) 
}

### Save results ###
save(ari_one, ari_all, weights,
     file = paste0("../results/ari-c-", j,".RData"))