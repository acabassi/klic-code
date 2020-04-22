#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

rm(list = ls())

library(klic)
library(coca)
library(mclust)

j <- 1 # as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### Data generation ### 
### B. Six clusters, each dataset has a different level of cluster ###
###    separability ###

n_experiments <- 100
n_separation_levels <- 4
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

Sigma <- diag(n_variables)

N <- n_obs_per_cluster*n_clusters
P <- n_variables

# True cluster labels 
uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

load("../synthetic-data/synthetic-data-b.RData")

### Clustering one dataset at a time ###

ari_one <- rep(NA, n_separation_levels)

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count_m <- n_clusters

for(i in 1:n_separation_levels){
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

n_subsets <- 4
n_datasets_per_subset <- 3

ari_all <- rep(NA, c(n_subsets))
weights <- array(NA, c(n_datasets_per_subset, n_subsets))

subsets <- rbind(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4))

data_for_klic <- list()
for(i in 1:n_subsets){
  
  # Build list of datasets, input for KLIC
  count_m <- 0
  datasets_in_subset <- subsets[i,]
  for(l in datasets_in_subset){
    count_m <- count_m + 1
    data_for_klic[[count_m]] <- data[,,l,j]
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

### COCA & iCluster ###
ari_coca <- ari_icluster <- rep(NA, n_subsets)
moc <- array(NA, c(dim(data)[1], n_clusters*n_datasets_per_subset))

for(i in 1:n_subsets){
  
  # Select datasets 
  datasets_in_subset <- subsets[i,]
  
  # For iCluster
  data_iCluster <- list()
  
  # Over-write the matrix-of-clusters each time
  count_m <- 0
  count_l <- 0
  
  # Find clusters in each dataset
  for(l in datasets_in_subset){
    count_l <- count_l + 1
    data_iCluster[[count_l]] <- data[,,l,j] 
    
    kmeans_cluster_labels <- kmeans(data[,,l,j], n_clusters)$cluster
    
    # Fill the matrix-of-clusters
    for(m in 1:n_clusters){
      count_m <- count_m + 1
      moc[, count_m] <- (kmeans_cluster_labels == m)*1
    }
  }
  
  # Use COCA to find final clusters
  coca_cluster_labels <- coca(moc, n_clusters)$clusterLabels
  # Compute ARI 
  ari_coca[i] <- adjustedRandIndex(cluster_labels, coca_cluster_labels)
  
  # Use iCluster to find final clusters
  icluster_labels <- iCluster2(data_iCluster, n_clusters)$clusters
  # Compute ARI
  ari_icluster[i] <- adjustedRandIndex(cluster_labels, icluster_labels)

}

### Save results ###
save(ari_one, ari_all, weights, ari_coca, ari_icluster,
     file = paste0("../results/ari-b-", j,".RData"))
