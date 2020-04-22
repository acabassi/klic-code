#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(mclust) # for adjustedRandIndex
library(klic)
library(coca)
library(iCluster)

j <- 1 # as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

### Load data ### 
### A. Six clusters, same cluster separability in each dataset ###

n_separation_levels <- 1
n_datasets_same_rho <- 4
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

Sigma <- diag(n_variables)

N <- n_obs_per_cluster*n_clusters
P <- n_variables

separation_level <- 4 # Same for all datasets

uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

separation_level <- as.integer(args[1]) # Same for all datasets
load(paste0("../data/synthetic-data-a-sep", separation_level,".RData"))

### Clustering one dataset at a time ###

ari_one <- rep(NA, n_datasets_same_rho)

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters

for(i in 1:n_datasets_same_rho){
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

### All together ###
data_for_klic <- list()

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
weights <- colMeans(klicOutput$weights)

# Compute ARI
ari_all <- adjustedRandIndex(klic_labels, cluster_labels) 

### COCA ###

moc <- array(NA, c(dim(data)[1],n_clusters*n_datasets_same_rho))

# Over-write the matrix-of-clusters each time
count <- 0

# Find clusters in each dataset
for(i in 1:n_datasets_same_rho){
  
  kmeans_cluster_labels <- kmeans(data[,,i,j], n_clusters)$cluster
  
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
ari_coca <- adjustedRandIndex(cluster_labels, coca_cluster_labels)

### iCluster ###

# Put data into a list of datasets
data_iCluster <- list()
for(i in 1:n_datasets_same_rho){
  data_iCluster[[i]] <- data[,,i,j]
}

# Use iCluster to find cluster labels
icluster_labels <- iCluster2(data_iCluster, n_clusters)$clusters
# Compute ARI
ari_icluster <- adjustedRandIndex(cluster_labels, icluster_labels)

### Save results ###
save(ari_one, ari_all, weights, ari_coca, ari_icluster,
     file = paste0("../results/ari-a-", j,".RData"))
