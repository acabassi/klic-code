#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(coca)
library(iCluster)
library(klic)
library(mclust) # for adjustedRandIndex
library(rdetools) # for rbfkernel

j <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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

uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

args <- commandArgs(trailingOnly=TRUE)
separation_level <- as.integer(args[1]) # Same for all datasets
load(paste0("../data/synthetic-data-a-sep", separation_level,".RData"))

### Clustering one dataset at a time ###

ari_one <- ari_one_rbfk <- rep(NA, n_datasets_same_rho)

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters
CM_rbfk <- array(NA, c(N, N, n_datasets_same_rho))

for(i in 1:n_datasets_same_rho){
    # Use consensus clustering to find kernel matrix
    CM_temp <- consensusCluster(data[,,i,j], n_clusters)
    # Use RBF kernel to obtain another kernel matrix
    CM_rbfk_temp <- rbfkernel(data[,,i,j], sigma = 1)
    # Make sure it is REALLY symmetric
    CM_rbfk_temp[lower.tri(CM_rbfk_temp)] = t(CM_rbfk_temp)[lower.tri(CM_rbfk_temp)]
    # Shift the eigenvalues of the kernel matrices so that they are positive
    # semi-definite
    CM_temp <- spectrumShift(CM_temp)
    CM_rbfk[,,i] <- spectrumShift(CM_rbfk_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeans(CM_temp, parameters)$clustering
    kkmeans_labels_rbfk <- kkmeans(CM_rbfk[,,i], parameters)$clustering
    # Compute ARI
    ari_one[i] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
    ari_one_rbfk[i] <- adjustedRandIndex(kkmeans_labels_rbfk, cluster_labels)
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

### Kernel k-means with RBF kernel ###
parameters$iteration_count <- 100
rbfk_cluster_labels <- lmkkmeans(CM_rbfk, parameters)$clustering
ari_all_rbfk <- adjustedRandIndex(rbfk_cluster_labels, cluster_labels)

### Save results ###
save(ari_one, ari_one_rbfk,
     ari_all,  ari_coca, ari_icluster, ari_all_rbfk,
     weights,
     file = paste0("../results/ari-a-", j,"-sep-", separation_level,".RData"))
