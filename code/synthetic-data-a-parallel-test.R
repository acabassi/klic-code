#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(coca)
library(klic)
library(mclust)
library(sparcl)

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
load("../data/RBF_sigma_values.RData")

### Clustering one dataset at a time ###

ari_one <- ari_one_binary <- ari_one_sparse <-
  coph <- coph_binary <- coph_sparse <- 
  rep(NA, n_datasets_same_rho)

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters
CM <- CM_binary <- CM_sparse <- array(NA, c(N, N, n_datasets_same_rho))

for(i in 1:n_datasets_same_rho){
    # Use consensus clustering to find kernel matrix
    CM_temp <-
      consensusCluster(data[,,i,j], n_clusters, clMethod = "kmeans")
    CM_sparse_temp <-
      consensusCluster(data[,,i,j], n_clusters, clMethod = "sparse-kmeans")
    
    kmeans_cluster_labels <- kmeans(data[,,i,j], 6)$cluster
    CM_binary_temp <- matrix(0, N, N)
    for(index_1 in 1:N){
      for(index_2 in 1:N){
        if(kmeans_cluster_labels[index_1] == kmeans_cluster_labels[index_2])
          CM_binary_temp[index_1, index_2] <- 1 
      }
    }
    
    # Shift the eigenvalues of the kernel matrices so that they are positive
    # semi-definite
    CM[,,i] <- spectrumShift(CM_temp) 
    CM_binary[,,i] <- spectrumShift(CM_binary_temp)
    CM_sparse[,,i] <- spectrumShift(CM_sparse_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <-
      kkmeans(CM[,,i], parameters)$clustering
    kkmeans_labels_binary <- 
      kkmeans(CM_binary[,,i], parameters)$clustering
    kkmeans_labels_sparse <- 
      kkmeans(CM_sparse[,,i], parameters)$clustering
    # Compute ARI
    ari_one[i] <-
      adjustedRandIndex(kkmeans_labels, cluster_labels)
    ari_one_binary[i] <-
      adjustedRandIndex(kkmeans_labels_binary, cluster_labels)
    ari_one_sparse[i] <-
      adjustedRandIndex(kkmeans_labels_sparse, cluster_labels)
    # Cophenetic correlation coefficients
    coph[i] <- copheneticCorrelation(CM[,,i])
    coph_binary[i] <- copheneticCorrelation(CM_binary[,,i])
    coph_sparse[i] <- copheneticCorrelation(CM_sparse[,,i])
}

### All together ###
data_for_klic <- list()

# Build list of datasets, input for KLIC
for(i in 1:n_datasets_same_rho){
  data_for_klic[[i]] <- data[,,i,j]
}

# Run KLIC
klicOutput_sparse <- tryCatch(klic(data_for_klic, n_datasets_same_rho,
                            individualK = rep(n_clusters, n_datasets_same_rho),
                            globalK = n_clusters, ccClMethods = "sparse-kmeans",
                            C = 1000),
                       error = function(err)
                         list(globalClusterLabels = rep(NA, 300),
                              weights = NA))

klicOutput <- tryCatch(klic(data_for_klic, n_datasets_same_rho,
                       individualK = rep(n_clusters, n_datasets_same_rho),
                       globalK = n_clusters, ccClMethods = "kmeans",
                       C = 1000),
                       error = function(err)
                         list(globalClusterLabels = rep(NA, 300),
                              weights = NA))

parameters$iteration_count <- 1000
klicOutput2 <- lmkkmeans(CM, parameters)
klicOutput_sparse2 <- lmkkmeans(CM_sparse, parameters)
klicOutput_binary <- tryCatch(lmkkmeans(CM_binary, parameters),
                              error = function(err)
                                list(clustering = rep(NA, 300),
                                     Theta = NA))

# Extract cluster labels and weights
klic_labels <- klicOutput$globalClusterLabels
klic_labels2 <- klicOutput2$clustering
klic_labels_binary <- klicOutput_binary$clustering
klic_labels_sparse <- klicOutput_sparse$globalClusterLabels
klic_labels_sparse2 <- klicOutput_sparse2$clustering
weights <- klicOutput$weights
weights_binary <- klicOutput_binary$Theta
weights_sparse <- klicOutput_sparse$weights

# Compute ARI
ari_all <- adjustedRandIndex(klic_labels, cluster_labels) 
all_all2 <- adjustedRandIndex(klic_labels2, cluster_labels)
ari_all_binary <- adjustedRandIndex(klic_labels_binary, cluster_labels)
ari_all_sparse <- adjustedRandIndex(klic_labels_sparse, cluster_labels)
ari_all_sparse2 <- adjustedRandIndex(klic_labels_sparse2, cluster_labels)

### COCA ###

moc <- moc_sparcl <- array(NA, c(dim(data)[1], n_clusters*n_datasets_same_rho))

# Over-write the matrix-of-clusters each time
count <- 0

ari_kmeans <- rep(NA, n_datasets_same_rho)
ari_kmeans_sparse <- rep(NA, n_datasets_same_rho)

# Find clusters in each dataset
for(i in 1:n_datasets_same_rho){
  
  kmeans_cluster_labels <- kmeans(data[,,i,j], n_clusters)$cluster
  sparse_kmeans_cluster_labels <- KMeansSparseCluster(data[,,i,j], n_clusters,
                                                      wbounds = sqrt(2))[[1]]$Cs
  
  ari_kmeans[i] <- adjustedRandIndex(cluster_labels, kmeans_cluster_labels)
  ari_kmeans_sparse[i] <- adjustedRandIndex(cluster_labels, sparse_kmeans_cluster_labels)
  
  # Fill the matrix-of-clusters
  for(l in 1:n_clusters){
    count <- count + 1
    moc[which(kmeans_cluster_labels == l), count] <- 1
    moc[which(kmeans_cluster_labels != l), count] <- 0
    
    moc_sparcl[which(sparse_kmeans_cluster_labels == l), count] <- 1
    moc_sparcl[which(sparse_kmeans_cluster_labels != l), count] <- 0
  }
}

# Use COCA to find final clusters
coca_cluster_labels <- coca(moc, n_clusters)$clusterLabels
coca_cluster_labels_sparcl <- coca(moc_sparcl, n_clusters)$clusterLabels
# Compute ARI 
ari_coca <- adjustedRandIndex(cluster_labels, coca_cluster_labels)
ari_coca_sparcl <- adjustedRandIndex(cluster_labels, coca_cluster_labels_sparcl)

### Weighted kernels ###

weighted_kernel <- weighted_kernel_binary <- weighted_kernel_sparse <-
  matrix(0, N, N) 
for (i in 1:n_datasets_same_rho) {
  if(!is.na(weights)){
    weighted_kernel <- weighted_kernel +
      (weights[, i] %*% t(weights[, i])) * CM[,, i]
  }
  
  if(!is.na(weights_binary)){
  weighted_kernel_binary <- weighted_kernel_binary +
    (weights_binary[, i] %*% t(weights_binary[, i])) * CM_binary[,, i]
  }
  
  if(!is.na(weights_sparse)){
    weighted_kernel_sparse <- weighted_kernel_sparse + 
      (weights_sparse[, i] %*% t(weights_sparse[, i])) * CM_sparse[,,i]
  }
}

### Save results ###
save(ari_one, ari_one_binary, ari_one_sparse,
     coph, coph_binary, coph_sparse,
     ari_all,  ari_all_binary, ari_all_sparse,
     ari_all2, ari_sparse2,
     ari_coca, ari_coca_sparcl,
     moc, moc_sparcl,
     weights, weights_binary, weights_sparse,
     weighted_kernel, weighted_kernel_binary, weighted_kernel_sparse,
     file = paste0("../results/ari-a-test-", j,"-sep-",
                   separation_level, ".RData"))
