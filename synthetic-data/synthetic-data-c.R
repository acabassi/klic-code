#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(MASS)
library(ggplot2)
library(GGally)
library(reshape)

### Data generation ### 
### C1. Six clusters ###

n_experiments <- 100
n_separation_levels <- 1
n_types <- 3
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

N <- n_obs_per_cluster*n_clusters
P <- n_variables
data <- array(NA, c(N, P, n_types, n_experiments))

Sigma <- diag(n_variables)

i = 3 # Medium separation level
for(j in 1:n_experiments){
  mu = rep(NA, N)
  for(k in 1:n_clusters){
    for(t in 1:(n_types-1)){
      for(k in 1:(n_clusters/t)){
          mu = rep(k*(i-1)/2, n_variables)
      data[((k-1)*n_obs_per_cluster*t+1):(k*n_obs_per_cluster*t),,t,j] <- mvrnorm(n = n_obs_per_cluster*t, mu, Sigma)
      }
    }
  }
}

i = 9 # Higher separation level
t = 1 # Only for dataset containing 6 clusters
for(j in 1:n_experiments){
  mu = rep(NA, N)
  for(k in 1:n_clusters){
      for(k in 1:(n_clusters/t)){
        mu = rep(k*(i-1)/2, n_variables)
        data[((k-1)*n_obs_per_cluster*t+1):(k*n_obs_per_cluster*t),,3,j] <- mvrnorm(n = n_obs_per_cluster*t, mu, Sigma)
      }
    }
  }

### Plots ###

t = 1 # 6 clusters
uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

dataset_example <- data[,,t,1]
example6_dataframe <- as.data.frame(dataset_example)
example6_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6_dataframe, # data.frame with variables
        title="6 clusters",  # title of the plot
        mapping=ggplot2::aes(colour = cluster))
ggsave("synthetic-data-c1.pdf", width = 7, height = 7)

t = 2 # 3 clusters
dataset_example <- data[,,t,1]
example3_dataframe <- as.data.frame(dataset_example)
example3_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example3_dataframe, # data.frame with variables
        title="3 clusters", # title of the plot
        mapping = ggplot2::aes(colour= cluster))
ggsave("synthetic-data-c2.pdf", width = 7, height = 7)

t = 3 # 6 clusters, well separated
dataset_example <- data[,,t,1]
example6plus_dataframe <- as.data.frame(dataset_example)
example6plus_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6plus_dataframe, # data.frame with variables
        title="6 clusters, well separated", # title of the plot
        mapping = ggplot2::aes(colour= cluster))
ggsave("synthetic-data-c3.pdf", width = 7, height = 7)

save(data, file = "synthetic-data-c.RData")

### Clustering one dataset at a time ###

library(coca)

ari_one <- array(NA, c(n_separation_levels, n_experiments))

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters

for(i in 1:n_types){
  for(j in 1:n_experiments){
    # Use consensus clustering to find kernel matrix
    CM_temp <- consensusCluster(data[,,i,j], n_clusters)
    # Shift the eigenvalues of the kernel matrix so that it is positive semi-definite
    CM_temp <- spectrumShift(CM_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeansTrain(CM_temp, parameters)$clustering
    # Compute ARI
    ari_one[i,j] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
  }
}

### Combining subsets of three datasets ###

n_subsets <- 2
n_datasets_per_subset <- 2

library(klic)
library(mclust)

ari_all <- array(NA, c(n_subsets, n_experiments))
weights <- array(NA, c(n_datasets_per_subset, n_subsets, n_experiments))

subsets <- rbind(c(1,2),c(2,3))

data_for_klic <- list()

for(j in 1:n_experiments){
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
    weights[,i,j] <- colMeans(klicOutput$weights)
    
    # Compute ARI
    ari_all[i,j] <- adjustedRandIndex(klic_labels, cluster_labels) 
  }
}

save(ari_one, ari_all, weights, file = "ari-c.RData")