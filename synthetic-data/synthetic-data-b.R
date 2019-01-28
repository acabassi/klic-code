#################################
####### Synthetic examples ######
### Different levels of noise ###
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

### Data generation ### 
### B. Six clusters, each dataset has a different level of cluster separability ###

n_experiments <- 100
n_separation_levels <- 4
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

Sigma <- diag(n_variables)

N <- n_obs_per_cluster*n_clusters
P <- n_variables
data <- array(NA, c(N, P, n_separation_levels, n_experiments))

set.seed(151)

for(i in 1:n_separation_levels){
  for(j in 1:n_experiments){
    mu = rep(NA, N)
    for(k in 1:n_clusters){
      mu = rep(k*(i-1)/2, n_variables)
      data[((k-1)*n_obs_per_cluster+1):(k*n_obs_per_cluster),,i,j] <- mvrnorm(n = n_obs_per_cluster, mu, Sigma)
    }
  }
}

# True cluster labels 
uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

### Plots ###

# Well separated clusters
dataset_example <- data[,,4,1]
example6_dataframe <- as.data.frame(dataset_example)
example6_dataframe$cluster <- as.factor(cluster_labels)
names(example6_dataframe) <- list("Var 1", "Var 2", "Cluster")
ggpairs(data=example6_dataframe, # data.frame with variables
        title="Well separated clusters",  # title of the plot
        mapping=ggplot2::aes(colour = Cluster))
ggsave("synthetic-data-b1.pdf", width = 7, height = 7)

# Not so well separated clusters
dataset_example <- data[,,2,1]
example2_dataframe <- as.data.frame(dataset_example) 
example2_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example2_dataframe, # data.frame with variables
        title="Not so well separated clusters", # title of the plot 
        mapping = ggplot2::aes(colour= cluster))
ggsave("synthetic-data-b2.pdf", width = 7, height = 7)

save(data, file = "synthetic-data-b.RData")

### Clustering one dataset at a time ###

ari_one <- array(NA, c(n_separation_levels, n_experiments))

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters

for(i in 1:n_separation_levels){
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

### Combining subsets of three datasets ###

n_subsets <- 4
n_datasets_per_subset <- 3

ari_all <- array(NA, c(n_subsets, n_experiments))
weights <- array(NA, c(n_datasets_per_subset, n_subsets, n_experiments))

subsets <- rbind(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4))

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

### COCA ###

ari_coca <- rep(NA, n_experiments)
moc <- array(NA, c(dim(data)[1], n_clusters*n_datasets_per_subset))

for(i in 1:n_experiments){
  for(j in 1:n_subsets){
    
  # Select datasets 
  datasets_in_subset <- subsets[i,]
  
  # Over-write the matrix-of-clusters each time
  count <- 0
  
  # Find clusters in each dataset
  for(l in datasets_in_subset){
    
    kmeans_cluster_labels <- kmeans(data[,,l,i], n_clusters)$cluster
    
    # Fill the matrix-of-clusters
    for(m in 1:n_clusters){
      count <- count + 1
      moc[, count] <- (kmeans_cluster_labels == m)*1
    }
  }
  
  # Use COCA to find final clusters
  coca_cluster_labels <- coca(moc, n_clusters)$clusterLabels
  # Compute ARI 
  ari_coca[j,i] <- adjustedRandIndex(cluster_labels, coca_cluster_labels)
  }
}

# Save results
save(ari_one, ari_all, weights, ari_coca, file = "ari-b.RData")

# Load results
load("ari-b.RData")

dim(ari_one)
dim(as.matrix(ari_all))
ari <- cbind(t(ari_one), t(ari_all))
colnames(ari) <- c("0", "1", "2", "3", "0+1+2", "0+1+3", "0+2+3", "1+2+3")

ari.m <- melt(ari)
ari.m # pasting some rows of the melted data.frame

ggplot(data = ari.m, aes(x=X2, y=value)) + geom_boxplot() + ylim(0,1)
ggsave("ari-b.pdf")

dim(weights)

dimnames(weights) <- list(c("1st", "2nd", "3rd"), c("0+1+2", "0+1+3", "0+2+3", "1+2+3"), 1:100)

labels <- matrix(data = c("0","1","2","0","1","3","0","2","3","1","2","3"), nrow = 4, ncol = 3, byrow = TRUE)
for(i in 1:4){
  weights_i <- weights[,i,]
  rownames(weights_i) <- labels[i,]
  weights.m <- melt(t(weights_i))
  ggplot(data = weights.m, aes(x=as.character(X2), y=value)) + geom_boxplot() + ylim(0,1)
  ggsave(paste("weights-b",as.character(i),".pdf", sep=""))
}

