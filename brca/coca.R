##############################################
### Breast cancer data analysis using COCA ###
##############################################

rm(list = ls())

library(devtools)
install_github("acabassi/coca")
library(coca)
library(cluster)

load("Breast.RData")

maxK <- 10
output <- list()
KM <- list()
coph <- list()
silh <- matrix(NA, 4, 9)
MOCmatrix <- matrix(NA, 0, 348)

# For each dataset
for(i in 1:4){
  
  # Transpose and scale datasets
  dataset_i <- scale(t(X[[i]]))
  
  # For all possible number of clusters
  for(j in 2:10){
    
    # Use k-means to find clustering 
    kmeans_j <- kmeans(dataset_i, j, iter.max = 1000, nstart = 20)
    # Compute silhouette 
    silh_ij <- silhouette(kmeans_j$cluster, dist = dist(dataset_i, method = "euclidean"))
    # Save average silhouette
    silh[i,j-1] <- mean(silh_ij[,3])
  }
  
  # Choose number of clusters that maximises the silhouette
  k_i <- which.max(silh[i,]) + 1
  # Re-do k-means with that number of clusters
  kmeans_k_i <- kmeans(dataset_i, k_i, iter.max = 1000, nstart = 20)
  
  # Add corresponding lines to the MOC matrix
  for(j in 1:k_i){
    MOCmatrix <- rbind(MOCmatrix, (kmeans_k_i$cluster==j)*1)
  }
  prod(colSums(MOCmatrix)) # just to check...
  
}

silh_final <- rep(0, 9)

# For each possible number of clusters
for(j in 2:10){
  
  # Find clustering using COCA
  COCAclusters_i <- coca(t(MOCmatrix), 2, B = 1000, pItem = 0.8, hclustMethod = 'average')  
  # Compute silhouette
  silh_final_j <- silhouette(kmeans_j$cluster, dist = dist(t(MOCmatrix), method = "euclidean"))
  # Save average silhouette
  silh_final[j-1] <- mean(silh_final_j[,3])
}

# Choose number of clusters that maximises the silhouette
k_final = which.max(silh_final) + 1
# Re-do COCA with chosen number of clusters
COCAclusters <- coca(t(MOCmatrix), k_final, B = 1000, pItem = 0.8, hclustMethod = 'average')  

# Save output (consensus matrix + cluster labels)
save(COCAclusters, file = "coca.RData")
