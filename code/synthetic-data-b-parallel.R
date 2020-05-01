#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

rm(list = ls())

library(coca)
library(ComplexHeatmap)
library(iCluster)
library(klic)
library(mclust) # for adjustedRandIndex
library(rdetools) # for rbfkernel

j <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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

load("../data/synthetic-data-b.RData")

### Clustering one dataset at a time ###

ari_one <- ari_one_rbfk <-rep(NA, n_separation_levels)

# Initialise parameters for kernel k-means
parameters <- list()
parameters$cluster_count <- n_clusters
parameters$iteration_count <- 100

# Set parameters for RBF kernels
RBFsigma <- c(0.001, 0.100, 1.000, 5.000)

CM_rbfk <- CM_cc <- array(NA, c(N, N, n_separation_levels))
dimnames(CM_rbfk) <- dimnames(CM_cc) <-
  list(as.character(1:N),
       as.character(1:N),
       as.character(1:n_separation_levels))

for(i in 1:n_separation_levels){
    # Use consensus clustering to find kernel matrix
    CM_cc_temp <- consensusCluster(data[,,i,j], n_clusters)
    # Use RBF kernel to obtain another kernel matrix
    CM_rbfk_temp <- rbfkernel(data[,,i,j], sigma = RBFsigma[i])
    # Make sure that it is actually symmetric
    CM_rbfk_temp[lower.tri(CM_rbfk_temp)] <-
      t(CM_rbfk_temp)[lower.tri(CM_rbfk_temp)]
    # Shift the eigenvalues of the kernel matrices so that they are positive
    # semi-definite
    CM_cc[,,i] <- spectrumShift(CM_cc_temp)
    CM_rbfk[,,i] <- spectrumShift(CM_rbfk_temp)
    # Use kernel k-means to find clusters
    kkmeans_labels <- kkmeans(CM_cc[,,i], parameters)$clustering
    kkmeans_labels_rbfk <- kkmeans(CM_rbfk[,,i], parameters)$clustering
    # Compute ARI
    ari_one[i] <- adjustedRandIndex(kkmeans_labels, cluster_labels)
    ari_one_rbfk[i] <- adjustedRandIndex(kkmeans_labels_rbfk, cluster_labels)
}

### Combining subsets of three datasets ###

n_subsets <- 4
n_datasets_per_subset <- 3

ari_all <- rep(NA, c(n_subsets))
weights <- array(NA, c(n_datasets_per_subset, n_subsets))

subsets <- rbind(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4))

data_for_klic <- list()
weighted_kernel <- list()
  
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
  all_weights <- klicOutput$weights
  weights[,i] <- colMeans(all_weights)
  
  weighted_kernel[[i]]<- matrix(0, N, N) 
  count <- 1
  for (l in datasets_in_subset) {
    weighted_kernel[[i]] <- weighted_kernel[[i]] +
      (all_weights[, count] %*% t(all_weights[, count])) * CM_cc[, , count]
    count <- count + 1
  }

  # Compute ARI
  ari_all[i] <- adjustedRandIndex(klic_labels, cluster_labels)
}

### COCA, iCluster and RBF kernels ###
ari_coca <- ari_icluster <- ari_all_rbfk <- rep(NA, n_subsets)
moc <- array(NA, c(dim(data)[1], n_clusters*n_datasets_per_subset))

weighted_kernel_rbfk <- list()

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

  ### Kernel k-means with RBF kernel ###
  lmkkmeans_rbfk <- lmkkmeans(CM_rbfk[,,datasets_in_subset], parameters)
  rbfk_cluster_labels <- lmkkmeans_rbfk$clustering
  weights_rbfk <- lmkkmeans_rbfk$Theta
  ari_all_rbfk[i] <- adjustedRandIndex(rbfk_cluster_labels, cluster_labels)

  ### Weighted kernels ###
  
  weighted_kernel_rbfk[[i]] <- matrix(0, N, N) 
  count <- 1
  for (l in datasets_in_subset) {
    weighted_kernel_rbfk[[i]] <- weighted_kernel_rbfk[[i]] +
      (weights_rbfk[, count] %*% t(weights_rbfk[, count])) * CM_rbfk[, , count]
    count <- count + 1
  }
}

### Save results ###
save(ari_one, ari_one_rbfk,
     ari_all, ari_coca, ari_icluster, ari_all_rbfk,
     weights, weighted_kernel, weighted_kernel_rbfk,
     file = paste0("../results/ari-b-", j, ".RData"))

### Plot kernels for paper ###
# if(j==1){
  library(ComplexHeatmap)
  library(circlize)
  col_fun = colorRamp2(c(0, 1), c("white","#003C71")) # Dark blue
  label_colors <- c("#6CACE4", # Light blue
                    "#E89CAE", # Light pink
                    "#F1BE48", # Light yellow
                    "#B7BF10", # Light green
                    "#85b09A", # Light cambridge blue
                    "#0072ce") # Core blue
  names(label_colors) <- as.character(1:n_clusters)
  row_annotation <- rowAnnotation(Label = as.character(cluster_labels),
                                 col = list(Label = label_colors),
                                 show_legend = FALSE,
                                 show_annotation_name = FALSE,
                                 annotation_width = unit(0.1, "cm"))

    H1 <- Heatmap(CM_cc[,,1],
          col = col_fun,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_heatmap_legend = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          heatmap_width = unit(5, "cm"),
          heatmap_height = unit(5, "cm"),
          right_annotation = row_annotation)
    H2 <- Heatmap(CM_cc[,,2],
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_heatmap_legend = FALSE,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 heatmap_width = unit(5, "cm"),
                 heatmap_height = unit(5, "cm"),
                 right_annotation = row_annotation)
    H3 <- Heatmap(CM_cc[,,3],
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_heatmap_legend = FALSE,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 heatmap_width = unit(5, "cm"),
                 heatmap_height = unit(5, "cm"),
                 right_annotation = row_annotation)
    H4 <- Heatmap(CM_cc[,,4],
                 col = col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_heatmap_legend = TRUE,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 heatmap_width = unit(5, "cm"),
                 heatmap_height = unit(5, "cm"),
                 right_annotation = row_annotation,
                 heatmap_legend_param = list(title = ""))

jpeg("../figures/heatmap-b.jpg",
     height = 5.5, width = 23, units = "cm", res = 1200)
  H1 + H2 + H3 + H4
dev.off()