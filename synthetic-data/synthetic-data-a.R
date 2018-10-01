#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(MASS)
library(ggplot2)
library(GGally)
library(reshape)

### Six clusters, same cluster separability in each dataset

n_experiments <- 100
n_separation_levels <- 1
n_datasets_same_rho <- 5
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

Sigma <- diag(n_variables)

N <- n_obs_per_cluster*n_clusters
P <- n_variables
data <- array(NA, c(N, P, n_datasets_same_rho, n_experiments))

separation_level <- 4 # Same for all datasets

set.seed(151)

for(i in 1:n_datasets_same_rho){
  for(j in 1:n_experiments){
    mu = rep(NA, N)
    for(k in 1:n_clusters){
      mu = rep((k*(separation_level-1)/2), n_variables)
      data[((k-1)*n_obs_per_cluster+1):(k*n_obs_per_cluster),,i,j] <- mvrnorm(n = n_obs_per_cluster, mu, Sigma)
    }
  }
}

uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

dataset_example <- data[,,1,1]
example6_dataframe <- as.data.frame(dataset_example)
example6_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6_dataframe, # data.frame with variables
        columns=c(1,2,3), # columns to plot, default to all.
        title="Datasets with the same cluster separability",  # title of the plot
        mapping=ggplot2::aes(colour = cluster))
ggsave("synthetic-data-a.pdf", width = 7, height = 7)

save(data, file = "synthetic-data-a.RData")
