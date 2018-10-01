#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(MASS)
library(ggplot2)
library(GGally)
library(reshape)

### Primo caso: sei cluster, sei diversi valori di rho

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

t = 1 # 6 clusters
uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

dataset_example <- data[,,t,1]
example6_dataframe <- as.data.frame(dataset_example)
example6_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6_dataframe, # data.frame with variables
        # columns=c(1,4,7), # columns to plot, default to all.
        title="6 clusters",  # title of the plot
        mapping=ggplot2::aes(colour = cluster))
ggsave("synthetic-data-c1.pdf", width = 7, height = 7)

t = 2
dataset_example <- data[,,t,1]
example3_dataframe <- as.data.frame(dataset_example)
example3_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example3_dataframe, # data.frame with variables
        # columns=c(1,2,6), # columns to plot, default to all.
        title="3 clusters", # title of the plot
        mapping = ggplot2::aes(colour= cluster))
ggsave("synthetic-data-c2.pdf", width = 7, height = 7)

t = 3
dataset_example <- data[,,t,1]
example6plus_dataframe <- as.data.frame(dataset_example)
example6plus_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6plus_dataframe, # data.frame with variables
        # columns=c(1,2,6), # columns to plot, default to all.
        title="6 clusters, better separated", # title of the plot
        mapping = ggplot2::aes(colour= cluster))
ggsave("synthetic-data-c3.pdf", width = 7, height = 7)

save(data, file = "synthetic-data-c.RData")
