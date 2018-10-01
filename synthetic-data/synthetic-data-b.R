#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

library(MASS)
library(ggplot2)
library(GGally)
library(reshape)

### Six clusters, each dataset has a different level of cluster separability

n_experiments <- 100
n_separation_levels <- 6
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

uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

dataset_example <- data[,,6,1]
example6_dataframe <- as.data.frame(dataset_example)
example6_dataframe$cluster <- as.factor(cluster_labels)
names(example6_dataframe) <- list("Var 1", "Var 2", "Cluster")
ggpairs(data=example6_dataframe, # data.frame with variables
        columns=c(1,2,3), # columns to plot, default to all.
        title="Well separated clusters",  # title of the plot
        mapping=ggplot2::aes(colour = Cluster))
ggsave("synthetic-data-b1.pdf", width = 7, height = 7)

dataset_example <- data[,,2,1]
example2_dataframe <- as.data.frame(dataset_example) 
example2_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example2_dataframe, # data.frame with variables
        columns=c(1,2,3), # columns to plot, default to all.
        title="Not so well separated clusters", # title of the plot 
        mapping = ggplot2::aes(colour= cluster))
ggsave("synthetic-data-b2.pdf", width = 7, height = 7)

save(data, file = "synthetic-data-b.RData")
