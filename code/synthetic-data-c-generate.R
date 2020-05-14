#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(MASS)
library(GGally)

# Define ggplot2 theme 
my_theme <-  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.major.x = element_blank() ,
    panel.grid.major.y = element_line(size=.1, color="black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
)

### Data generation ### 
### C1. Six clusters ###

n_experiments <- 100
# n_separation_levels <- 1
n_types <- 3
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6

N <- n_obs_per_cluster*n_clusters
P <- n_variables
data <- array(NA, c(N, P, n_types, n_experiments))

Sigma <- diag(n_variables)

i = 3 # Medium separation level
# (Change value to 8 just for the illustrative plot in the paper)
for(j in 1:n_experiments){
  for(t in 1:(n_types-1)){
    for(k in 1:(n_clusters/t)){
      mu = rep(k*(i-1)/3, n_variables)
      data[((k-1)*n_obs_per_cluster*t+1):(k*n_obs_per_cluster*t),,t,j] <-
        mvrnorm(n = n_obs_per_cluster*t, mu, Sigma)
    }
  }
}

i = 9 # Higher separation level
t = 1 # Only for dataset containing 6 clusters
for(j in 1:n_experiments){
  mu = rep(NA, N)
    for(k in 1:(n_clusters/t)){
      mu = rep(k*(i-1)/3, n_variables)
      data[((k-1)*n_obs_per_cluster*t+1):(k*n_obs_per_cluster*t),,3,j] <-
        mvrnorm(n = n_obs_per_cluster*t, mu, Sigma)
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
        mapping=ggplot2::aes(colour = cluster),
        lower=list(combo=wrap("facethist", binwidth=0.8)))
ggsave("../figures/synthetic-data-c1.jpg",
       device = "jpeg", width = 7, height = 7)

t = 2 # 3 clusters
dataset_example <- data[,,t,1]
example3_dataframe <- as.data.frame(dataset_example)
example3_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example3_dataframe, # data.frame with variables
        title="3 clusters", # title of the plot
        mapping = ggplot2::aes(colour= cluster),
        lower=list(combo=wrap("facethist", binwidth=0.8)))
ggsave("../figures/synthetic-data-c2.jpg",
       device = "jpeg", width = 7, height = 7)

t = 3 # 6 clusters, well separated
dataset_example <- data[,,t,1]
example6plus_dataframe <- as.data.frame(dataset_example)
example6plus_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example6plus_dataframe, # data.frame with variables
        title="6 clusters, well separated", # title of the plot
        mapping = ggplot2::aes(colour= cluster),
        lower=list(combo=wrap("facethist", binwidth=0.8)))
ggsave("../figures/synthetic-data-c3.jpg",
       device = "jpeg", width = 7, height = 7)

save(data, file = "../data/synthetic-data-c.RData")
