#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

rm(list = ls())

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
### B. Six clusters, each dataset has a different level of cluster ###
### separability ###

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
      mu = rep(k*(i-1), n_variables)
      data[((k-1)*n_obs_per_cluster+1):(k*n_obs_per_cluster),,i,j] <-
        mvrnorm(n = n_obs_per_cluster, mu, Sigma)
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
        mapping=ggplot2::aes(colour = Cluster),
        lower=list(combo=wrap("facethist", binwidth=0.8)))
ggsave("../figures/synthetic-data-b1.pdf", width = 7, height = 7)

# Not so well separated clusters
dataset_example <- data[,,2,1]
example2_dataframe <- as.data.frame(dataset_example) 
example2_dataframe$cluster <- as.factor(cluster_labels)
ggpairs(data=example2_dataframe, # data.frame with variables
        title="Not so well separated clusters", # title of the plot 
        mapping = ggplot2::aes(colour= cluster),
        lower=list(combo=wrap("facethist", binwidth=0.8)))
ggsave("../figures/synthetic-data-b2.pdf", width = 7, height = 7)

save(data, file = "../data/synthetic-data-b.RData")