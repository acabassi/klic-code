#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

rm(list = ls())

library(coca)
library(ComplexHeatmap)
library(klic)
library(ggplot2)
library(mclust) # for adjustedRandIndex
library(rdetools) # for rbfkernel
library(reshape)

### B. Six clusters, each dataset has a different level of cluster ###
###    separability ###

n_experiments <- 100
n_separation_levels <- 10
n_variables <- 2
n_obs_per_cluster <- 50
n_clusters <- 6
n_sigma_values <- 10

N <- n_obs_per_cluster*n_clusters
P <- n_variables

# Sigma values
sigmas <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50)
              
# True cluster labels 
uno <- rep(1, n_obs_per_cluster)
cluster_labels <- c(uno, uno*2, uno*3, uno*4, uno*5, uno*6)

new_data <- array(NA, c(N, P, n_separation_levels, n_experiments))
for(i in 1:n_separation_levels){
  load(paste0("../data/synthetic-data-a-sep", i,".RData"))
  new_data[,,i,] <- data[,,1,]
}

ari <- array(NA, c(n_experiments, n_separation_levels, n_sigma_values))
dimnames(ari) <- list(as.character(1:n_experiments),
                      c("s = 0", "s = 0.5", "s = 1", "s = 1.5", "s = 2",
                        "s = 2.5", "s = 3", "s = 3.5", "s = 4", "s = 4.5"),
                      as.character(sigmas))
parameters <- list()
parameters$cluster_count <- n_clusters
parameters$iteration_count <- 100

### Choice of the parameter sigma for the RBF kernels ###

tryClustering <- function(data, sigmas, i, j, k){
  CM_rbfk <- rbfkernel(data[,,i,j], sigma = sigmas[k])
  CM_rbfk[lower.tri(CM_rbfk)] <- t(CM_rbfk)[lower.tri(CM_rbfk)]
  CM_rbfk <- spectrumShift(CM_rbfk)
  clusters <- kkmeans(CM_rbfk, parameters)$clustering
  ari <- adjustedRandIndex(clusters, cluster_labels)
}


for(i in 1:n_separation_levels){
  for (j in 1:n_experiments){
    for(k in 1:n_sigma_values){
    ari[j, i, k] <- tryCatch(tryClustering(new_data, sigmas, i, j, k),
                             error = function(err) NA)
    }
  }
}

# Define ggplot2 theme
my_theme_rotated_labels <- theme(
  panel.background = element_rect(fill = NA),
  panel.grid.major = element_line(colour = "grey50"),
  panel.grid.major.x = element_blank() ,
  panel.grid.major.y = element_line(size=.1, color="black"),
  axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1)
)

ari.m <- melt(ari)
head(ari.m)
colnames(ari.m) <- c("Experiment", "Separation", "Sigma", "ARI")
ari.m$Sigma <- as.factor(ari.m$Sigma)
ggplot(data=ari.m, aes(x=Sigma, y=ARI)) + geom_boxplot(outlier.size = 0.35) +
  ylim(0, 1) + my_theme_rotated_labels + facet_wrap(~ Separation, nrow = 2)
ggsave("../figures/choiceRBFsigma.jpg", device = "jpeg",
       width = 18, height = 20,
       units = "cm")

# Choose most appropriate value of sigma for each separation level
RBF_sigmas <- 
  sigmas[apply(apply(ari, MARGIN = c(2,3), FUN = mean, na.rm = TRUE),
               MARGIN = 1,
               FUN = which.max)]
save(RBF_sigmas, file = "../data/RBF_sigma_values.RData")

# This gives 
# Sep 0: 1e-03
# Sep 0.5: 1e-01
# Sep 1: 1e-01
# Sep 1.5: 5e-01
# Sep 2: 1e+00
# Sep 2.5: 5e+00
# Sep 3: 5e+00
# Sep 3.5: 5e+00
# Sep 4: 1e+01
# Sep 4.5: 1e+01
