#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

rm(list = ls())

library(ggplot2)
library(reshape)
library(reshape2)

# Define ggplot2 theme 
my_basic_theme <-  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.major.x = element_blank() ,
    panel.grid.major.y = element_line(size=.1, color="black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
)

my_theme_rotated_labels <- theme(
  panel.background = element_rect(fill = NA),
  panel.grid.major = element_line(colour = "grey50"),
  panel.grid.major.x = element_blank() ,
  panel.grid.major.y = element_line(size=.1, color="black"),
  axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1)
)

### Data generation ### 
### B. Six clusters, each dataset has a different level of cluster ###
### separability ###

n_experiments <- 100

all_weights <- array(NA, c(3, 4, n_experiments))
all_ari_all <- all_ari_coca <- all_ari_one <- all_ari_one_rbfk <-
  all_ari_one_rbfk_fixed <- all_ari_icluster <- all_ari_clusternomics <-
  all_ari_all_rbfk <- all_ari_all_rbfk_fixed <- matrix(NA, 4, n_experiments)

# Load results
for(j in 1:n_experiments){
  load(paste0("../results/ari-b-", j,".RData"))
  if(!is.na(weights)){
    all_weights[,,j] <- weights
  }
  all_ari_all[,j] <- ari_all
  all_ari_coca[,j] <- ari_coca
  all_ari_one[,j] <- ari_one
  all_ari_icluster[,j] <- ari_icluster
  all_ari_clusternomics[,j] <- ari_clusternomics
  all_ari_all_rbfk[,j] <- ari_all_rbfk
  all_ari_all_rbfk_fixed[,j] <- ari_all_rbfk_fixed
  all_ari_one_rbfk[,j] <- ari_one_rbfk
  all_ari_one_rbfk_fixed[,j] <- ari_one_rbfk_fixed
}

### Plot ARI of KLIC: individual + groups ###
ari <- cbind(t(all_ari_one), t(all_ari_all))
colnames(ari) <- c("0", "1", "2", "3", "0+1+2", "0+1+3", "0+2+3", "1+2+3")

ari.m <- melt(ari)
head(ari.m)
colnames(ari.m) <- c("Experiment", "Dataset", "ARI")
ari.m$Dataset <- factor(ari.m$Dataset,
                         levels = c("0","1", "2", "3", "4", "0+1+2", "0+1+3",
                                    "0+2+3", "1+2+3"),
                        ordered = TRUE)

ggplot(data = ari.m, aes(x=Dataset, y=ARI)) + geom_boxplot(outlier.size = 0.3) +
  ylim(0,1) +
  my_basic_theme
# ggsave("../figures/ari-b.jpg", device = "jpeg", width = 10.5, height = 8,
#        units = "cm")

### Same thing but with RBF kernel ###
ari <- cbind(t(all_ari_one_rbfk), t(all_ari_all_rbfk))
colnames(ari) <- c("0", "1", "2", "3", "0+1+2", "0+1+3", "0+2+3", "1+2+3")
ari.m <- melt(ari)
head(ari.m)
colnames(ari.m) <- c("Experiment", "Dataset", "ARI")
ari.m$Dataset <- factor(ari.m$Dataset,
                        levels = c("0","1", "2", "3", "4", "0+1+2", "0+1+3",
                                   "0+2+3", "1+2+3"),
                        ordered = TRUE)

ggplot(data = ari.m, aes(x=Dataset, y=ARI)) + geom_boxplot(outlier.size = 0.3) +
  ylim(0,1) + my_basic_theme
# ggsave("../figures/ari-b-rbf.jpg", device = "jpeg", width = 15, height = 10,
#        units = "cm")

dimnames(all_weights) <- list(c("1st", "2nd", "3rd"), c("0+1+2", "0+1+3",
                                                        "0+2+3", "1+2+3"),
                              1:n_experiments)

labels <- matrix(data = c("0","1","2","0","1","3","0","2","3","1","2","3"),
                 nrow = 4,
                 ncol = 3,
                 byrow = TRUE)
weights.m <- melt(all_weights)
head(weights.m)
colnames(weights.m) <- c("Dataset", "Combination", "Experiment", "Weight")
ggplot(data = weights.m, aes(x=Dataset, y=Weight)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_basic_theme +
  facet_grid(cols = vars(Combination))
# ggsave("../figures/weights-b.jpg", device = "jpeg", width = 10.5, height = 8,
#        units = "cm")

# Plot comparison
ari_comparison <- array(c(all_ari_all,
                          all_ari_coca,
                          all_ari_all_rbfk,
                          all_ari_all_rbfk_fixed,
                          all_ari_icluster,
                          all_ari_clusternomics),
                        dim = c(4, n_experiments, 6))
dimnames(ari_comparison) <- list(c("0+1+2", "0+1+3", "0+2+3", "1+2+3"),
                                 1:n_experiments,
                                 c("KLIC", "COCA", "RBF opt.", "RBF fixed",
                                   "iCluster", "Clusternomics"))
ari_comparison.m <- melt(ari_comparison)

head(ari_comparison.m)
colnames(ari_comparison.m) <- c("Combination", "Experiment", "Method", "ARI")

ggplot(data = ari_comparison.m, aes(x=Method, y=ARI)) +
  geom_boxplot(outlier.size = 0.35) +
  ylim(0,1) + my_theme_rotated_labels + facet_grid(cols=vars(Combination))
ggsave(paste0("../figures/coca-b-nstart20.jpg"), device = "jpeg",
       width = 14, height = 8,
       units = "cm")                           

