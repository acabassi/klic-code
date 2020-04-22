#################################
####### Synthetic examples ######
### Different levels of noise ###
#################################

rm(list = ls())

library(MASS)
library(ggplot2)
library(reshape)
library(reshape2)

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

all_weights <- array(NA, c(3, 4, n_experiments))
all_ari_all <- all_ari_coca <- all_ari_one <- matrix(NA, 4, n_experiments)

# Load results
for(j in 1:n_experiments){
  load(paste0("../results/ari-b-", j,".RData"))
  all_weights[,,j] <- weights
  all_ari_all[,j] <- ari_all
  all_ari_coca[,j] <- ari_coca
  all_ari_one[,j] <- ari_one
}

ari <- cbind(t(all_ari_one), t(all_ari_all))
colnames(ari) <- c("0", "1", "2", "3", "0+1+2", "0+1+3", "0+2+3", "1+2+3")

ari.m <- melt(ari)
head(ari.m)
colnames(ari.m) <- c("Experiment", "Dataset", "ARI")
ari.m$Dataset <- factor(ari.m$Dataset,
                         levels = c("0","1", "2", "3", "4", "0+1+2", "0+1+3", "0+2+3", "1+2+3"), ordered = TRUE)

ggplot(data = ari.m, aes(x=Dataset, y=ARI)) + geom_boxplot() + ylim(0,1) + my_theme
ggsave("../figures/ari-b.pdf", width = 15, height = 10, units = "cm")

dimnames(all_weights) <- list(c("1st", "2nd", "3rd"), c("0+1+2", "0+1+3", "0+2+3", "1+2+3"), 1:n_experiments)

labels <- matrix(data = c("0","1","2","0","1","3","0","2","3","1","2","3"), nrow = 4, ncol = 3, byrow = TRUE)
weights.m <- melt(all_weights)
head(weights.m)
colnames(weights.m) <- c("Dataset", "Combination", "Experiment", "Weight")
ggplot(data = weights.m, aes(x=Dataset, y=Weight)) + geom_boxplot() + ylim(0,1) + my_theme + facet_grid(cols = vars(Combination))
ggsave("../figures/weights-b.pdf", width = 15, height = 11, units = "cm")

# Plot comparison
ari_comparison <- array(c(all_ari_all, all_ari_coca), dim = c(4, n_experiments, 2))
dimnames(ari_comparison) <- list(c("0+1+2", "0+1+3", "0+2+3", "1+2+3"), 1:n_experiments, c("KLIC", "COCA"))
ari_comparison.m <- melt(ari_comparison)

head(ari_comparison.m)
colnames(ari_comparison.m) <- c("Combination", "Experiment", "Method", "ARI")

ggplot(data = ari_comparison.m, aes(x=Method, y=ARI)) + geom_boxplot() + ylim(0,1) + my_theme + facet_grid(cols=vars(Combination))
ggsave("../figures/coca-b.pdf", width = 15, height = 10, units = "cm")                                  
