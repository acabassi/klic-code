#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(ggplot2)
library(GGally)
library(reshape)
library(reshape2)

### Define ggplot2 theme ###
my_theme <-  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.major.x = element_blank() ,
    panel.grid.major.y = element_line(size=.1, color="black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
)


n_experiments <- 100
all_ari_one <- all_weights <- matrix(NA, 4, n_experiments)
all_ari_all <- all_ari_coca <- rep(NA, n_experiments)

###  Load results ###
for(j in 1:n_experiments){
  load(paste0("../results/ari-a-", j,".RData"))
  all_ari_one[,j] <- ari_one 
  all_weights[,j] <- weights
  all_ari_all[j] <- ari_all
  all_ari_coca[j] <- ari_coca
}


# Plot ARI

dim(all_ari_one)
dim(as.matrix(all_ari_all))
ari <- cbind(t(all_ari_one), as.matrix(all_ari_all))
colnames(ari) <- c("A", "B", "C", "D", "A+B+C+D")

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")
ari.m$Datasets <- factor(ari.m$Datasets,
                         levels = c("A", "B", "C", "D", "A+B+C+D"),
                         ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot() + ylim(0,1) + my_theme
ggsave("../figures/ari-a.pdf", width = 15, height = 10, units = "cm")

# Plot weights

dim(all_weights)
rownames(all_weights) <- c("A", "B", "C", "D")
weights.m <- melt(t(all_weights))
head(weights.m)
colnames(weights.m) <- c("Experiment", "Dataset", "Weight")
ggplot(data = weights.m, aes(x=Dataset, y=Weight)) + 
  geom_boxplot() + ylim(0,1) + my_theme
ggsave("../figures/weights-a.pdf", width = 12, height = 10, units = "cm")

# Plot comparison

ari_comparison <- cbind(all_ari_all, all_ari_coca)
colnames(ari_comparison) <- c("KLIC", "COCA")
ari_comparison.m <- melt(ari_comparison)
head(ari_comparison.m)
colnames(ari_comparison.m) <- c("Experiment", "Method", "ARI")

ggplot(data = ari_comparison.m, aes(x=Method, y=ARI)) +
  geom_boxplot() + ylim(0,1) + my_theme
ggsave("../figures/coca-a.pdf", width = 6, height = 10, units = "cm")
