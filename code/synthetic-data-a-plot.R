#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(ggplot2)
library(GGally)
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

##


n_experiments <- 100
all_ari_one <- all_ari_one_rbfk <- all_weights <- matrix(NA, 4, n_experiments)
all_ari_all <- all_ari_coca <- all_ari_icluster <- all_ari_all_rbfk <-
  rep(NA, n_experiments)

separation_level <- 4 # Must be an integer between 1 and 10
# In the main paper, we are showing separation_level = 4

###  Load results ###
for(j in 1:n_experiments){
  load(paste0("../results/ari-a-", j,"-sep-", separation_level,".RData"))
  all_ari_one[,j] <- ari_one
  all_ari_one_rbfk[,j] <- ari_one_rbfk
  all_weights[,j] <- weights
  all_ari_all[j] <- ari_all
  all_ari_coca[j] <- ari_coca
  all_ari_icluster[j] <- ari_icluster
  all_ari_all_rbfk[j] <- ari_all_rbfk
}


### Plot ARI of KLIC ###

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
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_basic_theme
ggsave(paste0("../figures/ari-a-sep", separation_level,".jpg"),
       device = "jpeg", width = 15, height = 10,
       units = "cm")

### Plot weights of KLIC ###

dim(all_weights)
rownames(all_weights) <- c("A", "B", "C", "D")
weights.m <- melt(t(all_weights))
head(weights.m)
colnames(weights.m) <- c("Experiment", "Dataset", "Weight")
ggplot(data = weights.m, aes(x=Dataset, y=Weight)) + 
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_basic_theme
ggsave(paste0("../figures/weights-a-sep", separation_level, ".jpg"),
       device = "jpeg", width = 12, height = 10,
       units = "cm")

### Plot ARI of RBF ###
### 
ari <- cbind(t(all_ari_one_rbfk), as.matrix(all_ari_all_rbfk))
colnames(ari) <- c("A", "B", "C", "D", "A+B+C+D")

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")
ari.m$Datasets <- factor(ari.m$Datasets,
                         levels = c("A", "B", "C", "D", "A+B+C+D"),
                         ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_basic_theme
ggsave(paste0("../figures/ari-a-sep", separation_level,"-rbf.jpg"),
       device = "jpeg", width = 15, height = 10,
       units = "cm")

# Plot comparison

ari_comparison <- cbind(all_ari_all, all_ari_coca,
                        all_ari_icluster, all_ari_all_rbfk)
colnames(ari_comparison) <- c("KLIC", "COCA", "iCluster", "RBF")
ari_comparison.m <- melt(ari_comparison)
head(ari_comparison.m)
colnames(ari_comparison.m) <- c("Experiment", "Method", "ARI")

ggplot(data = ari_comparison.m, aes(x=Method, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme_rotated_labels
ggsave(paste0("../figures/coca-a-sep",separation_level,".jpg"),
       device = "jpeg", width = 6, height = 10,
       units = "cm")
