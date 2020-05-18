#################################
####### Synthetic examples ######
######## Similar datasets #######
#################################

rm(list = ls())

library(ComplexHeatmap)
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
all_ari_one <- all_ari_one_sparse <- all_ari_one_binary <-
  matrix(NA, 4, n_experiments)
all_ari_all <- all_ari_coca <- all_ari_all_sparse <- all_ari_all_binary <-
  all_ari_coca_sparse <- 
  rep(NA, n_experiments)

separation_level <- 10 # Must be an integer between 1 and 10
# In the main paper, we are showing separation_level = 4

###  Load results ###
for(j in 1:n_experiments){
  load(paste0("../results/ari-a-test-", j,"-sep-", separation_level,".RData"))
  all_ari_one[,j] <- ari_one
  all_ari_one_sparse[,j] <- ari_one_sparse
  all_ari_one_binary[,j] <- ari_one_binary
  all_ari_all[j] <- ari_all
  all_ari_coca[j] <- ari_coca
  all_ari_coca_sparse[j] <- ari_coca_sparcl
  all_ari_all_binary[j] <- ari_all_binary
  all_ari_all_sparse[j] <- ari_all_sparse
}

### Plot weighted kernel matrices of one of the experiments ###

# Heatmap(weighted_kernel)
# Heatmap(weighted_kernel_binary)
# Heatmap(weighted_kernel_sparse)
# Heatmap(moc, cluster_rows = FALSE, cluster_columns = FALSE)
# Heatmap(moc_sparcl, cluster_rows = FALSE, cluster_columns = FALSE)

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
ggsave(paste0("../figures/ari-a-test-sep", separation_level,".jpg"),
       device = "jpeg", width = 15, height = 10,
       units = "cm")

### Plot ARI of KLIC with sparse k-means ###
### 
ari <- cbind(t(all_ari_one_sparse), as.matrix(all_ari_all_sparse))
colnames(ari) <- c("A", "B", "C", "D", "A+B+C+D")

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")
ari.m$Datasets <- factor(ari.m$Datasets,
                         levels = c("A", "B", "C", "D", "A+B+C+D"),
                         ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_basic_theme
ggsave(paste0("../figures/ari-a-test-sep", separation_level,"-sparse.jpg"),
       device = "jpeg", width = 15, height = 10,
       units = "cm")

### Plot ARI of KLIC with sparse k-means ###
### 
ari <- cbind(t(all_ari_one_binary), as.matrix(all_ari_all_binary))
colnames(ari) <- c("A", "B", "C", "D", "A+B+C+D")

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")
ari.m$Datasets <- factor(ari.m$Datasets,
                         levels = c("A", "B", "C", "D", "A+B+C+D"),
                         ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_basic_theme
ggsave(paste0("../figures/ari-a-test-sep", separation_level,"-binary.jpg"),
       device = "jpeg", width = 15, height = 10,
       units = "cm")

# Plot comparison

ari_comparison <- cbind(all_ari_all, all_ari_all_binary, all_ari_all_sparse,
                        all_ari_coca, all_ari_coca_sparse)
colnames(ari_comparison) <- c("KLIC", "KLIC binary matrix", "KLIC sparse k-means",
                              "COCA", "COCA sparse k-means")
ari_comparison.m <- melt(ari_comparison)
head(ari_comparison.m)
colnames(ari_comparison.m) <- c("Experiment", "Method", "ARI")

ggplot(data = ari_comparison.m, aes(x=Method, y=ARI)) +
  geom_boxplot(outlier.size = 0.3) + ylim(0,1) + my_theme_rotated_labels
ggsave(paste0("../figures/coca-a-test-sep",separation_level,".jpg"),
       device = "jpeg", width = 6, height = 10,
       units = "cm")
