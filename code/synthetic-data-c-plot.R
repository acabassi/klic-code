#################################
####### Synthetic examples ######
######## Nested clusters ########
#################################

rm(list=ls())

library(MASS)
library(ggplot2)
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

### Load results ###
n_experiments <- 100
all_weights <- array(NA, c(2, 2, n_experiments))
all_ari_all <- array(NA, c(2, n_experiments))
all_ari_one <- array(NA, c(3, n_experiments))

for(j in 1:n_experiments){
  load(paste0("../results/ari-c-", j,".RData"))
  all_weights[,,j] <- weights
  all_ari_all[,j] <- ari_all
  all_ari_one[,j] <- ari_one
}

ari <- cbind(t(all_ari_one), t(all_ari_all))
colnames(ari) <- c("6", "3", "6*", "3+6", "3+6*")

ari.m <- melt(ari)
head(ari.m) # pasting some rows of the melted data.frame
colnames(ari.m) <- c("Experiment", "Datasets", "ARI")

ari.m$Datasets <- factor(ari.m$Datasets,
                       levels = c("3","6", "6*", "3+6", "3+6*"), ordered = TRUE)

ggplot(data = ari.m, aes(x=Datasets, y=ARI)) +
  geom_boxplot() + ylim(0,1) + my_theme
ggsave("../figures/ari-c.pdf", width = 15, height = 10, units = "cm")

dim(weights)
 
weights_1 <- all_weights[,1,]
weights_2 <- all_weights[,2,]

rownames(weights_1) <- c("6", "3")
rownames(weights_2) <- c("3", "6*")
colnames(weights_1) <- colnames(weights_2) <- 1:n_experiments

weights_1.m <- melt(t(weights_1))
head(weights_1.m) 
colnames(weights_1.m) <- c("Experiments", "Dataset", "Weight")
weights_1.m$Dataset <- factor(weights_1.m$Dataset,
                         levels = c("3","6"), ordered = TRUE)


ggplot(data = weights_1.m, aes(x=Dataset, y=Weight)) +
  geom_boxplot() + ylim(0,1) + my_theme
ggsave("../figures/weights-c1.pdf",
       width = 6, height = 10, units = "cm")

weights_2.m <- melt(t(weights_2))
head(weights_2.m) 
colnames(weights_2.m) <- c("Experiments", "Dataset", "Weight")

ggplot(data = weights_2.m, aes(x=Dataset, y=Weight)) +
  geom_boxplot() + ylim(0,1) + my_theme
ggsave("weights-c2.pdf",
       width = 6, height = 10, units = "cm")
