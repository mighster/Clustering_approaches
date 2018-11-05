library(dplyr) # for data cleaning
library(ISLR) # for college dataset
library(cluster) # for gower similarity and pam
library(Rtsne) # for t-SNE plot
library(ggplot2) # for visualization
library(WGCNA)
########################################
#https://www.r-bloggers.com/clustering-mixed-data-types-in-r/
#######################################

set.seed(1680) # for reproducibility

data<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/GWAS/KGF_AdjustedBLUPsAllDays_thinned_Oct19.csv")
metadata <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Randomizations Origin Data GWAS Names/Meta_data.csv", sep=",", header=T, check.names = FALSE)

Merged_df <- left_join(metadata, data, by="Entry")

Merged_df$Country=as.factor(Merged_df$Country)
Merged_df$Region=as.factor(Merged_df$Region)
Merged_df$Diversity=as.factor(Merged_df$Diversity)
Merged_df$Diversity=as.factor(Merged_df$Diversity)
Merged_df$`Stem Termination`=as.factor(Merged_df$`Stem Termination`)
Merged_df$`seed coat color`=as.factor(Merged_df$`seed coat color`)
Merged_df$`hilum color`=as.factor(Merged_df$`hilum color`)
Merged_df$MG=as.factor(Merged_df$MG)
Merged_df$Cluster.6=as.factor(Merged_df$Cluster.6)
Merged_df$Cluster.8=as.factor(Merged_df$Cluster.8)
Merged_df$Cluster.9=as.factor(Merged_df$Cluster.9)


Merged_df <- Merged_df[-c(9,98,167,182,190,202,204,273),]
new_df <- Merged_df[,c(3:7,9:13,99:136)]
new_df[1:292,1:2]

gower_dist <- daisy(new_df[, -c(1:3)],
                    metric = "gower",
                    type = list(logratio = 3))

gower_mat <- as.matrix(gower_dist)

new_df[
  which(gower_mat == min(gower_mat[gower_mat != min(gower_mat)]),
        arr.ind = TRUE)[1, ], ]

new_df[
  which(gower_mat == max(gower_mat[gower_mat != max(gower_mat)]),
        arr.ind = TRUE)[1, ], ]

sil_width <- c(NA)

for(i in 2:10){
    pam_fit <- pam(gower_dist,
                 diss = TRUE,
                 k = i)
    sil_width[i] <- pam_fit$silinfo$avg.width
}

plot(1:10, sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, sil_width)

pam_fit <- pam(gower_dist, diss = TRUE, k = 2)

pam_results <- new_df %>%
  mutate(cluster = pam_fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))

pam_results$the_summary

new_df[pam_fit$medoids, ]

tsne_obj <- Rtsne(gower_dist, is_distance = TRUE)

tsne_data_Day12 <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         name = Merged_df$Name)


tiff("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Gower/tSNE_Day12.tiff", units = "in" ,height = 8, width = 10, res=150)
ggplot(aes(x = X, y = Y), data = tsne_data_Day12) +
  geom_point(aes(color = cluster, size = 0.25))+
  ggtitle("tSNE Clustering at 12 Days")
dev.off()

write.csv(tsne_data,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Gower/tsne_data_Day12.csv", row.names = T)

Merged_tsne <- left_join(tsne_data_Day6, tsne_data_Day9,by="name")
Merged_tsne <- left_join(Merged_tsne, tsne_data_Day12,by="name")

write.csv(Merged_tsne,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Gower/tsne_data_mergedDays.csv", row.names = T)
