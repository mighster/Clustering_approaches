#Script to try to cluster and plot discrete values (like SNP data) for Kevin
setwd("~/Desktop/KevinFalk/Re_Heatmap_fun/")
library(dplyr) # for data cleaning
library(ISLR) # for college dataset
library(cluster) # for gower similarity and pam
library(Rtsne) # for t-SNE plot
library(ggplot2) # for visualization


###based off online examples using 'flower' as the dataset:
#data(flower)
#head(flower) # this gives you a dataframe that is 18 objects of 8 variables
#str(flower)  # this prints out the 8 variables with the different associated facts about each one
### We will want to cluster the SNPs using the Gower distance because Euclidian distances are only applicable to continuous datasets.  
mydata = read.csv("~/Desktop/KevinFalk/Re_Heatmap_fun/SignificantSNPsForKevin2.csv")
mydata = read.csv("~/Desktop/KevinFalk/Re_Heatmap_fun/SignificantSNPsForKevin3a.csv")
head(mydata)

df<-data.frame(mydata)
m <- df[,2:31]
rownames(m)<-make.names(df[,1], unique="True")
head(m)
rownames(m)
m_matrix <- data.matrix(m)
head(m_matrix)

gower_dist <- daisy(m_matrix, metric = "gower", stand = FALSE)
summary(gower_dist)
head(gower_dist)

hac <-agnes(gower_dist, diss = TRUE)
hac$order
plot(hac)
dev.off()

library(CrossClustering)
cross.clust <- CrossClustering(grower_dist, k.w.min = 2, k.w.max = 5,
                              k.c.max = 6, out = TRUE)
