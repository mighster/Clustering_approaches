####Construction of phylogenetic tree
#https://popgen.nescent.org/2015-05-18-Dist-SNP.html
setwd("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff/")
install.packages("ape")
install.packages("ggtree")
library(ape)
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
## biocLite("BiocUpgrade") ## you may need this
biocLite("ggtree")
biocLite("treeio")
library(ggtree)
library(treeio)
library(ggplot2)

#Input SNP File
a<-read.table("KGF_KAP_GWAS_300_50k_SNP_maf_.05_GD_data.txt",header = T)
germans<-read.csv("a.csv",header = F)
orig<-a

country <- read.csv("Genotype Origin Data_Basic.csv", header=F)
head(b[1:5,1:5])

b <- cbind(germans[,3],a)

#colnames(b) <- c
row.names(b)<-b$PI
colnames(b)
dm = dist(b,method="euclidean")
fit<-hclust(dm,method='ward.D')





#Number of cluster
#Elbow Method for finding the optimal number of clusters
#https://www.r-bloggers.com/finding-optimal-number-of-clusters/
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- dm
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

cluster<-kmeans(data,centers = 6, nstart=50,iter.max = 15 )


#Assign cluster number to genotype
clusters<-data.frame(cbind(cluster$cluster,orig$PI))

write.table(a$PI,"a.csv")
write.table(a$PI,"a.csv")

label=a[,1]

#Use base plotting function
color=c("firebrick",'blue','purple',"black","darkgreen","tomato4")
nclus=6
color_list=rep(color,nclus/length(color))
clus=cutree(fit,nclus)

tiff("GWAS_Phylo.tiff", height = 8, width = 8, res = 600,units = "in")
plot(as.phylo(fit),type='fan',label.offset=0.1,no.margin=TRUE,cex=0.5,show.tip.label = T,tip.color =color_list[clus] )
dev.off()






