# Phenomic Diversity
#population structure analyses reproducing the results of the widely-used computer program structure

#Install packages
install.packages(c("fields","RColorBrewer","mapplots"))
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
install.packages('adegenet')
install.packages("adegenet", repos="http://cran.rstudio.com/", dependencies=TRUE)
install.packages('ggplot2')
install.packages('stringi')
library(stringi)
library(adegenet)
remove.packages("adegenet")
library(dplyr)
library(ggplot2)
#Parallel computing

library(snow)
library(doSNOW)
library(parallel)
detectCores()
cl<-makeCluster(4,type="SOCK")
registerDoSNOW(cl)



#Read in csv data
setwd('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Genotypic Data stuff')
GD = read.table("GWAS_GD.txt", sep = '\t',header = T)

### Population Structure Analysis ###

#Find optimal number of clusters using iterative kmeans approach
## to try different values of k (interactive) using kmeans
#Make genind object to be used in further analysis
      obj <- df2genind(GD, ploidy=2,sep = '/t')
      grp <- find.clusters(obj, max.n=20, n.pca=200, scale=FALSE)
      ## number of accessions per group
      table(grp$grp)
      
      grouping = data.frame(GD$name,grp$grp)
      colnames(grouping)[1] = 'name'
      colnames(grouping)[2] = 'subpop'
      
      #Write out grouping of genotype
      write.csv(grouping, "Population_Clustering.csv",row.names = F)
      grouping <- read.csv("Population_Clustering.csv")
      metadata <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/2017 Data - Growth Chamber/Randomizations Origin Data GWAS Names/Meta_data.csv")
      
#Read in BLUP file with explanatory information
      #setwd('/Volumes/GoogleDrive/My Drive/Phenomic_Diversity/Results')
      dat2 = df<-read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/GWAS/KGF_AdjustedBLUPsAllDays_thinned_Oct19.csv")
      table(dat2$Country)
#Merge grouping and explanatory info
      dat3 = merge(metadata, dat2, by = 'Entry', sort = F)
      
#Modify dat3 so only keep; China, Japan, U.S., and make all other 'other'
      dat3$Country[dat3$Country == 'Algeria'] <- 'Other';dat3$Country[dat3$Country == 'France'] <- 'Other'
      dat3$Country[dat3$Country == 'Georgia'] <- 'Other';dat3$Country[dat3$Country == 'Korea'] <- 'Other'
      dat3$Country[dat3$Country == 'Morocco'] <- 'Other';dat3$Country[dat3$Country == 'Poland'] <- 'Other'
      dat3$Country[dat3$Country == 'Portugal'] <- 'Other';dat3$Country[dat3$Country == 'Taiwan'] <- 'Other'
      dat3$Country[dat3$Country == 'Turkey'] <- 'Other';dat3$Country[dat3$Country == 'Uzbekistan'] <- 'Other'
      dat3$Country[dat3$Country == 'Vietnam'] <- 'Other';dat3$Country[dat3$Country == 'Yugoslavia'] <- 'Other'
      #dat3$Country[dat3$Country == 'Russia'] <- 'Other';dat3$Country[dat3$Country == 'South Korea'] <- 'Other'
      #dat3$Country[dat3$Country == 'North Korea'] <- 'Other'
      table(dat3$Country)

            
### Compute LD, Fst, and genetic distance parameters using NAM package  ###
      install.packages('NAM')
      library(NAM)
      library(phylogram)      
      #Convert GD into matrix form 
      geno = as.matrix(GD[,5:35581])
      rownames(geno) = GD[,1]
      
      #Compute Fst among subpopulations found using the Population Structure Analysis Methods
      #Make fam vector
      fam = as.vector(grouping[,2])
      fst = Fst(geno,fam)
      plot(fst$fst)
      
      #Compute genetic distance using Nei Distance
      gdist = Gdist(geno, method = 1)
      
      fit<-hclust(gdist,method='ward.D')
      plot(fit)
      plot(as.phylo(fit),cex = 0.5,show.tip.label = T)
      tree = cutree(fit, h = 2.18)
      table(tree)
      library(ape)
      plot(as.phylo(fit),type='fan',label.offset=0.1,no.margin=TRUE,cex=0.5,show.tip.label = T )
      
write.csv(pca1_loading,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/PCA.csv", row.names = T)
      
###   PCA of Genotypic Information    ###
      
# PCA with function prcomp
      #pca1 = prcomp(GD[,2:ncol(GD)], scale. = TRUE)
      pca1 <- prcomp(geno, scale. = TRUE)
      # loadings
      pca1_loading = pca1$x
      # add cluster info to pca
      ########################################################################################this doesn't work because of PI names
      #### change underscores in PI names
      str(pca1_loading)
      str(dat3)
      pca1_loading <- cbind(GD[,1:4], pca1_loading)
      names(pca1_loading)[1] <- "Name"
      pca1_loading <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/PCA.csv")
      dat4 = merge(dat3[,1:5], pca1_loading, by = 'Entry', sort = T)
      
      #percent variance explaine
      summary(pca1)
PCA <- as.matrix(pca1)

write.csv(pca1,"C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/PCA.csv", row.names = T)
      
# make figure colored by pca
setwd('C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#2/Genetic Distance/')
tiff("02_Subpopulation_pca.tiff", width = 4, height = 3, units = 'in', res = 300)
ggplot(dat4, aes(x = PC1, y = PC2, color = GenetBackground)) +
  geom_point(alpha=0.5) +
  labs(x = "PC1 (11.4%)", y = "PC2 (6.2%)")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dev.off()

#make figure of subpopuation assignment by country
  dat5 = data.frame(table(dat3$Country,dat3$Subpop))
  tiff("02_Subpopulation_composition.tiff", width = 6, height = 6, units = 'in', res = 300)
  ggplot(dat5, aes(Var2))+
    geom_bar(aes(weight=Freq,fill=Var1), width = 0.5) + 
    theme(axis.text.x = element_text(angle=65, vjust=0.6))+
    theme_classic()+
    scale_color_brewer(palette="Set1")+
    theme(axis.text=element_text(size=8),
          axis.title=element_text(size=8,face="bold"))+
    coord_flip()+
    labs(x = "Cluster", y = 'Count ')+
    guides(fill=guide_legend(title="Country"))
  dev.off()

######## Not In Use ##########
# Not in use:  Alternative method using snapclust implements fast maximum-likelihood (ML) genetic clustering 
#Find optimal number of clusters
#x = snapclust.choose.k(20,obj, IC = BIC, IC.only = T)
#plot(1:20, x, xlab = "Number of clusters (k)",
# ylab = "AIC", type = "b", pch = 20, cex = 3)

## run EM algo with defined number of clusters
res <- snapclust(obj, 5, pop.ini = grp$grp ,hybrids = F)
# names(res)
#res$converged
#res$n.iter
## plot result
#compoplot(res)

#d.dapc <- dapc(obj, n.pca = 20, n.da = 2)
#scatter(d.dapc, clab = 0.85, col = funky(24),
#posi.da="topleft", posi.pca = "bottomleft", scree.pca = TRUE)
  
  
  
  ### Compute LD using synbreed package###
  #Import GM file
  #install.packages("synbreed",repos="http://r-forge.r-project.org")
  library(synbreed)
  #Preprocessing
  #GM = read.table("GWAS_GM.txt", sep = '\t',header = T)
  #colnames(GM)[2] = 'chr'
  #map = GM[,2:3]
  #row.names(map) = GM[,1]
  #geno = as.matrix(GD[,2:35581])
  #rownames(geno) = GD[,1]
  #Convert to gpData format
  #gpData<- create.gpData(pheno=NULL, geno=geno, map=map, covar= NULL, reorderMap=F,map.unit="bp")
  #gData_coded <- codeGeno(gpData, maf= 0.05, nmiss= 0.1, impute=F,
                          #impute.type = "random",  label.heter="alleleCoding",cores = 4) 
  
  #Compute pairwaise LD
  #pLD = pairwiseLD(gData_coded, chr = NULL, type = "matrix",use.plink=FALSE,
                   #ld.threshold=0, ld.window=99999, rm.unmapped = TRUE, cores=4)
  
  #plot(pLD)
  

