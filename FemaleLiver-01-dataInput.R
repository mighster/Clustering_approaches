#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
install.packages("/Users/julied/Documents/R programming/WGCNA_1.49.tgz", repos = NULL, lib = .Library)

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "/Users/julied/Documents/R programming/WCGNA_Tutorial1/Data/FemaleLiver-Data";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");
# Take a quick look at what is in the data set:
dim(femData);
names(femData);


#=====================================================================================
#
#  Code chunk 2 Create a data frame and store names
#
#=====================================================================================


datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH;
rownames(datExpr0) = names(femData)[-c(1:8)];


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (gsg$allOK) {
  datExpr = datExpr0
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
}

#=====================================================================================
#
#  Code chunk 7 Read in trait data, select a subset
#
#=====================================================================================


traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Mice);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


#=====================================================================================
#
#  Code chunk 8 Look for correlation between expression data and trait values
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


#=====================================================================================
#
#  Code chunk 9 save data files
#
#=====================================================================================


save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")


