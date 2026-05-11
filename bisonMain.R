#setwd("")

library(dplyr)
library(diptest)
library(Amelia)
library(cluster)
library(ggfortify)
library(ggplot2)
library(ggforce)
library(factoextra)
library(FactoMineR)
library(missMDA)
library(VIM)
library(cowplot)
library(plyr)
library(ggrepel)
library(gridExtra)
library(mclust)
library(GGally)
library(BayesFactor)
library(pheatmap)
library(ComplexHeatmap)
#library(circlize)
library(mixtools)
library(patchwork)

bison <- read.csv("Bison_Measurements.csv", check.names = F) #load dataframe
bisonAdults <- bison[-20, ] #remove 68350, specimen with adult dentition still erupting
rownames(bisonAdults) <- bisonAdults[,1] #rename rows as specimen numbers
adultSex <- bisonAdults[,2] #make a vector of sex for adult specimens
bisonAdults <- bisonAdults[,-c(1:4)] #remove helper variables
for (i in 1:ncol(bisonAdults)) {
  bisonAdults[,i] <- as.numeric(bisonAdults[,i])
}

str(bisonAdults)

bisonAdults <- bisonAdults[,-4] #remove maxMax to retain 8 variables

colnames(bisonAdults)

columns <- c("Width between tips of frontal Horns", "Width bewteen bases of Frontal Horns", 
             "Width between exorbitale",
             "Length from prosthion to nasion", "Length from nasion to nuchal crest", 
             "Length of frontal horn from base to tip", "Length of skull from orbit to prosthion", 
             "Width of frontal horn at base")
NVariables <- length(columns)
Measurement <- vector()
n <- vector()
wStatistic <- vector()
pValue <- vector()
dipStatistic <- vector() 

Stats <- data.frame(Measurement, n, wStatistic, pValue, dipStatistic, pValue)

#RUN for-loop that plots histograms for each variable, runs Shapiro-Wilk and 
#Hartigan's Dip Test, and stores these variables in "Stats" dataframe.
#Plot all 8 histograms together

par(mfrow = c(4,2))
par(mar=c(2,2,2,2))

for (i in 1:NVariables) {
  histi <- hist(bisonAdults[,i][!is.na(bisonAdults[,i])], main = paste(columns[i]),
                col = "skyblue1", xlab = paste(columns[i], "(cm)"))
  lengthi <- unique(bisonAdults[,i][!is.na(bisonAdults[,i])])
  xfit <- seq(min(bisonAdults[,i][!is.na(bisonAdults[,i])], na.rm = T), 
              max(bisonAdults[,i][!is.na(bisonAdults[,i])], na.rm = T), 
              length = length(lengthi)) 
  yfit <- dnorm(xfit, mean = mean(bisonAdults[,i][!is.na(bisonAdults[,i])], na.rm = T), 
                sd = sd(bisonAdults[,i][!is.na(bisonAdults[,i])], na.rm = T)) 
  yfit <- yfit * diff(histi$mids[1:2]) * length(bisonAdults[,i][!is.na(bisonAdults[,i])]) 
  lines(xfit, yfit, col = "red", lwd = 3)
  Stats[i,1] <- columns[i]
  Stats[i,2] <- length(bisonAdults[,i][!is.na(bisonAdults[,i])])
  Stats[i,3] <- shapiro.test(bisonAdults[,i])[1]
  Stats[i,4] <- shapiro.test(bisonAdults[,i])[2]
  Stats[i,5] <- dip.test(bisonAdults[,i])[1]
  Stats[i,6] <- dip.test(bisonAdults[,i])[2]
  mtext(paste("w-statistic: ", Stats[i,3], "; p-value: ", Stats[i,4], "\n",
              "dip-statistic: ", Stats[i,5], "; p-value: ", Stats[i,6],  sep = ""), 
        padj = .3, side = 4, cex = 0.45)
}

par(mfrow = c(1,1))

#Generate Violin Plots for Different Measurements

par(mfrow = c(4,2))
par(mar=c(2,2,2,2))

dims <- list(length(colnames(bisonAdults)))
dims[[8]] <- ggplot(data = bisonAdults,
                  mapping = aes(x="", y= 'width between tips of frontal Horns' )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Width between tips of frontal Horns")


for (i in 1:length(colnames(bisonAdults))) {
  dims[[i]] <-  ggplot(data = bisonAdults,
                       mapping = aes(x="", y=.data[[colnames(bisonAdults)[i]]])) +
    geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
    labs(y = "Length (cm)", x="", title = colnames(bisonAdults)[i])
}

#may want to apply log transformation to data above

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]],
             ncol = 3, nrow = 3)

set.seed(123)

MyMin <- min(bisonAdults$tipFro, na.rm = T)
MyMax <- max(bisonAdults$tipFro, na.rm = T)
MyMean <- mean(bisonAdults$tipFro, na.rm = T)
MyVar <- var(bisonAdults$tipFro, na.rm = T)
MySD <- sqrt(MyVar)
MyNorm <- rnorm(22, mean = MyMean, sd = MySD)
MyData <- data.frame(random = 1:22)
MyData$random <- MyNorm

random <- ggplot(data = MyData,
       mapping = aes(x="", y= `random` )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Random Normal Data")

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]], random,
             ncol = 3, nrow = 3)

colnames(bisonAdults)

mixtureModels <- data.frame(
  metric = character(),
  numberComponents = numeric(),
  BIC = numeric()
)

for (i in 1:length(colnames(bisonAdults))) {
  mixMod <- Mclust(bisonAdults[,i][!is.na(bisonAdults[,i])])
  assign(paste(colnames(bisonAdults)[i], "mixMod"), mixMod)
  mixtureModels[i,1] <- colnames(bisonAdults)[i]
  mixtureModels[i,2] <- mixMod$G
  mixtureModels[i,3] <- mixMod$bic
}

write.csv(Stats, "bisonStats.csv")
write.csv(mixtureModels, "bisonMixtureModels.csv")

par(mfrow=c(1,1))


bisonAdults <- bisonAdults[,2:7] #remove extra horn variables

bis.ncp <- estim_ncpPCA(bisonAdults, method.cv= "Kfold")
plot(bis.ncp$criterion~names(bis.ncp$criterion),xlab="number of dimensions", 
     ylab = "ncp criterion")
ncp <- as.numeric(bis.ncp$ncp[1]) #The optimal number of dimensions(components) which should considered when imputing data, as obtained by k-fold cross validation

bis.impute <- imputePCA(bisonAdults, ncp = ncp)
bis.impute <- bis.impute$completeObs
bis.impute <- as.data.frame(bis.impute)

bis.impute #Missing values are now imputed according to the PCA cross-validation.


bis.impute$sex <- adultSex

bis.impute$sex <- as.factor(bis.impute$sex)

bisNum <- bis.impute[,1:6]

bisScale <- scale(bisNum) #need to scale data before running PCA

bisScale <- as.data.frame(bisScale)

#First, let's just make a heat map(hierarchical cluster plot) to see how skull tend to cluster together

pheatmap(bisScale, clustering_method = "ward.D2", show_rownames = T, show_colnames = T,
         main = "Heatmap with Dendrogram")

#Making a fancier heatmap

row_colors <- ifelse(adultSex == "F", "red", 
                     ifelse(adultSex == "M", "blue", "black"))
Heatmap(bisScale,
        name = "Z-score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_names_rot = 65,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(col = row_colors, fontsize = 15),
        row_dend_width = unit(100, "points"))  # <- color the labels

heatMapClusters <- data.frame(
  specimen = character(23),
  sex = character(23),
  cluster = character(23)
)

heatMapClusters$specimen <- rownames(bisScale)
heatMapClusters$sex <- adultSex

cluster1 <- c(75033, 140290, 21906, 16211, 2414, 92698, 117340, 2022)
cluster1 <- as.character(cluster1)
cluster2 <- heatMapClusters$specimen[!heatMapClusters$specimen %in% cluster1]
for (i in 1:nrow(heatMapClusters)) {
  if(heatMapClusters$specimen[i] %in% cluster1){
    heatMapClusters$cluster[i] <- "1"
  } else if (heatMapClusters$specimen[i] %in% cluster2){
    heatMapClusters$cluster[i] <- "2"
  }
}

contTable <- table(heatMapClusters$sex, heatMapClusters$cluster)
fisher.test(contTable)

normPCA <- PCA(bisScale)

#write.csv(normPCA$eig, "newBisonEigVals.csv")
#write.csv(normPCA$var$coord, "neBisonPCLoad.csv")
# write.csv(normPCA$var$cos2, "bisonPCContVar.csv")
# write.csv(normPCA$var$contrib, "bisonVarContPc.csv")
# 
fviz_pca_ind(normPCA)
fviz_pca_biplot(normPCA)
fviz_pca_var(normPCA, col.var = "contrib")


coords <- normPCA$ind$coord
coords <- as.data.frame(coords)
bisonSex <- bisonAdults
bisonSex$sex <- adultSex
isFem <- bisonSex$sex=="F"
coords$Fem <- isFem
coords$sex <- adultSex
coords$sex <- as.factor(coords$sex)

#Now, let's rerun our statistical tests for our new PC's

pcCols <- c("wStat", "pVal", "dipStat", "pVal1")
pcStats <- matrix(nrow = 2, ncol = 4)
rownames(pcStats) <- colnames(coords)[1:2]
colnames(pcStats) <- pcCols
pcStats <- as.data.frame(pcStats)

for (i in 1:2) {
  wTest <- shapiro.test(coords[,i])
  pcStats[i,]$wStat <- wTest$statistic
  pcStats[i,]$pVal <- wTest$p.value
  dipTest <- dip.test(coords[,i])
  pcStats[i,]$dipStat <- dipTest$statistic
  pcStats[i,]$pVal1 <- dipTest$p.value
}

#Back to plotting

hull_sex <- coords %>%
  group_by(sex) %>%
  slice(chull(Dim.1, Dim.2))

coords$Sex <- as.factor(coords$sex)

MyPlot <- ggplot(coords, aes(Dim.1, Dim.2)) + geom_point(shape = 21)
MyPlot + aes(fill = factor(sex)) + geom_polygon(data = hull_sex, alpha = 0.5) +
  labs(x = "PC1", y = "PC2", title = "PCA of Bison Skulls")+ geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme(legend.title = element_blank()) +
  geom_label_repel(aes(label = rownames(coords))) +
  scale_fill_manual(values = c("#00BA38", "#F8766D", "#619CFF"))

fviz_nbclust(bisScale, kmeans, k.max = 2, method = "gap_stat")

# write.csv(normPCA$eig, "eigVals.csv")
# write.csv(normPCA$var$coord, "pcLoad.csv")
# write.csv(normPCA$var$cos2, "pcContVar.csv")
# write.csv(normPCA$var$contrib, "varContPc.csv")

gapStat <- fviz_nbclust(bisScale[,1:6], kmeans, k.max = 22, method = "gap_stat")
wss <- fviz_nbclust(bisScale[,1:6], kmeans, k.max = 22, method = "wss")
silh <- fviz_nbclust(bisScale[,1:6], kmeans, k.max = 22, method = "silhouette")

wss / silh / gapStat

kmBison <- kmeans(bisScale[,1:6], centers = 2)
bisClust <- bisScale
clustVec <- kmBison$cluster
bisClust <- as.data.frame(bisClust)
bisClust$cluster = clustVec
bisClust$cluster <- as.factor(bisClust$cluster)

ClusteredPCA <- PCA(bisClust[,1:6])

clusterCoords <-ClusteredPCA$ind$coord
clusterCoords <- as.data.frame(coords)
clusterCoords <- clusterCoords[,-c(6,7)]
clusterCoords$cluster <- bisClust$cluster
clusterCoords$Cluster <- as.factor(clusterCoords$cluster)

hull_cluster <- clusterCoords %>%
  group_by(cluster) %>%
  slice(chull(Dim.1, Dim.2))

MyClusterPlot <- ggplot(clusterCoords, aes(x = Dim.1, y = Dim.2)) + geom_point(shape = 21)

MyClusterPlot + aes(fill = cluster) + geom_polygon(data = hull_cluster, alpha = 0.5) +
  labs(x = "PC1", y = "PC2", colour = "Cluster", title = "PCA of Bison Skulls") + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_label_repel(aes(label = rownames(clusterCoords))) +
  scale_fill_manual(values = c("#619CFF", "#F8766D"))

hull_sex$Sex <- as.factor(hull_sex$sex)
hull_cluster$Cluster <- as.factor(hull_cluster$cluster)

bothPC <- grid.arrange(MyPlot + geom_polygon(data = hull_sex, alpha = 0.5) + aes(fill = Sex) +
                         labs(x = "PC1 (74.54%)", y = "PC2 (18.47%)", 
                              title = "PCA of Bison Skulls")+ geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
#                         theme(legend.title = element_blank()) +
                         geom_label_repel(size = 3,
                                          force = 30,
                                          aes(label = rownames(coords))) +
                         scale_fill_manual(values = c("#B8B8B8", "#F8766D", "#619CFF")),
                       MyClusterPlot + aes(fill = Cluster) + geom_polygon(data = hull_cluster, alpha = 0.5) +
                         labs(x = "PC1 (74.54%)", y = "PC2 (18.47%)", 
                              title = "PCA of Bison Skulls (Clustered)") + 
                         geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
                         geom_label_repel(size = 3,
                                          force = 30,
                                          aes(label = rownames(coords))) +
                         scale_fill_manual(values = c("#C77CFF", "#00BA38")),
                       ncol = 1)

bisClust$sex <- coords$sex
# bisonFisher <- bisClust[bisClust$sex != "?",]
# bisonFisher <- bisonFisher[,7:8]
# contTable <- table(bisonFisher$sex, bisonFisher$cluster)
# fisher.test(contTable)

contTable <- table(bisClust$sex, bisClust$cluster)
fisher.test(contTable)
