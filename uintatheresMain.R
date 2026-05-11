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
library(circlize)
library(mixtools)
library(patchwork)

Uintatheres <- read.csv("Uintathere_Measurements.csv", check.names = F) #load dataframe
Uintatheres <- Uintatheres[Uintatheres$Specimen != "YPM VP 011039L",] #remove damaged view of YPM VP 011039
rownames(Uintatheres) <- Uintatheres[,1] #rename rows as specimen numbers
uinFem <- Uintatheres$Specimen[Uintatheres$sex=="F"]
uinMale <- Uintatheres$Specimen[Uintatheres$sex=="M"]
uinNB <- Uintatheres$Specimen[Uintatheres$sex=="?"]
Uintatheres <- Uintatheres[,4:12] #remove helper variables

columns <- colnames(Uintatheres)
columns <- c("Length from Nasals to Occipital Crest", 
             "Length from Nasals to Occipital Condyles", 
             "Length from Premaxillaries to Occipital Condyles", 
             "Height from Occipital Crest to Occipital Condyles",
             "Height from Parietal Protuberance to tip of Postglenoid Process",
             "Height of Maxillary protuberance from Canine Alveolus", 
#             "Height from Maxillary Protuberance to Tip of Canine", #For now, I have commented out this measurement
             "Length of Canine from Alveolus to Tip",
             "Diastemal Length")
NVariables <- length(columns)

Uintatheres <- Uintatheres[,-7] #remove old measurement, length from top of maxillary protuberance to tip of canine

Measurement <- vector()
n <- vector()
wStatistic <- vector()
pValue <- vector()
dipStatistic <- vector() 

Stats <- data.frame(Measurement, n, wStatistic, pValue, dipStatistic, pValue)

#RUN for-loop that plots histograms for each variable, runs Shapiro-Wilk and 
#Hartigan's Dip Test, and stores these variables in "Stats" dataframe.

for (i in 1:NVariables) {
  histi <- hist(Uintatheres[,i][!is.na(Uintatheres[,i])], main = paste(columns[i]),
                col = "skyblue1", xlab = columns[i], cex.main = 2, cex.lab=1.5)
  lengthi <- unique(Uintatheres[,i][!is.na(Uintatheres[,i])])
  xfit <- seq(min(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              max(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              length = length(lengthi)) 
  yfit <- dnorm(xfit, mean = mean(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
                sd = sd(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T)) 
  yfit <- yfit * diff(histi$mids[1:2]) * length(Uintatheres[,i][!is.na(Uintatheres[,i])]) 
  lines(xfit, yfit, col = "red", lwd = 3)
  Stats[i,1] <- columns[i]
  Stats[i,2] <- length(Uintatheres[,i][!is.na(Uintatheres[,i])])
  Stats[i,3] <- shapiro.test(Uintatheres[,i])[1]
  Stats[i,4] <- shapiro.test(Uintatheres[,i])[2]
  Stats[i,5] <- dip.test(Uintatheres[,i])[1]
  Stats[i,6] <- dip.test(Uintatheres[,i])[2]
  mtext(paste("w-statistic: ", Stats[i,2], "; p-value: ", Stats[i,3], "\n",
              "dip-statistic: ", Stats[i,4], "; p-value: ", Stats[i,5],  sep = ""), 
        padj = .3, side = 4, cex = 2)
  #text(x=mean(Uintatheres[,i]), y=-0.5, paste("w-statistic: ", Stats[i,2], "; p-value: ", 
  #             Stats[i,3], "\n","dip-statistic: ", Stats[i,4], "; p-value: ", Stats[i,5],  sep = ""),
  #     cex = 15)
}

#Plot all 8 histograms together

par(mfrow = c(4,2))
par(mar=c(2,2,2,2))

for (i in 1:NVariables) {
  histi <- hist(log(Uintatheres[,i])[!is.na(Uintatheres[,i])], main = paste(columns[i]),
                col = "skyblue1", xlab = paste(columns[i], "(cm)"))
  lengthi <- unique(Uintatheres[,i][!is.na(Uintatheres[,i])])
  xfit <- seq(min(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              max(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              length = length(lengthi)) 
  yfit <- dnorm(xfit, mean = mean(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
                sd = sd(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T)) 
  yfit <- yfit * diff(histi$mids[1:2]) * length(Uintatheres[,i][!is.na(Uintatheres[,i])]) 
  lines(xfit, yfit, col = "red", lwd = 3)
  Stats[i,1] <- columns[i]
  Stats[i,2] <- length(Uintatheres[,i][!is.na(Uintatheres[,i])])
  Stats[i,3] <- shapiro.test(Uintatheres[,i])[1]
  Stats[i,4] <- shapiro.test(Uintatheres[,i])[2]
  Stats[i,5] <- dip.test(Uintatheres[,i])[1]
  Stats[i,6] <- dip.test(Uintatheres[,i])[2]
  mtext(paste("w-statistic: ", Stats[i,2], "; p-value: ", Stats[i,3], "\n",
              "dip-statistic: ", Stats[i,4], "; p-value: ", Stats[i,5],  sep = ""), 
        padj = .3, side = 4, cex = 0.45)
} #this is log trnasformed

par(mfrow = c(1,1))

#Generate Violin Plots for Different Measurements

par(mfrow = c(4,2))
par(mar=c(2,2,2,2))

dims <- list(length(colnames(Uintatheres)))
dims[[8]] <- ggplot(data = Uintatheres,
                  mapping = aes(x="", y= `nasals to occipital crest` )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Nasals to Occipital Crests")


for (i in 1:length(colnames(Uintatheres))) {
  dims[[i]] <-  ggplot(data = Uintatheres,
                       mapping = aes(x="", y=.data[[colnames(Uintatheres)[i]]])) +
    geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
    labs(y = "Length (cm)", x="", title = colnames(Uintatheres)[i])
}

#may want to apply log transformation to data above

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]],
             ncol = 3, nrow = 3)

set.seed(123)

MyMin <- min(Uintatheres$parPro, na.rm = T)
MyMax <- max(Uintatheres$parPro, na.rm = T)
MyMean <- mean(Uintatheres$parPro, na.rm = T)
MyVar <- var(Uintatheres$parPro, na.rm = T)
MySD <- sqrt(MyVar)
MyNorm <- rnorm(19, mean = MyMean, sd = MySD)
MyData <- data.frame(random = 1:19)
MyData$random <- MyNorm

random <- ggplot(data = MyData,
       mapping = aes(x="", y= `random` )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Random Normal Data")

random

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]], random,
             ncol = 3, nrow = 3)

mixtureModels <- data.frame(
  metric = character(),
  numberComponents = numeric(),
  BIC = numeric()
)

for (i in 1:length(colnames(Uintatheres))) {
  mixMod <- Mclust(Uintatheres[,i][!is.na(Uintatheres[,i])])
  assign(paste(colnames(Uintatheres)[i], "mixMod"), mixMod)
  mixtureModels[i,1] <- colnames(Uintatheres)[i]
  mixtureModels[i,2] <- mixMod$G
  mixtureModels[i,3] <- mixMod$bic
}

par(mfrow=c(1,1))


aggr(Uintatheres)
uinTrim <- Uintatheres[rowSums(is.na(Uintatheres)) < 4,]
aggr(uinTrim)
uinScale <- scale(uinTrim)
rownames(uinScale)[18] <- "YPM VPPU 010298"


ut.ncp <- estim_ncpPCA(uinScale, method.cv= "Kfold")
plot(ut.ncp$criterion~names(ut.ncp$criterion),xlab="number of dimensions", 
     ylab ="ncp criterion")
ncp <- as.numeric(ut.ncp$ncp[1]) 

uinImp <- imputePCA(uinScale, ncp = ncp)
uinImp <- uinImp$completeObs
uinImp <- as.data.frame(uinImp)

uinSex <- uinImp
uinSex$sex <- NA

uinFem <- c("AMNH FM 1671", "AMNH FM 1693", "YPM VP 011202", "YPM VPPU 010298")

uinSex$sex <- ifelse(rownames(uinSex) %in% uinFem, "F", 
                      ifelse(rownames(uinSex) %in% uinMale, "M", "?"))

row_colors <- ifelse(uinSex$sex == "F", "red", 
                     ifelse(uinSex$sex == "M", "blue", "black"))

Heatmap(uinImp,
        name = "Z-score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_names_rot = 65,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(col = row_colors, fontsize = 15),
        row_dend_width = unit(100, "points")
        )

uinImp$sex <- uinSex$sex
uinImp$sex <- as.factor(uinImp$sex)

uinImpPCA <- PCA(uinImp[,1:8], scale.unit=F)

write.csv(uinImpPCA$eig, "eigVals.csv")
write.csv(uinImpPCA$var$coord, "pcLoad.csv")
write.csv(uinImpPCA$var$cos2, "pcContVar.csv")
write.csv(uinImpPCA$var$contrib, "varContPc.csv")

fviz_pca_var(uinImpPCA, col.var = "contrib") 
fviz_pca_ind(uinImpPCA)
fviz_pca_biplot(uinImpPCA)

#Now, I'll revisualize the same PCA, this time grouping "females" and "not females" into
#two convex hulls to make the morphospatial pattern a bit clearer.

coords <- uinImpPCA$ind$coord
coords <- as.data.frame(coords)
isFem <- rownames(coords) %in% uinFem
coords$Fem <- isFem
coords$sex <- uinSex$sex
coords$sex <- as.factor(coords$sex)

#Now, let's rerun our statistical tests for our new PC's

pcCols <- c("wStat", "pVal", "dipStat", "pVal1")
pcStats <- matrix(nrow = 5, ncol = 4)
rownames(pcStats) <- colnames(coords)[1:5]
colnames(pcStats) <- pcCols
pcStats <- as.data.frame(pcStats)

for (i in 1:5) {
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

MyPlot <- ggplot(coords, aes(Dim.1, Dim.2)) + geom_point(shape = 21)
MyPlot + aes(fill = factor(sex)) + geom_polygon(data = hull_sex, alpha = 0.5) +
  labs(x = "PC1", y = "PC2", title = "PCA of Uintatherium Skulls")+ geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme(legend.title = element_blank()) +
  geom_label_repel(aes(label = rownames(coords))) +
  scale_fill_manual(values = c("#00BA38", "#F8766D", "#619CFF"))

gapStat <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 8, method = "gap_stat")
fviz_nbclust(uinImp[,1:8], kmeans, k.max = 8, method = "gap_stat")

gapStat <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 17, method = "gap_stat")
wss <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 17, method = "wss")
silh <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 17, method = "silhouette")

wss / silh / gapStat

km.UinTrimImp <- kmeans(uinImp[,1:8], centers = 2)
UintatheresClustered <- uinImp
clustVec <- km.UinTrimImp$cluster
UintatheresClustered <- as.data.frame(UintatheresClustered)
UintatheresClustered$cluster = clustVec
UintatheresClustered$cluster <- as.factor(UintatheresClustered$cluster)

ClusteredPCA <- PCA(UintatheresClustered[,1:8])

clusterCoords <-ClusteredPCA$ind$coord
clusterCoords <- as.data.frame(coords)
clusterCoords <- clusterCoords[,-c(6,7)]
clusterCoords$Cluster <- UintatheresClustered$cluster

hull_cluster <- clusterCoords %>%
  group_by(Cluster) %>%
  slice(chull(Dim.1, Dim.2))

MyClusterPlot <- ggplot(clusterCoords, aes(x = Dim.1, y = Dim.2)) + geom_point(shape = 21)

MyClusterPlot + aes(fill = Cluster) + geom_polygon(data = hull_cluster, alpha = 0.5) +
  labs(x = "PC1", y = "PC2", colour = "Cluster", title = "PCA of Uintatherium Skulls") + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  geom_label_repel(aes(label = rownames(clusterCoords))) +
  scale_fill_manual(values = c("#619CFF", "#F8766D"))

grid.arrange(MyPlot + aes(fill = factor(sex)) + geom_polygon(data = hull_sex, alpha = 0.5) +
               labs(x = "PC1", y = "PC2", title = "PCA of Uintathere Skulls", colour = "Sex")+ 
               geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
               theme(legend.title = element_blank()) +
               geom_label_repel(aes(label = rownames(coords))) +
               scale_fill_manual(values = c("#00BA38", "#F8766D", "#619CFF")),
             MyClusterPlot + aes(fill = Cluster) + geom_polygon(data = hull_cluster, alpha = 0.5) +
               labs(x = "PC1", y = "PC2", colour = "Cluster", title = "PCA of Uintathere Skulls (Clustered)") + 
               geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
               geom_label_repel(aes(label = rownames(clusterCoords))) +
               scale_fill_manual(values = c("#619CFF", "#F8766D")),
             ncol = 1)

Uintatheres <- read.csv("Uintathere_Measurements_With_Turnbull.csv", check.names = F) #load dataframe
Uintatheres <- Uintatheres[Uintatheres$Specimen != "YPM VP 011039L",] #remove damaged view of YPM 11039
rownames(Uintatheres) <- Uintatheres[,1] #rename rows as specimen numbers
uinFem <- Uintatheres$Specimen[Uintatheres$sex=="F"]
uinMale <- Uintatheres$Specimen[Uintatheres$sex=="M"]
uinNB <- Uintatheres$Specimen[Uintatheres$sex=="?"]
Uintatheres <- Uintatheres[,4:12] #remove helper variables
#Uintatheres[1,]$alvCan <- NA #This canine is potentially pathologic so it's inclusion is questionable

str(Uintatheres)

colnames(Uintatheres)

columns <- colnames(Uintatheres)
columns <- c("Length from Nasals to Occipital Crest", 
             "Length from Nasals to Occipital Condyles", 
             "Length from Premaxillaries to Occipital Condyles", 
             "Height from Occipital Crest to Occipital Condyles",
             "Height from Parietal Protuberance to tip of Postglenoid Process",
             "Height of Maxillary protuberance from Canine Alveolus", 
             #             "Height from Maxillary Protuberance to Tip of Canine", #For now, I have commented out this measurement
             "Length of Canine from Alveolus to Tip",
             "Diastemal Length")
NVariables <- length(columns)

Uintatheres <- Uintatheres[,-7] #get rid of maxCan measurement

Measurement <- vector()
n <- vector()
wStatistic <- vector()
pValue <- vector()
dipStatistic <- vector() 

Stats <- data.frame(Measurement, n, wStatistic, pValue, dipStatistic, pValue)

#RUN for-loop that plots histograms for each variable, runs Shapiro-Wilk and 
#Hartigan's Dip Test, and stores these variables in "Stats" dataframe.

for (i in 1:NVariables) {
  histi <- hist(Uintatheres[,i][!is.na(Uintatheres[,i])], main = paste(columns[i]),
                col = "skyblue1", xlab = columns[i], cex.main = 2, cex.lab=1.5)
  lengthi <- unique(Uintatheres[,i][!is.na(Uintatheres[,i])])
  xfit <- seq(min(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              max(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              length = length(lengthi)) 
  yfit <- dnorm(xfit, mean = mean(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
                sd = sd(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T)) 
  yfit <- yfit * diff(histi$mids[1:2]) * length(Uintatheres[,i][!is.na(Uintatheres[,i])]) 
  lines(xfit, yfit, col = "red", lwd = 3)
  Stats[i,1] <- columns[i]
  Stats[i,2] <- length(Uintatheres[,i][!is.na(Uintatheres[,i])])
  Stats[i,3] <- shapiro.test(Uintatheres[,i])[1]
  Stats[i,4] <- shapiro.test(Uintatheres[,i])[2]
  Stats[i,5] <- dip.test(Uintatheres[,i])[1]
  Stats[i,6] <- dip.test(Uintatheres[,i])[2]
  mtext(paste("w-statistic: ", Stats[i,2], "; p-value: ", Stats[i,3], "\n",
              "dip-statistic: ", Stats[i,4], "; p-value: ", Stats[i,5],  sep = ""), 
        padj = .3, side = 4, cex = 2)
  #text(x=mean(Uintatheres[,i]), y=-0.5, paste("w-statistic: ", Stats[i,2], "; p-value: ", 
  #             Stats[i,3], "\n","dip-statistic: ", Stats[i,4], "; p-value: ", Stats[i,5],  sep = ""),
  #     cex = 15)
}

#Plot all 8 histograms together

par(mfrow = c(4,2))
par(mar=c(2,2,2,2))

for (i in 1:NVariables) {
  histi <- hist(log(Uintatheres[,i])[!is.na(Uintatheres[,i])], main = paste(columns[i]),
                col = "skyblue1", xlab = paste(columns[i], "(cm)"))
  lengthi <- unique(Uintatheres[,i][!is.na(Uintatheres[,i])])
  xfit <- seq(min(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              max(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
              length = length(lengthi)) 
  yfit <- dnorm(xfit, mean = mean(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T), 
                sd = sd(Uintatheres[,i][!is.na(Uintatheres[,i])], na.rm = T)) 
  yfit <- yfit * diff(histi$mids[1:2]) * length(Uintatheres[,i][!is.na(Uintatheres[,i])]) 
  lines(xfit, yfit, col = "red", lwd = 3)
  Stats[i,1] <- columns[i]
  Stats[i,2] <- length(Uintatheres[,i][!is.na(Uintatheres[,i])])
  Stats[i,3] <- shapiro.test(Uintatheres[,i])[1]
  Stats[i,4] <- shapiro.test(Uintatheres[,i])[2]
  Stats[i,5] <- dip.test(Uintatheres[,i])[1]
  Stats[i,6] <- dip.test(Uintatheres[,i])[2]
  mtext(paste("w-statistic: ", Stats[i,2], "; p-value: ", Stats[i,3], "\n",
              "dip-statistic: ", Stats[i,4], "; p-value: ", Stats[i,5],  sep = ""), 
        padj = .3, side = 4, cex = 0.45)
} #this is log trnasformed, just so's ya know

par(mfrow = c(1,1))

#Generate Violin Plots for Different Measurements

par(mfrow = c(4,2))
par(mar=c(2,2,2,2))

dims <- list(length(colnames(Uintatheres)))
dims[[8]] <- ggplot(data = Uintatheres,
                    mapping = aes(x="", y= `nasals to occipital crest` )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Nasals to Occipital Crests")


for (i in 1:length(colnames(Uintatheres))) {
  dims[[i]] <-  ggplot(data = Uintatheres,
                       mapping = aes(x="", y=.data[[colnames(Uintatheres)[i]]])) +
    geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
    labs(y = "Length (cm)", x="", title = colnames(Uintatheres)[i])
}

#may want to aqpply log transformation to data above

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]],
             ncol = 3, nrow = 3)

set.seed(123)

MyMin <- min(Uintatheres$parPro, na.rm = T)
MyMax <- max(Uintatheres$parPro, na.rm = T)
MyMean <- mean(Uintatheres$parPro, na.rm = T)
MyVar <- var(Uintatheres$parPro, na.rm = T)
MySD <- sqrt(MyVar)
MyNorm <- rnorm(19, mean = MyMean, sd = MySD)
MyData <- data.frame(random = 1:19)
MyData$random <- MyNorm

random <- ggplot(data = MyData,
                 mapping = aes(x="", y= `random` )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Random Normal Data")

random

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]], random,
             ncol = 3, nrow = 3)

#Attempting some mixture modeling

colnames(Uintatheres)

mixtureModels <- data.frame(
  metric = character(),
  numberComponents = numeric(),
  BIC = numeric()
)

for (i in 1:length(colnames(Uintatheres))) {
  mixMod <- Mclust(Uintatheres[,i][!is.na(Uintatheres[,i])])
  assign(paste(colnames(Uintatheres)[i], "mixMod"), mixMod)
  mixtureModels[i,1] <- colnames(Uintatheres)[i]
  mixtureModels[i,2] <- mixMod$G
  mixtureModels[i,3] <- mixMod$bic
}

# write.csv(Stats, "uintaStats.csv")
# write.csv(mixtureModels, "uintaMixtureModels.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(1,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PART 2: Run a PCA to determine/visualize how data look when multiple variables 
#are considered at once.

#First, I have to impute missing data for PCA. Since I'm only interested in
#filling in single point values, single imputation will work fine.

######WORKING WITH PROPPER IMPUTATION PROTOCOLS (Impute if 3 or fewer missing, otherwise, remove specimen)

uinTrim <- Uintatheres[rowSums(is.na(Uintatheres)) < 4,]
aggr(uinTrim)
uinScale <- scale(uinTrim)
rownames(uinScale)[18] <- "YPM VPPU 010298"
ut.ncp <- estim_ncpPCA(uinScale, method.cv= "Kfold")
plot(ut.ncp$criterion~names(ut.ncp$criterion),xlab="number of dimensions", 
     ylab ="ncp criterion")
ncp <- as.numeric(ut.ncp$ncp[1]) 

uinImp <- imputePCA(uinScale, ncp = ncp)
uinImp <- uinImp$completeObs
uinImp <- as.data.frame(uinImp)

uinFem <- c("AMNH FM 1671", "AMNH FM 1693", "YPM VP 011202", "YPM VPPU 010298")

uinSex <- uinImp
uinSex$sex <- NA

uinSex$sex <- ifelse(rownames(uinSex) %in% uinFem, "F", 
                     ifelse(rownames(uinSex) %in% uinMale, "M", "?"))

uinSexTrim <- uinSex[rownames(uinSex) != "YPM 11039R",]
row_colors <- ifelse(uinSexTrim$sex == "F", "red", 
                     ifelse(uinSexTrim$sex == "M", "blue", "black"))

Heatmap(uinSexTrim[,1:8],
        name = "Z-score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_names_rot = 65,
        column_names_gp = gpar(fontsize = 12),
        row_names_gp = gpar(col = row_colors, fontsize = 20))


uinImp$sex <- uinSex$sex
uinImp$sex <- as.factor(uinImp$sex)

uinImpPCA <- PCA(uinImp[,1:8], scale.unit=F)

# write.csv(uinImpPCA$eig, "uintaEigVals.csv")
# write.csv(uinImpPCA$var$coord, "uintaPCLoad.csv")
# write.csv(uinImpPCA$var$cos2, "uintaPCContVar.csv")
# write.csv(uinImpPCA$var$contrib, "uintaVarContPc.csv")

uinImpPCA$eig

fviz_pca_var(uinImpPCA, col.var = "contrib") 
fviz_pca_ind(uinImpPCA)
fviz_pca_biplot(uinImpPCA)

#Now, I'll revisualize the same PCA, this time grouping "females" and "not females" into
#two convex hulls to make the morphospatial pattern a bit clearer.

coords <- uinImpPCA$ind$coord
coords <- as.data.frame(coords)
isFem <- rownames(coords) %in% uinFem
coords$Fem <- isFem
coords$sex <- uinSex$sex
coords$sex <- as.factor(coords$sex)

hull_sex <- coords %>%
  group_by(sex) %>%
  slice(chull(Dim.1, Dim.2))

coords$Sex <- as.factor(coords$sex)

MyPlot <- ggplot(coords, aes(Dim.1, Dim.2)) + geom_point(shape = 21)

MyPlot +
  geom_polygon(data = hull_sex, aes(group = sex, fill = sex), alpha = 0.5) +
  labs(x = "PC1", y = "PC2", title = "PCA of Uintatherium Skulls") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_label_repel(aes(label = rownames(coords))) +
  scale_fill_manual(
    name = "Sex",  # legend title
    values = c(
      "F" = "#F8766D",
      "M" = "#619CFF",
      "?" = "#B8B8B8"
    )
  )

#Now, let's rerun our statistical tests for our new PC's

pcCols <- c("wStat", "pVal", "dipStat", "pVal1")
pcStats <- matrix(nrow = 3, ncol = 4)
rownames(pcStats) <- colnames(coords)[1:3]
colnames(pcStats) <- pcCols
pcStats <- as.data.frame(pcStats)

for (i in 1:3) {
  wTest <- shapiro.test(coords[,i])
  pcStats[i,]$wStat <- wTest$statistic
  pcStats[i,]$pVal <- wTest$p.value
  dipTest <- dip.test(coords[,i])
  pcStats[i,]$dipStat <- dipTest$statistic
  pcStats[i,]$pVal1 <- dipTest$p.value
}


#SO, it seems like the "females" occupy a different morphospace than the "not females", scoring
#negative PC1 and PC2 values. Cool! So we have evidence for dimorphism after all right?? Well...
#if that were the case, that would mean that only 3/17 skulls in this dataset are female. There's no reason 
#to say that uintathere sex ratios favored males, so let's look at a normal distribution to see
#how rare a sample like this would be from a 50/50 population.

ratioNorm <- rnorm(1000000, 0.5, 0.121)
hist(ratioNorm, main = "Draws from a 50:50 Pouplation", xlab = "Proportion of Females",
     ylab = "Frequency")
abline(v=3/23, col = "red", lwd = 2)
(sum(ratioNorm > 20/23) + sum(ratioNorm<3/23))/length(ratioNorm)
(sum(ratioNorm > 11/13) + sum(ratioNorm<3/13))/length(ratioNorm)
#The probability of observing a 14:3 male/female sample from a 50/50 population is only about 0.00223.

#Let's do a K-means clustering analysis based on the gap statistic to see how many clusters 
#actually seem to be present in the data.

gapStat <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 18, method = "gap_stat")
wss <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 18, method = "wss")
silh <- fviz_nbclust(uinImp[,1:8], kmeans, k.max = 18, method = "silhouette")

wss / silh / gapStat

#As is, there's no statistical evidence that the data forms more than one cluster.

#Let's just see how the PCA looks when the data is forced to form 2 clusters...

km.UinTrimImp <- kmeans(uinImp[,1:8], centers = 2)
UintatheresClustered <- uinImp
clustVec <- km.UinTrimImp$cluster
UintatheresClustered <- as.data.frame(UintatheresClustered)
UintatheresClustered$cluster = clustVec
UintatheresClustered$cluster <- as.factor(UintatheresClustered$cluster)

ClusteredPCA <- PCA(UintatheresClustered[,1:8])

clusterCoords <-ClusteredPCA$ind$coord
clusterCoords <- as.data.frame(coords)
clusterCoords <- clusterCoords[,-c(6,7)]
clusterCoords$Cluster <- UintatheresClustered$cluster

clusterCoords$Cluster <- UintatheresClustered$cluster

hull_cluster <- clusterCoords %>%
  group_by(Cluster) %>%
  slice(chull(Dim.1, Dim.2))

MyClusterPlot <- ggplot(clusterCoords, aes(x = Dim.1, y = Dim.2)) + geom_point(shape = 21)

hull_sex$sex <- as.factor(hull_sex$sex)
hull_sex$Sex <- hull_sex$sex

bothPC <- grid.arrange(MyPlot + aes(fill = Sex) + geom_polygon(data = hull_sex, alpha = 0.5) +
                         labs(x = "PC1", y = "PC2", title = "PCA of Uintathere Skulls")+ 
                         geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
                         #               theme(legend.title = element_blank()) +
                         geom_label_repel(size = 3,
                                          force = 30,
                                          aes(label = rownames(coords))) +
                         scale_fill_manual(values = c("#B8B8B8", "#F8766D", "#619CFF")),
                       MyClusterPlot + aes(fill = Cluster) + geom_polygon(data = hull_cluster, alpha = 0.5) +
                         labs(x = "PC1", y = "PC2", title = "PCA of Uintathere Skulls (Clustered)") + 
                         geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
                         geom_label_repel(size = 3,
                                          force = 30,
                                          aes(label = rownames(coords))) +
                         scale_fill_manual(values = c("#C77CFF", "#00BA38")),
                       ncol = 1)

######## testing for correlation between sex and cluster

uintatheresFisher <- UintatheresClustered[UintatheresClustered$sex != "?",]
uintatheresFisher <- uintatheresFisher[,9:10]
contTable <- table(uintatheresFisher$sex, uintatheresFisher$cluster)
fisher.test(contTable)

statsRound <- Stats
statsRound[,3:6] <- round(statsRound[,3:6], 3)

write.csv(statsRound, "uintaStatsRound.csv")
