setwd("")

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
library(EnvStats)
library(dimorph)
library(effectsize)
library(pvclust)
library(ggnewscale)
library(MASS)
library(devtools)
# library(ggord)
library(psych)
library(klaR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########Load and manipulate data

Uintatheres <- read.csv("Uintathere_Measurements_With_Canine_Width.csv", check.names = F) #load dataframe
Uintatheres <- Uintatheres[Uintatheres$Specimen != "YPM VP 011039L",] #remove damaged view of YPM 11039
rownames(Uintatheres) <- Uintatheres[,1] #rename rows as specimen numbers
uinFem <- Uintatheres$Specimen[Uintatheres$sex=="F"]
uinMale <- Uintatheres$Specimen[Uintatheres$sex=="M"]
uinNB <- Uintatheres$Specimen[Uintatheres$sex=="?"]
Uintatheres <- Uintatheres[,4:13] #remove helper variables
#Uintatheres[1,]$alvCan <- NA #This canine is potentially pathologic so it's inclusion is questionable

str(Uintatheres)

colnames(Uintatheres)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PART 1: Generate histograms for each variable and test
# for normality and unimodality

#Setting up "Stats" dataframe to store values of various test statistics of each measurement

columns <- colnames(Uintatheres)
columns <- c("Length from Nasals to Occipital Crest", 
             "Length from Nasals to Occipital Condyles", 
             "Length from Premaxillaries to Occipital Condyles", 
             "Height from Occipital Crest to Occipital Condyles",
             "Height from Parietal Protuberance to tip of Postglenoid Process",
             "Height of Maxillary protuberance from Canine Alveolus", 
             #             "Height from Maxillary Protuberance to Tip of Canine", #For now, I have commented out this measurement
             "Length of Canine from Alveolus to Tip",
             "Diastemal Length",
             "Canine Width at Base")
NVariables <- length(columns)

Uintatheres <- Uintatheres[,-7] #Remove old variable

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
  mtext(paste("w-statistic: ", Stats[i,3], "; p-value: ", Stats[i,4], "\n",
              "dip-statistic: ", Stats[i,5], "; p-value: ", Stats[i,6],  sep = ""), 
        padj = .3, side = 4, cex = 2)
}

par(mfrow = c(1,1))

#Generate Violin Plots for Each Measurements

par(mfrow = c(3,3))
par(mar=c(2,2,2,2))

dims <- list(length(colnames(Uintatheres)))
dims[[9]] <- ggplot(data = Uintatheres,
                    mapping = aes(x="", y= `nasals to occipital crest` )) +
  geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
  labs(y = "Length (cm)", x="", title = "Nasals to Occipital Crests")


for (i in 1:length(colnames(Uintatheres))) {
  dims[[i]] <-  ggplot(data = Uintatheres,
                       mapping = aes(x="", y=.data[[colnames(Uintatheres)[i]]])) +
    geom_violin( fill = "gray80", trim = F) + geom_point() + geom_boxplot(width = 0.1, fill = "gray70") +
    labs(y = "Length (cm)", x="", title = colnames(Uintatheres)[i])
}


grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]], dims[[9]],
             ncol = 3, nrow = 3)

set.seed(123) #generating random normal data for reference

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

#plot with random normal data

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]], dims[[9]], random,
             ncol = 2, nrow = 5)


#Mixture modeling


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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(1,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PART 2: Multivariate analysis

#First, missing data must be imputed for PCA. 
######WORKING WITH PROPPER IMPUTATION PROTOCOLS (Impute if 3 or fewer missing, otherwise, remove specimen)

aggr(Uintatheres)
uinTrim <- Uintatheres[rowSums(is.na(Uintatheres)) < 5,]
aggr(uinTrim)
uinScale <- scale(uinTrim)


ut.ncp <- estim_ncpPCA(uinScale, method.cv= "Kfold")
plot(ut.ncp$criterion~names(ut.ncp$criterion),xlab="number of dimensions", 
     ylab ="ncp criterion")
ncp <- as.numeric(ut.ncp$ncp[1]) 

uinImp <- imputePCA(uinScale, ncp = ncp)
uinImp <- uinImp$completeObs
uinImp <- as.data.frame(uinImp)

####Hierarchical Clustering

uinSex <- uinImp
uinSex$sex <- NA

uinSex$sex <- ifelse(rownames(uinSex) %in% uinFem, "F", 
                     ifelse(rownames(uinSex) %in% uinMale, "M", "?"))

row_colors <- ifelse(uinSex$sex == "F", "red", 
                     ifelse(uinSex$sex == "M", "blue", "black"))

Heatmap(uinImp,
        name = "Z-score",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        column_names_rot = 65,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(col = row_colors, fontsize = 15),
        row_dend_width = unit(75, "points")
)

###ALTERNATIVE CLUSTERING WITH BOOTSTRAP SUPPORT VALUES

bootMap <- pvclust(t(uinImp), method.dist="euclidean", method.hclust="complete", nboot=1000)
plot(bootMap)
pvrect(bootMap, alpha=0.95)

par(mar = c(4.5, 2, 2, 2))
plot(bootMap, labels = FALSE, hang = -5, main = "Hierarchical Clustering of Uintathere Skulls", print.num = F,
     xlab = "", sub="")
tipLabels <- bootMap$hclust$labels[bootMap$hclust$order]
tipColors <- ifelse(bootMap$hclust$labels[bootMap$hclust$order] %in% uinFem, "red",
                    ifelse(bootMap$hclust$labels[bootMap$hclust$order] %in% uinMale, "blue", "black"))
text(x = 1:length(tipLabels), y = -0.02, labels = tipLabels, col = tipColors,
     srt = 90,adj = 1,xpd = TRUE, cex = 0.75, font = 2)
dev.off()


uinImp$sex <- uinSex$sex
uinImp$sex <- as.factor(uinImp$sex)

###MANOVA AND R2 FOR HEATMAP

uinHeatMan <- manova(cbind(nasOcc, nasCon, preCon, occCon, parPro, maxAlv, 
                           alvCan, diaLen, canWid) ~ sex, data = uinImp)
summary(uinHeatMan, test = "Pillai")
summary.aov(uinHeatMan)
eta_squared(uinHeatMan, partial = T)

#####PCA

uinImpPCA <- PCA(uinImp[,1:9], scale.unit=F)

# write.csv(uinImpPCA$eig, "eigVals.csv")
# write.csv(uinImpPCA$var$coord, "pcLoad.csv")
# write.csv(uinImpPCA$var$cos2, "pcContVar.csv")
# write.csv(uinImpPCA$var$contrib, "varContPc.csv")

fviz_pca_var(uinImpPCA, col.var = "contrib") 
fviz_pca_ind(uinImpPCA)
fviz_pca_biplot(uinImpPCA)

#Designate purported males and females to be separated into different hulls in PCA

coords <- uinImpPCA$ind$coord
coords <- as.data.frame(coords)
isFem <- rownames(coords) %in% uinFem
coords$Fem <- isFem
coords$sex <- uinSex$sex
coords$sex <- as.factor(coords$sex)

hull_sex <- coords %>%
  group_by(sex) %>%
  slice(chull(Dim.1, Dim.2))

hull_sex$sex <- as.factor(hull_sex$sex)
hull_sex$Sex <- hull_sex$sex

coords$Sex <- as.factor(coords$sex)

MyPlot <- ggplot(coords, aes(Dim.1, Dim.2)) + geom_point(shape = 21)
MyPlot + aes(fill = factor(sex)) + geom_polygon(data = hull_sex, alpha = 0.5) +
  labs(x = "PC1", y = "PC2", title = "PCA of Uintatherium Skulls")+ geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme(legend.title = element_blank()) +
  geom_label_repel(aes(label = rownames(coords))) +
  scale_fill_manual(values = c("#B8B8B8", "#F8766D", "#619CFF"))

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

#Determine number of clusters data appears to form

gapStat <- fviz_nbclust(uinImp[,1:9], kmeans, k.max = 8, method = "gap_stat")
fviz_nbclust(uinImp[,1:9], kmeans, k.max = 8, method = "gap_stat")

gapStat <- fviz_nbclust(uinImp[,1:9], kmeans, k.max = 17, method = "gap_stat")
wss <- fviz_nbclust(uinImp[,1:9], kmeans, k.max = 17, method = "wss")
silh <- fviz_nbclust(uinImp[,1:9], kmeans, k.max = 17, method = "silhouette")

wss / silh / gapStat

#Use k-means clustering to form two most favorable clusters

km.UinTrimImp <- kmeans(uinImp[,1:9], centers = 2)
UintatheresClustered <- uinImp
clustVec <- km.UinTrimImp$cluster
UintatheresClustered <- as.data.frame(UintatheresClustered)
UintatheresClustered$cluster = clustVec
UintatheresClustered$cluster <- as.factor(UintatheresClustered$cluster)

ClusteredPCA <- PCA(UintatheresClustered[,1:9])

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

hull_sex$sex <- as.factor(hull_sex$sex)
hull_sex$Sex <- hull_sex$sex

###Plot sex and cluster on the same PCA

MyClusterPlot +
  geom_polygon(
    data = hull_cluster,
    aes(x = Dim.1, y = Dim.2, group = Cluster, fill = factor(Cluster)),
    alpha = 0.25,
    linewidth = 0.3,
    colour = "gray25",
    inherit.aes = FALSE) +
  scale_fill_manual(name = "Cluster", values = c("1" = "#4DB6AC", "2" = "#D4A72C")) +
  new_scale_fill() +
  geom_point(aes(fill = Sex), shape = 21, size = 3) +
  geom_label_repel(aes(label = rownames(clusterCoords), fill = Sex),colour = "black" ) +
  scale_fill_manual(values = c("?" = "#B8B8B8","F" = "#F8766D","M" = "#619CFF")) +
  labs(x = "PC1 (60.05%)", y = "PC2 (19.61%)", title = "PCA of Uintathere Skulls") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)

###########################RUNNING MANOVA TO SEE IF PC SCORE IS AFFECTED BY SEX

sigCols <- uinImpPCA$ind$coord[,1:4]
sigCols <- as.data.frame(sigCols)
sigCols$sex <- NA
sigCols$sex <- ifelse(rownames(sigCols) %in% uinFem, "F", 
                      ifelse(rownames(sigCols) %in% uinMale, "M", "?"))
uinPCAMan <- manova(cbind(Dim.1, Dim.2, Dim.3, Dim.4) ~ sex, data = sigCols)
summary(uinPCAMan)
eta_squared(uinPCAMan, partial = TRUE)

################### ESTIMATING DEGREE OF DIMORPHISM IN EACH 
################### METRIC THROUGH FINITE MIXTURE MODELING AND METHOD OF MEANS

estimateTable <- matrix(nrow = 9, ncol = 2)
rownames(estimateTable) <- colnames(Uintatheres)
colnames(estimateTable) <- c("FMA", "MM")
estimateTable <- as.data.frame(estimateTable)

pearsonTable <- read.csv("pearsonTable.csv")

for (i in 1:nrow(estimateTable)) {
  metric <- na.omit(Uintatheres[,i])
  mean <- mean(metric)
  OR <- max(metric) - min(metric)
  n <- length(metric)
  k <- pearsonTable$stDevsInRange[pearsonTable$SampleSize == n]
  sigmaTot <- OR/k
  expSub <- k*(2^.5) #expected number of subpop st devs in total OR
  twoSubPercent <- 2/expSub #portion of OR expected to be occupied by 2 sigmaSub
  oneSubPercent <- twoSubPercent/2 #portion of OR expected to be occupied by 1 sigmaSub
  distanceToSub <- oneSubPercent*OR #numeric value of one sigmaSub <- distance from overall mean to mean of each subpop
  bigMean <- mean + distanceToSub
  smallMean <- mean - distanceToSub
  fmaEstimate <- bigMean/smallMean
  estimateTable$FMA[i] <- fmaEstimate
  bigs <- metric[metric > mean(metric)]
  bigMean <- mean(bigs)
  smalls <- metric[metric < mean(metric)]
  smallMean <- mean(smalls)
  mmEstimate <- bigMean/smallMean
  estimateTable$MM[i] <- mmEstimate
}
