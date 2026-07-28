setwd("")

library(ggplot2)
library(dplyr)
library(effectsize)
library(cowplot)
library(gridExtra)

test <- read.csv("bisonTest.csv")
test$angle <- as.factor(test$angle)
test$scalebarDepth <- as.factor(test$scalebarDepth)
test$zoom <- as.factor(test$zoom)

notZoomed <- test[test$zoom != "zoomed",] #only testing effects of angle and scalebarDepth

#plotting each variable for all unzoomed to see patterns in data
plot(notZoomed$tipFro ~ notZoomed$basFro)
text(notZoomed$basFro, notZoomed$tipFro, labels = notZoomed$Trial, pos=3, cex=0.5)
plot(notZoomed$tipFro)
text(row_number(notZoomed), notZoomed$tipFro, labels = notZoomed$Trial, pos=3, cex=0.5)
plot(notZoomed$basFro)
text(row_number(notZoomed), notZoomed$basFro, labels = notZoomed$Trial, pos=3, cex=0.5)
plot(notZoomed$ectEct)
text(row_number(notZoomed), notZoomed$ectEct, labels = notZoomed$Trial, pos=3, cex=0.5)

#Alternative plotting technique
#Also use ANOVA to quantify effects of scalebar placement and camera angle

#basFro
basFroPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = basFro, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
basFroTest <- aov(basFro ~ angle*scalebarDepth, data = notZoomed)
summary(basFroTest)
eta_squared(basFroTest, partial = TRUE)

#ectEct
ectEctPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = ectEct, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
ectEctTest <- aov(ectEct~ angle*scalebarDepth, data = notZoomed)
summary(ectEctTest)
eta_squared(ectEctTest, partial = TRUE)

#proNas
proNasPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = proNas, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
proNasTest <- aov(proNas~ angle*scalebarDepth, data = notZoomed)
summary(proNasTest)
eta_squared(proNasTest, partial = TRUE)

#nasNuc
nasNucPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = nasNuc, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
nasNucTest <- aov(nasNuc~ angle*scalebarDepth, data = notZoomed)
summary(nasNucTest)
eta_squared(nasNucTest, partial = TRUE)

#basTip
basTipPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = basTip, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
basTipTest <- aov(basTip~ angle*scalebarDepth, data = notZoomed)
summary(basTipTest)
eta_squared(basTipTest, partial = TRUE)

#ectPro
ectProPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = ectPro, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
ectProTest <- aov(ectPro ~ angle*scalebarDepth, data = notZoomed)
summary(ectProTest)
eta_squared(ectProTest, partial = TRUE)

#basWid
basWidPlot <- ggplot(notZoomed) +
  aes(x = scalebarDepth, y = basWid, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
basWidTest <- aov(basWid ~ angle*scalebarDepth, data = notZoomed)
summary(basWidTest)
eta_squared(basWidTest, partial = TRUE)

legend <- get_legend(basWidPlot)

basFroPlot <- basFroPlot + theme(legend.position = "none")
ectEctPlot <- ectEctPlot + theme(legend.position = "none")
proNasPlot <- proNasPlot + theme(legend.position = "none")
nasNucPlot <- nasNucPlot + theme(legend.position = "none")
basTipPlot <- basTipPlot + theme(legend.position = "none")
ectProPlot <- ectProPlot + theme(legend.position = "none")
basWidPlot <- basWidPlot + theme(legend.position = "none")

grid.arrange(basFroPlot, ectEctPlot, proNasPlot, nasNucPlot,
             basTipPlot, ectProPlot, basWidPlot, legend, nrow = 4, ncol = 2)

###Visualizing and quantifying the effects of camera zoom
zoomTest <- test[test$scalebarDepth == "mid",]

#basFro
basFroZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = basFro, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
basFroZoomTest <- aov(basFro ~ angle*zoom, data = zoomTest)
summary(basFroZoomTest)
eta_squared(basFroZoomTest, partial = TRUE)

#ectEct
ectEctZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = ectEct, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
ectEctZoomTest <- aov(ectEct ~ angle*zoom, data = zoomTest)
summary(ectEctZoomTest)
eta_squared(ectEctZoomTest, partial = TRUE)

#proNas
proNasZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = proNas, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
proNasZoomTest <- aov(proNas ~ angle*zoom, data = zoomTest)
summary(proNasZoomTest)
eta_squared(proNasZoomTest, partial = TRUE)

#nasNuc
nasNucZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = nasNuc, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
nasNucZoomTest <- aov(nasNuc ~ angle*zoom, data = zoomTest)
summary(nasNucZoomTest)
eta_squared(nasNucZoomTest, partial = TRUE)

#basTip
basTipZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = basTip, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
basTipZoomTest <- aov(basTip ~ angle*zoom, data = zoomTest)
summary(basTipZoomTest)
eta_squared(basTipZoomTest, partial = TRUE)

#ectPro
ectProZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = ectPro, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
ectProZoomTest <- aov(ectPro ~ angle*zoom, data = zoomTest)
summary(ectProZoomTest)
eta_squared(ectProZoomTest, partial = TRUE)

#basWid
basWidZoomPlot <- ggplot(zoomTest) +
  aes(x = zoom, y = basWid, fill = angle) +
  geom_boxplot() + 
  scale_fill_manual(values = c(
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73"   # green
  ))
basWidZoomTest <- aov(basWid ~ angle*zoom, data = zoomTest)
summary(basWidZoomTest)
eta_squared(basWidZoomTest, partial = TRUE)

legend <- get_legend(basWidZoomPlot)

basFroZoomPlot <- basFroZoomPlot + theme(legend.position = "none")
ectEctZoomPlot <- ectEctZoomPlot + theme(legend.position = "none")
proNasZoomPlot <- proNasZoomPlot + theme(legend.position = "none")
nasNucZoomPlot <- nasNucZoomPlot + theme(legend.position = "none")
basTipZoomPlot <- basTipZoomPlot + theme(legend.position = "none")
ectProZoomPlot <- ectProZoomPlot + theme(legend.position = "none")
basWidZoomPlot <- basWidZoomPlot + theme(legend.position = "none")

grid.arrange(basFroZoomPlot, ectEctZoomPlot, proNasZoomPlot, nasNucZoomPlot,
             basTipZoomPlot, ectProZoomPlot, basWidZoomPlot, legend, nrow = 4, ncol = 2)

### Calculating maximum variance within treatments

angles <- unique(test$angle)
depths <- unique(test$scalebarDepth)
zooms <- unique(test$zoom)
metrics <- colnames(test)[2:9]

ranges <- matrix(nrow=144, ncol = 5)
colnames(ranges) <- c("treatment", "min", "max", "mean", "maxVar")
ranges <- as.data.frame(ranges)

allTreats <- vector(mode = "character", length = 144)
treat <- 1
for (column in 1:length(metrics)) {
  for (i in 1:length(angles)) {
    for (j in 1:length(depths)){
      for (k in 1: length(zooms)){
        values <- paste(metrics[column], angles[i], depths[j], zooms[k])
        print(values)
        allTreats[treat] <- values
        treat <- treat + 1
      }
    }
  }
}

ranges$treatment <- allTreats
rows <- 1

for (met in 1:length(metrics)) {
  for (i in 1:length(angles)) {
    for (j in 1:length(depths)){
      for (k in 1: length(zooms)){
        column <- test[,colnames(test) == metrics[met]]
        values <- column[test$angle == levels(angles)[i] 
                          & test$scalebarDepth == levels(depths)[j]  
                          & test$zoom == levels(zooms)[k]]
        if (is.numeric(values)) {
          ranges[rows,]$min <- min(values)
          ranges[rows,]$max <- max(values)
          ranges[rows,]$mean <- mean(values)
          ranges[rows,]$maxVar <- ((max(values)/min(values))*100)-100
        }
        rows <- rows+1
      }
    }
  }
}

max(ranges$maxVar[!is.na(ranges$maxVar)]) ###Most a measurement using the same zoom, angle, and scalebar depth can vary by
mean(ranges$maxVar[!is.na(ranges$maxVar)])

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

### Testing if observed differences w/in treatment affect bison and uintathere results

#########Load and manipulate uintathere data

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

set.seed(123)

Uintatheres <- Uintatheres * runif(length(unlist(Uintatheres)), 0.865, 1.135)

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

#############
### Same analysis with Bison now

#########Load and manipulate data

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

bisonSex <- bisonAdults
bisonSex$sex <- adultSex

colnames(bisonAdults)

set.seed(123)

bisonAdults <- bisonAdults * runif(length(unlist(bisonAdults)), 0.865, 1.135)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PART 1: Generate histograms for each variable and test
# for normality and unimodality

#Setting up "Stats" dataframe to store values of various test statistics of each measurement

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

#The 'width  between tips of frontal horns' measurement is not considered in final analyses

Stats <- data.frame(Measurement, n, wStatistic, pValue, dipStatistic, pValue)

#RUN for-loop that plots histograms for each variable, runs Shapiro-Wilk and 
#Hartigan's Dip Test, and stores these variables in "Stats" dataframe.
#Plot all 7 histograms together

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

grid.arrange(dims[[1]], dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]],
             ncol = 3, nrow = 3)

#generating random normal data for reference

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

#plot with random normal data

grid.arrange(dims[[2]], dims[[3]], dims[[4]], dims[[5]], dims[[6]], dims[[7]], dims[[8]], random,
             ncol = 2, nrow = 4)

#Mixture modeling

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

# write.csv(Stats, "bisonStats.csv")
# write.csv(mixtureModels, "bisonMixtureModels.csv")

par(mfrow=c(1,1))

#PART 2: Multivariate analysis

#First, missing data must be imputed for PCA. 
######WORKING WITH PROPPER IMPUTATION PROTOCOLS (Impute if 3 or fewer missing, otherwise, remove specimen)

bisonAdults <- bisonAdults[,-1] #remove extra horn variable

bis.ncp <- estim_ncpPCA(bisonAdults, method.cv= "Kfold")
plot(bis.ncp$criterion~names(bis.ncp$criterion),xlab="number of dimensions", 
     ylab = "ncp criterion")
ncp <- as.numeric(bis.ncp$ncp[1]) #The optimal number of dimensions(components) which should considered when imputing data, as obtained by k-fold cross validation

bis.impute <- imputePCA(bisonAdults, ncp = ncp)
bis.impute <- bis.impute$completeObs
bis.impute <- as.data.frame(bis.impute)

bis.impute #Missing values are now imputed according to the PCA cross-validation.

####Hierarchical Clustering

bis.impute$sex <- adultSex

bis.impute$sex <- as.factor(bis.impute$sex)

bisNum <- bis.impute[,1:7]

bisScale <- scale(bisNum) #need to scale data before running PCA

bisScale <- as.data.frame(bisScale)

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

###ALTERNATIVE CLUSTERING WITH BOOTSTRAP SUPPORT VALUES

bootMap <- pvclust(t(bisScale), method.dist="euclidean", method.hclust="complete", nboot=1000)
par(mar=c(1.25,2,2,2))
plot(bootMap,labels = FALSE, hang = -1, main = "Hierarchical Clustering of Bison Skulls", print.num = F)
tipLabels <- bootMap$hclust$labels[bootMap$hclust$order]
bisFem <- rownames(bisonSex)[bisonSex$sex=="F"]
bisMale <- rownames(bisonSex)[bisonSex$sex=="M"]
bisNB <- rownames(bisonSex)[bisonSex$sex=="?"]
tipColors <- ifelse(bootMap$hclust$labels[bootMap$hclust$order] %in% bisFem, "red",
                    ifelse(bootMap$hclust$labels[bootMap$hclust$order] %in% bisMale, "blue", "black"))
text(x = 1:length(tipLabels), y = -0.02, labels = tipLabels, col = tipColors,
     srt = 90,adj = 1,xpd = TRUE,cex = 0.75)

###MANOVA AND R2 FOR HEATMAP
bisScale$sex <- bisonSex$sex

bisHeatMan <- manova(cbind(basFro, ectEct, proNas, nasNuc, basTip, ectPro, basWid) ~ sex, data = bisScale)
summary(bisHeatMan)
summary.aov(bisHeatMan)
eta_squared(bisHeatMan, partial = TRUE)

###PERFORMING PCA

normPCA <- PCA(bisScale[,1:7])

#write.csv(normPCA$eig, "newBisonEigVals.csv")
#write.csv(normPCA$var$coord, "neBisonPCLoad.csv")
# write.csv(normPCA$var$cos2, "bisonPCContVar.csv")
# write.csv(normPCA$var$contrib, "bisonVarContPc.csv")
# 
fviz_pca_ind(normPCA)
fviz_pca_biplot(normPCA)


coords <- normPCA$ind$coord
coords <- as.data.frame(coords)
bisonSex <- bisonAdults
bisonSex$sex <- adultSex
isFem <- bisonSex$sex=="F"
coords$Fem <- isFem
coords$sex <- adultSex
coords$sex <- as.factor(coords$sex)

hull_sex <- coords %>%
  group_by(sex) %>%
  slice(chull(Dim.1, Dim.2))

coords$Sex <- as.factor(coords$sex)

MyPlot <- ggplot(coords, aes(Dim.1, Dim.2)) + geom_point(shape = 21)
MyPlot + aes(fill = factor(sex)) + geom_polygon(data = hull_sex, alpha = 0.5) +
  labs(x = "PC1", y = "PC2", title = "PCA of Bison Skulls")+ geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  theme(legend.title = element_blank()) +
  geom_label_repel(aes(label = rownames(coords))) +
  scale_fill_manual(values = c("#B8B8B8", "#F8766D", "#619CFF"))

# write.csv(normPCA$eig, "eigVals.csv")
# write.csv(normPCA$var$coord, "pcLoad.csv")
# write.csv(normPCA$var$cos2, "pcContVar.csv")
# write.csv(normPCA$var$contrib, "varContPc.csv")

#Determine number of clusters data appears to form

gapStat <- fviz_nbclust(bisScale[,1:7], kmeans, k.max = 22, method = "gap_stat")
wss <- fviz_nbclust(bisScale[,1:7], kmeans, k.max = 22, method = "wss")
silh <- fviz_nbclust(bisScale[,1:7], kmeans, k.max = 22, method = "silhouette")

wss / silh / gapStat

#Use k-means clustering to form two most favorable clusters

kmBison <- kmeans(bisScale[,1:7], centers = 2)
bisClust <- bisScale
clustVec <- kmBison$cluster
bisClust <- as.data.frame(bisClust)
bisClust$cluster = clustVec
bisClust$cluster <- as.factor(bisClust$cluster)

ClusteredPCA <- PCA(bisClust[,1:7])

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
  labs(x = "PC1 (75.50%)", y = "PC2 (16.47%)", title = "PCA of Bison Skulls") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)

###########################RUNNING MANOVA TO SEE IF PC SCORE IS AFFECTED BY SEX

sigCols <- normPCA$ind$coord[,1:3] #only include PC's accounting for first 95% of variance
sigCols <- as.data.frame(sigCols)
sigCols$sex <- NA
bisFem <- rownames(bisonSex)[bisonSex$sex == "F"]
bisMale <- rownames(bisonSex)[bisonSex$sex == "M"]
bisNB <- rownames(bisonSex)[bisonSex$sex == "?"]
sigCols$sex <- ifelse(rownames(sigCols) %in% bisFem, "F", 
                      ifelse(rownames(sigCols) %in% bisMale, "M", "?"))
bisPCAMan <- manova(cbind(Dim.1, Dim.2, Dim.3) ~ sex, data = sigCols)
summary(bisPCAMan)
eta_squared(bisPCAMan, partial = TRUE)

################### ESTIMATING DEGREE OF DIMORPHISM IN EACH 
################### METRIC THROUGH FINITE MIXTURE MODELING AND METHOD OF MEANS

estimateTable <- matrix(nrow = 7, ncol = 2)
rownames(estimateTable) <- colnames(bisonAdults)
colnames(estimateTable) <- c("FMA", "MM")
estimateTable <- as.data.frame(estimateTable)

pearsonTable <- read.csv("pearsonTable.csv")

for (i in 1:nrow(estimateTable)) {
  metric <- na.omit(bisonAdults[,i])
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


