setwd("/Users/kevin/Desktop/plosRound3Files/CodeAndData/")

library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)

### Working with uintatheres first

byHand <- read.csv("Uintathere_By_Hand_Meausurements.csv")
imageJ <- read.csv("Uintathere_ImageJ_Measurements.csv", check.names = F) #load dataframe

#trim and manage dataframes to be directly comparable

imageJ <- imageJ[imageJ$Specimen != "YPM VP 011039L",] #remove damaged view of YPM 11039
rownames(imageJ) <- imageJ[,1] #rename rows as specimen numbers
imageJ <- imageJ[,4:12]
rownames(imageJ) <- gsub("[ -]", ".", rownames(imageJ))
rownames(imageJ)[26] <- "YPM.VPPU.010298"
rownames(imageJ)[12] <- "YPM.VPPU.010076"

byHand <- byHand[,-1]
byHandTrans <- t(byHand)
colnames(byHandTrans) <- byHandTrans[1,]
byHandTrans <- byHandTrans[-1,]
byHandTrans <- as.data.frame(byHandTrans)
byHandShort <- byHandTrans[,c(1,3,2,4,18,28,19,24,22)]

imageJshort <- imageJ[rownames(imageJ) %in% rownames(byHandShort),]
byHandShort <- byHandShort[rownames(byHandShort) %in% rownames(imageJshort),]
byHandSorted <- byHandShort[rownames(imageJshort), ]


byHandClean <- byHandSorted
rownames(byHandClean) <- rownames(imageJshort)

byHandClean[4,1] <- NA
byHandClean[12,1] <- NA
byHandClean[15,1] <- 62
byHandClean[1,2] <- 77
byHandClean[2,2] <- NA
byHandClean[3,2] <- 71
byHandClean[4,2] <- NA
byHandClean[6,2] <- NA
byHandClean[10,2] <- NA
byHandClean[11,2] <- 68
byHandClean[12,2] <- NA
byHandClean[13,2] <- NA
byHandClean[14,2] <- 68
byHandClean[15,2] <- 68
byHandClean[2,3] <- NA
byHandClean[4,3] <- NA
byHandClean[6,3] <- NA
byHandClean[7,3] <- NA
byHandClean[8,3] <- NA
byHandClean[10,3] <- NA
byHandClean[11,3] <- 64
byHandClean[12,3] <- NA
byHandClean[13,3] <- NA
byHandClean[14,3] <- 62.5
byHandClean[15,3] <- 57.5
byHandClean[2,4] <- NA
byHandClean[4,4] <- NA
byHandClean[6,4] <- NA
byHandClean[10,4] <- NA
byHandClean[12,4] <- NA
byHandClean[13,4] <- NA
byHandClean[14,4] <- 29
byHandClean[2,5] <- NA
byHandClean[4,5] <- NA
byHandClean[6,5] <- NA
byHandClean[7,5] <- NA
byHandClean[8,5] <- 40.5
byHandClean[10,5] <- NA
byHandClean[12,5] <- NA
byHandClean[13,5] <- NA
byHandClean[6,6] <- NA
byHandClean[7,6] <- 19
byHandClean[8,6] <- 17
byHandClean[13,6] <- NA
byHandClean[2,7] <- NA
byHandClean[4,7] <- 45
byHandClean[5,7] <- NA
byHandClean[6,7] <- NA
byHandClean[7,7] <- NA
byHandClean[8,7] <- NA
byHandClean[9,7] <- NA
byHandClean[10,7] <- NA
byHandClean[13,7] <- NA
byHandClean[15,7] <- NA
byHandClean[2,8] <- NA
byHandClean[4,8] <- 17.5
byHandClean[5,8] <- NA
byHandClean[6,8] <- NA
byHandClean[7,8] <- NA
byHandClean[8,8] <- NA
byHandClean[9,8] <- NA
byHandClean[10,8] <- NA
byHandClean[13,8] <- NA
byHandClean[15,8] <- NA
byHandClean[6,9] <- NA
byHandClean[7,9] <- NA
byHandClean[8,9] <- NA
byHandClean[10,9] <- NA
byHandClean[13,9] <- NA


metNames <- colnames(imageJshort)

plots <- list()

#generate Bland-Altman and Error plots

for (i in 1:length(metNames)) {
  df <- data.frame(
    specimen = rownames(byHandClean),
    hand = byHandClean[,i],
    imageJ = imageJshort[,i]
  )
  
  df$hand <- as.numeric(df$hand)
  df$imageJ <- as.numeric(df$imageJ)
  df <- na.omit(df)
  
  df <- df %>%
    mutate(
      difference = imageJ - hand,
      average = (hand+imageJ) / 2,
      percentBias = difference / average*100
    )
  
  meanBias <- mean(df$percentBias)
  meanDiff <- mean(df$difference)
  correlation <- cor(df$hand, df$imageJ)
  
  cat("Mean percent bias: ", paste(metNames[i]), round(meanBias, 2), "%\n")
  cat("Mean difference: ", paste(metNames[i]), round(meanDiff, 2), "mm\n")
  cat("Correlation: ", paste(metNames[i]), round(correlation, 4), "\n")
  
  diffPlot <- ggplot(df, aes(x = imageJ, y = hand)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste("Hand vs. ImageJ Measurements", metNames[i]),
      x = "ImageJ Measurement (cm)",
      y = "Hand Measurement (cm)"
    )
  
  BAPlot <- ggplot(df, aes(x = average, y = difference)) +
    geom_point(size = 3) +
    geom_hline(yintercept = meanDiff, color = "red", linetype = "dashed") +
    geom_hline(yintercept = meanDiff + 1.96 * sd(df$difference), linetype = "dotted") +
    geom_hline(yintercept = meanDiff - 1.96 * sd(df$difference), linetype = "dotted") +
    theme_minimal() +
    labs(
      title = paste(metNames[i], "Bland–Altman Plot"),
      x = "Average of Hand and ImageJ (cm)",
      y = "Difference (Hand - ImageJ) (cm)"
    )
  
  lm_model <- lm(diff ~ avg, data = data.frame(diff = df$hand - df$imageJ, 
                                               avg = (df$hand + df$imageJ)/2))
  
  correction_factor <- 1 + (meanBias / 100)
  df <- df %>%
    mutate(correctedImageJ = imageJ * correction_factor)
  
  assign(paste(metNames[i]), df)
  assign(paste(metNames[i], "diffPlot"), diffPlot)
  assign(paste(metNames[i], "BAPlot"), BAPlot)
  assign(paste(metNames[i], "linear model"), lm_model)
}

colnames(byHandClean) <- metNames
bigDF <- cbind(byHandClean, imageJshort)

for (i in 1:length(metNames)) {
  avg <- paste(metNames[i], "avg")
  bigDF[,avg] <- NA
}

for (i in 1:length(metNames)) {
  dif <- paste(metNames[i], "dif")
  bigDF[dif] <- NA
}

for (i in 1:ncol(bigDF)) {
  bigDF[,i] <- as.numeric(bigDF[,i])
}

for (j in 1:nrow(bigDF)) {
  for (i in 1:9) {
    if (is.na(bigDF[j,i]) == F && is.na(bigDF[j,i+9]) == F){
      bigDF[j,i+18] <- mean(c(bigDF[j,i], bigDF[j,i+9]))
    }
    if (is.na(bigDF[j,i]) == F && is.na(bigDF[j,i+9]) == F){
      bigDF[j,i+27] <- bigDF[j,i] - bigDF[j,i+9]
    }
  }
}

plot(bigDF$`diaLen avg`, bigDF$`diaLen dif`)

#running linear models in Uintatherium

summary(lm(nasOcc$difference ~ nasOcc$average))
summary(lm(nasOcc$percentBias ~ nasOcc$average))
summary(lm(nasCon$difference ~ nasCon$average))
summary(lm(nasCon$percentBias ~ nasCon$average))
summary(lm(preCon$difference ~ preCon$average))
summary(lm(preCon$percentBias ~ preCon$average))
summary(lm(occCon$difference ~ occCon$average))
summary(lm(occCon$percentBias ~ occCon$average))
summary(lm(parPro$difference ~ parPro$average))
summary(lm(parPro$percentBias ~ parPro$average))
summary(lm(maxAlv$difference ~ maxAlv$average))
summary(lm(maxAlv$percentBias ~ maxAlv$average))
summary(lm(alvCan$difference ~ alvCan$average))
summary(lm(alvCan$percentBias ~ alvCan$average))
summary(lm(diaLen$difference ~ diaLen$average))
summary(lm(diaLen$percentBias ~ diaLen$average))

#######MAKING % DIFFERENCE PLOTS
yMarks <- seq(-40, 40, by = 10)

plot(nasOcc$percentBias ~ nasOcc$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "nasOcc Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(nasOcc$percentBias)
mean(nasOcc$percentBias)
sd(nasOcc$percentBias)

plot(nasCon$percentBias ~ nasCon$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "nasCon Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(nasCon$percentBias)
mean(nasCon$percentBias)
sd(nasCon$percentBias)

plot(preCon$percentBias ~ preCon$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "preCon Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(preCon$percentBias)
mean(preCon$percentBias)
sd(preCon$percentBias)

plot(occCon$percentBias ~ occCon$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "occCon Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(occCon$percentBias)
mean(occCon$percentBias)
sd(occCon$percentBias)

plot(parPro$percentBias ~ parPro$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "parPro Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(parPro$percentBias)
mean(parPro$percentBias)
sd(parPro$percentBias)

plot(maxAlv$percentBias ~ maxAlv$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "maxAlv Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(maxAlv$percentBias)
mean(maxAlv$percentBias)
sd(maxAlv$percentBias)

plot(alvCan$percentBias ~ alvCan$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "alvCan Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(alvCan$percentBias)
mean(alvCan$percentBias)
sd(alvCan$percentBias)

plot(diaLen$percentBias ~ diaLen$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "diaLen Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(diaLen$percentBias)
mean(diaLen$percentBias)
sd(diaLen$percentBias)

grid.arrange(`nasOcc BAPlot`, `nasCon BAPlot`, `preCon BAPlot`, `occCon BAPlot`, `parPro BAPlot`, `maxAlv BAPlot`,
             `alvCan BAPlot`, `diaLen BAPlot`, ncol = 2, nrow = 4)

### Testing if canine length is associated with canine width
uintatheres <- read.csv("Uintathere_Measurements_With_Canine_Width.csv")
plot(uintatheres$canWid ~ uintatheres$alvCan)
canModel <- lm(uintatheres$canWid ~ uintatheres$alvCan)
summary(canModel)

##### WORKING WITH BISON NOW
##### FIRST, SHOW THAT SKULL LENGTH IS CONSISTENT ACROSS PLANES

bisonScans <- read.csv("bisonMorphosource.csv")
bisonScans <- bisonScans[1:5, 1:3]
rownames(bisonScans) <- bisonScans[,1]
bisonScans <- bisonScans[,2:3]
colnames(bisonScans) <- c("Anterior", "Lateral")
bisonScans$Diff <- NA
for (i in 1:nrow(bisonScans)) {
  bisonScans$Diff[i] <- ((bisonScans$Anterior[i]-bisonScans$Lateral[i])/
                            mean(c(bisonScans$Anterior[i], bisonScans$Lateral[i]))) * 100
}
bisonScans$Average <- NA
for (i in 1:nrow(bisonScans)) {
  bisonScans$Average[i] <- mean(c(bisonScans$Anterior[i], bisonScans$Lateral[i]))
}
plot(bisonScans$Diff ~ bisonScans$Average, ylim = c(-40, 40))
range(bisonScans$Diff)
mean(bisonScans$Diff)
sd(bisonScans$Diff)
model <- lm(bisonScans$Anterior ~ bisonScans$Lateral)
summary(model)
dev.off()
plot(bisonScans$Anterior ~ bisonScans$Lateral, xlab = "Lateral measurement (m)", ylab = "Anterior measurement (m)", 
     main= "Bison skull length: Lateral vs Anterior")
abline(model, col = "red", lwd = 2)
rect(0.5, 0.567, 0.520, 0.572, col = "lightblue")
text(0.51, 0.57, expression(R^2 == 0.9148))

##GENERATING ERROR AND BA PLOTS FOR BISON

bisonImageJ <- read.csv("Bison_Measurements.csv")
bisonByHand <- read.csv("bisonByHand.csv")
bisonImageJ <- bisonImageJ[-20,]#remove 68350, specimen with adult dentition still erupting
bisonByHand <- bisonByHand[-20,]#remove 68350, specimen with adult dentition still erupting
rownames(bisonImageJ) <- bisonImageJ[,1]
rownames(bisonByHand) <- bisonByHand[,1]
bisonImageJ <- bisonImageJ[,-c(1:4)]
bisonByHand <- bisonByHand[,-1]
bisonByHand <- bisonByHand[,colnames(bisonByHand) %in% colnames(bisonImageJ)]
bisonImageJ <- bisonImageJ[,colnames(bisonImageJ) %in% colnames(bisonByHand)]
bisonByHand <- bisonByHand %>% relocate(basFro, .after = tipFro)

metNames <- colnames(bisonImageJ)

plots <- list()

for (i in 1:length(metNames)) {
  df <- data.frame(
    specimen = rownames(bisonByHand),
    hand = bisonByHand[,i],
    imageJ = bisonImageJ[,i]
  )
  
  df$hand <- as.numeric(df$hand)
  df$imageJ <- as.numeric(df$imageJ)
  df <- na.omit(df)
  
  df <- df %>%
    mutate(
      difference = imageJ - hand,
      average = (hand+imageJ) / 2,
      percentBias = difference / average*100
    )
  
  meanBias <- mean(df$percentBias)
  meanDiff <- mean(df$difference)
  correlation <- cor(df$hand, df$imageJ)
  
  cat("Mean percent bias: ", paste(metNames[i]), round(meanBias, 2), "%\n")
  cat("Mean difference: ", paste(metNames[i]), round(meanDiff, 2), "mm\n")
  cat("Correlation: ", paste(metNames[i]), round(correlation, 4), "\n")
  
  diffPlot <- ggplot(df, aes(x = imageJ, y = hand)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste("Hand vs. ImageJ Measurements", metNames[i]),
      x = "ImageJ Measurement (cm)",
      y = "Hand Measurement (cm)"
    )
  
  BAPlot <- ggplot(df, aes(x = average, y = difference)) +
    geom_point(size = 3) +
    geom_hline(yintercept = meanDiff, color = "red", linetype = "dashed") +
    geom_hline(yintercept = meanDiff + 1.96 * sd(df$difference), linetype = "dotted") +
    geom_hline(yintercept = meanDiff - 1.96 * sd(df$difference), linetype = "dotted") +
    theme_minimal() +
    labs(
      title = paste(metNames[i], "Bland–Altman Plot"),
      x = "Average of Hand and ImageJ (cm)",
      y = "Difference (Hand - ImageJ) (cm)"
    )
  
  lm_model <- lm(diff ~ avg, data = data.frame(diff = df$hand - df$imageJ, 
                                               avg = (df$hand + df$imageJ)/2))
  
  correction_factor <- 1 + (meanBias / 100)
  df <- df %>%
    mutate(correctedImageJ = imageJ * correction_factor)
  
  assign(paste(metNames[i]), df)
  assign(paste(metNames[i], "diffPlot"), diffPlot)
  assign(paste(metNames[i], "BAPlot"), BAPlot)
  assign(paste(metNames[i], "linear model"), lm_model)
}

for (i in 1:length(metNames)) {
  diffPlot <- paste(metNames[i], "diffPlot")
  baPlot <- paste(metNames[i], "BAPlot")
  print(get(diffPlot))
  print(get(baPlot))
}

for (i in 1:length(metNames)) {
  plotName <- paste(metNames[i], "BAPlot")
  print(get(plotName))
}

bigDF <- cbind(bisonByHand, bisonImageJ)

for (i in 1:length(metNames)) {
  avg <- paste(metNames[i], "avg")
  bigDF[,avg] <- NA
}

for (i in 1:length(metNames)) {
  dif <- paste(metNames[i], "dif")
  bigDF[dif] <- NA
}

for (i in 1:ncol(bigDF)) {
  bigDF[,i] <- as.numeric(bigDF[,i])
}

for (j in 1:nrow(bigDF)) {
  for (i in 1:8) {
    if (is.na(bigDF[j,i]) == F && is.na(bigDF[j,i+8]) == F){
      bigDF[j,i+16] <- mean(c(bigDF[j,i], bigDF[j,i+8]))
    }
    if (is.na(bigDF[j,i]) == F && is.na(bigDF[j,i+8]) == F){
      bigDF[j,i+24] <- bigDF[j,i] - bigDF[j,i+8]
    }
  }
}

#Running linear models in bison

summary(lm(basFro$difference ~ basFro$average))
summary(lm(basFro$percentBias ~ basFro$average))
summary(lm(ectEct$difference ~ ectEct$average))
summary(lm(ectEct$percentBias ~ ectEct$average))
summary(lm(proNas$difference ~ proNas$average))
summary(lm(proNas$percentBias ~ proNas$average))
summary(lm(nasNuc$difference ~ nasNuc$average))
summary(lm(nasNuc$percentBias ~ nasNuc$average))
summary(lm(basTip$difference ~ basTip$average))
summary(lm(basTip$percentBias ~ basTip$average))
summary(lm(ectPro$difference ~ ectPro$average))
summary(lm(ectPro$percentBias ~ ectPro$average))

#######MAKING ERROR PLOTS

yMarks <- seq(from = -40, to = 40, by = 10)

plot(basFro$percentBias ~ basFro$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "basFro Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(basFro$percentBias)
mean(basFro$percentBias)
sd(basFro$percentBias)

plot(ectEct$percentBias ~ ectEct$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "ectEct Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(ectEct$percentBias)
mean(ectEct$percentBias)
sd(ectEct$percentBias)

plot(proNas$percentBias ~ proNas$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "proNas Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(proNas$percentBias)
mean(proNas$percentBias)
sd(proNas$percentBias)

plot(nasNuc$percentBias ~ nasNuc$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "nasNuc Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(nasNuc$percentBias)
mean(nasNuc$percentBias)
sd(nasNuc$percentBias)

plot(basTip$percentBias ~ basTip$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "basTip Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(basTip$percentBias)
mean(basTip$percentBias)
sd(basTip$percentBias)

plot(ectPro$percentBias ~ ectPro$average, ylim = c(-40, 40),
     xlab = "Average of Hand and ImageJ (cm)",
     ylab = "% Difference",
     main = "ectPro Percentage Difference")
axis(side = 2, at = yMarks)
abline(h = yMarks, col = "lightgray", lty = "dotted")
range(ectPro$percentBias)
mean(ectPro$percentBias)
sd(ectPro$percentBias)

