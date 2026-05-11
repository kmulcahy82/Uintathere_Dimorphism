library(ggplot2)
library(dplyr)

byHand <- read.csv("Uintathere_By_Hand_Meausurements.csv")
imageJ <- read.csv("Uintathere_ImageJ_Measurements.csv", check.names = F) #load dataframe

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

for (i in 1:length(metNames)) {
  df <- data.frame(
    specimen = rownames(byHandSorted),
    hand = byHandClean[,i],
    imageJ = imageJshort[,i]
  )
  
  df$hand <- as.numeric(df$hand)
  df$imageJ <- as.numeric(df$imageJ)
  df <- na.omit(df)
  
  df <- df %>%
    mutate(
      difference = hand - imageJ,
      percentBias = difference / hand*100,
      average = (hand+imageJ) / 2
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
      x = "ImageJ Measurement (mm)",
      y = "Hand Measurement (mm)"
    )
  
  BAPlot <- ggplot(df, aes(x = average, y = difference)) +
    geom_point(size = 3) +
    geom_hline(yintercept = meanDiff, color = "red", linetype = "dashed") +
    geom_hline(yintercept = meanDiff + 1.96 * sd(df$difference), linetype = "dotted") +
    geom_hline(yintercept = meanDiff - 1.96 * sd(df$difference), linetype = "dotted") +
    theme_minimal() +
    labs(
      title = paste(metNames[i], "Bland–Altman Plot"),
      x = "Average of Hand and ImageJ (mm)",
      y = "Difference (Hand - ImageJ) (mm)"
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

nasOccTest <- cor.test(bigDF$`nasOcc avg`, bigDF$`nasOcc dif`)
nasConTest <- cor.test(bigDF$`nasCon avg`, bigDF$`nasCon dif`)
preConTest <- cor.test(bigDF$`preCon avg`, bigDF$`preCon dif`)
occConTest <- cor.test(bigDF$`occCon avg`, bigDF$`occCon dif`)
parProTest <- cor.test(bigDF$`parPro avg`, bigDF$`parPro dif`)
maxAlvTest <- cor.test(bigDF$`maxAlv avg`, bigDF$`maxAlv dif`)
alvCanTest <- cor.test(bigDF$`alvCan avg`, bigDF$`alvCan dif`)
diaLenTest <- cor.test(bigDF$`diaLen avg`, bigDF$`diaLen dif`)

usedMets <- metNames[-7]

stats <- data.frame(
  metric = character(),
  statistic = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

for (i in usedMets) {
  testObj <- get(paste0(i, "Test"))  # e.g., "nasOccTest"
  
  stats <- rbind(
    stats,
    data.frame(
      metric = i,
      statistic = unname(testObj$statistic),
      p.value = testObj$p.value
    )
  )
}

write.csv(stats, "pearsonTests.csv")
stats
