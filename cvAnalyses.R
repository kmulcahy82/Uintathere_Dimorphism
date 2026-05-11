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
library(EnvStats)
library(pwr)

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
bisonAdults$sex <- adultSex

bisonAdults <- bisonAdults[bisonAdults$sex != "?",]
bisonMale <- bisonAdults[bisonAdults$sex == "M",]
bisonFem <- bisonAdults[bisonAdults$sex == "F",]

sexVar <- data.frame(1, 1:8)
rownames(sexVar) <- colnames(bisonAdults[1:8])
colnames(sexVar) <- c("SDI", "CV")

for (i in 1:8) {
  SDI <- mean(bisonMale[,i][!is.na(bisonMale[,i])])/
    mean(bisonFem[,i][!is.na(bisonFem[,i])])
  sexVar[i,1] <- SDI
  allSpecs <- c(bisonMale[,i][!is.na(bisonMale[,i])], 
                bisonFem[,i][!is.na(bisonFem[,i])])
  myCV <- cv(allSpecs)
  sexVar[i,2] <- myCV
}

bisonLM <- lm(sexVar$SDI ~ sexVar$CV)

Uintatheres <- read.csv("Uintathere_Measurements_With_Turnbull.csv", check.names = F) #load dataframe
Uintatheres <- Uintatheres[Uintatheres$Specimen != "YPM 11039L",] #remove damaged view of YPM 11039
rownames(Uintatheres) <- Uintatheres[,1] #rename rows as specimen numbers
uinFem <- Uintatheres$Specimen[Uintatheres$sex=="F"]
uinMale <- Uintatheres$Specimen[Uintatheres$sex=="M"]
uinNB <- Uintatheres$Specimen[Uintatheres$sex=="?"]
Uintatheres <- Uintatheres[,4:12] #remove helper variables
Uintatheres <- Uintatheres[,-7] #remove old variable
#Uintatheres[1,]$alvCan <- NA #This canine is potentially pathologic so it's inclusion is questionable

uintaSDI <- data.frame(1, 1:8)
rownames(uintaSDI) <- colnames(Uintatheres)
colnames(uintaSDI) <- c("CV", "SDI")

for (i in 1:8) {
  myCV <- cv(Uintatheres[,i][!is.na(Uintatheres[,i])])
  uintaSDI[i,1] <- myCV
  mySDI <- 1.8441*myCV + 0.9441
  uintaSDI[i,2] <- mySDI
}

ggplot(sexVar, aes(x = CV, y = SDI)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
  coord_cartesian(xlim = c(0.05, 0.275)) +
  geom_point(data = uintaSDI, color = "red", size = 2.5) +
  geom_text_repel(data = uintaSDI, aes(label = rownames(uintaSDI)), 
                  segment.linetype = "dashed",
                  nudge_x = -0.01,
                  nudge_y = 0.1,
                  size = 5) + 
  geom_text_repel(data = sexVar, aes(label = rownames(sexVar)), 
                  segment.linetype = "dashed", size = 4) +
  theme(axis.title = element_text(size = 18))

summary(bisonLM)
