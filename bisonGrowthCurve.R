library(FSA)
library(propagate)
library(minpack.lm)
library(nlstools)

#setwd()

set.seed(123)

bison <- read.csv("Bison_Ages.csv", check.names = F) #load dataframe
rownames(bison) <- bison[,1] #rename rows as specimen numbers
bison$Specimen <- as.character(bison$Specimen)

for (i in 5:14) {
  bison[,i] <- as.numeric(bison[,i])
}

bison <- bison[!is.na(bison$`AgeGroup`), ]

bison$year <- NA

bisIndet <- bison$Specimen[bison$Sex=="?"]
bisMale <- bison$Specimen[bison$Sex=="M"]
bisFem <- bison$Specimen[bison$Sex=="F"]

for (i in 1:nrow(bison)){
  if (bison$AgeGroup[i] == 1){
    bison$year[i] <- runif(1, 1.5, 4)
  }
  else if (bison$AgeGroup[i] == 2){
    bison$year[i] <- runif(1, 4, 6)
  }
  else if (bison$AgeGroup[i] == 3){
    bison$year[i] <- runif(1, 6, 12)
  }
  else if (bison$AgeGroup[i] == 4){
    bison$year[i] <- runif(1, 12, 25)
  }
}

babyBis <- matrix(nrow = 1, ncol = ncol(bison))
babyBis <- data.frame(babyBis)
colnames(babyBis) <- colnames(bison)
babyBis <- babyBis[,c(1, 5:13)]
babyBis[1,] <- c("Newborn Estimate", 
                 1.05*(mean(bison$ectEct)*0.6), #newborn tipFro appears to be about 5% longer than newborn ectEct
                 0.80*(mean(bison$ectEct)*0.6), #newborn basFro appears to be about 0.8x newborn ectEct
                 mean(bison$ectEct)*0.6, #newborn ectEct appears to be about 0.6x adult ectEct
                 mean(na.omit(bison$maxMax))*0.6, #newborn maxMax appears to be about 0.6x adult maxMax
                 mean(na.omit(bison$proNas))*0.6, #newborn proNas appears to be about 0.6x adult proNas
                 mean(bison$nasNuc)*0.6, #newborn nasNuc appears to be about 0.6x adult nasNuc
                 0.15*(mean(bison$ectEct)*0.6), #newborn basTip appears to be about 0.15x newborn ectEct
                 mean(bison$ectEct)*0.6, #newborn ectPro appears to be about equivalent to newborn ectEct
                 0.12*(mean(bison$ectEct)*0.6)) #newborn basWid appears to be about 0.12x newborn ectEct


set.seed(123)

tipFro <- bison[!is.na(bison$tipFro),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$tipFro <- runif(1, 
                                            min = 0.75 * as.numeric(babyBis$tipFro), #minimum set as newborn basFro
                                            max = 1.25 * as.numeric(babyBis$tipFro)) #maximum set as 1.4x basFro
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  tipFro <- rbind(tipFro, baby)
}

plot(tipFro$year, tipFro$tipFro, xlim = c(0,25), ylim = c(10,75))
#fit Gompertz curve to whole data
x<-as.vector(tipFro$year)
y<-as.vector(tipFro$tipFro)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(tipFro), control = nls.lm.control(maxiter = 100))
gompertz.bis
BIC(gompertz.bis)
logistic.bis <- nlsLM(y~l/(1+2.71828^(q+(k*x))), start=list(q=5,k=0.5,l=75), 
                      data=as.list(tipFro), control = nls.lm.control(maxiter = 100))
logistic.bis
BIC(logistic.bis)


label_colors <- ifelse(tipFro$Specimen %in% bisFem, "red",
                       ifelse(tipFro$Specimen %in% bisMale, "blue", "black"))
text(tipFro$year, tipFro$tipFro, 
     labels = tipFro$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = label_colors)
                       
points(x=seq(from=0, to=38, by=0.1), y=(51.7211)*2.71828^(-2.71828^(-(0.3832)*(seq(from=0, to=38, by=0.1)-(0.8459)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(tipFro$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(tipFro$year, bis.resid, labels = tipFro$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(tipFro, bis.resid>0)
bis.female<-subset(tipFro, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$tipFro, col="blue", xlab="Age estimate (years)", 
     ylab="Width Between Tips of Frontal Horns", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(10,75), pch=16)
points(bis.female$year, bis.female$tipFro, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$tipFro, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$tipFro, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$tipFro)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(63.2000)*2.71828^(-2.71828^(-(0.3259)*(seq(from=0, to=40, by=0.1)-(0.8964)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$tipFro)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(42.9689)*2.71828^(-2.71828^(-(0.4898)*(seq(from=0, to=40, by=0.1)-(0.8125)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

tipFroEffectSize <- ((63.2000/42.9689)-1)*100
tipFroNum <- length(tipFro$tipFro)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(tipFroEffectSize, 2), 
                     "\nN = ", tipFroNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

tipFroPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 2: Width between widest margins of orbits

set.seed(123)

ectEct <- bison[!is.na(bison$ectEct),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$ectEct <- runif(1, 
                                                      min = 0.75*as.numeric(babyBis$ectEct), #minumum set as 0.75 estimate of newborn ectEct
                                                      max = 1.25*as.numeric(babyBis$ectEct)) #maximum set as 1.25 estimate of newborn ectEct
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  ectEct <- rbind(ectEct, baby)
}

plot(ectEct$year, ectEct$ectEct, xlim = c(0,25), ylim = c(5,45))
#fit Gompertz curve to whole data
x<-as.vector(ectEct$year)
y<-as.vector(ectEct$ectEct)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(ectEct), control = nls.lm.control(maxiter = 100))
gompertz.bis
BIC(gompertz.bis)
logistic.bis <- nlsLM(y~l/(1+2.71828^(q+(k*x))), start=list(q=5,k=0.5,l=75), 
                      data=as.list(ectEct), control = nls.lm.control(maxiter = 100))
logistic.bis
BIC(logistic.bis)
label_colors <- ifelse(ectEct$Specimen %in% bisFem, "red",
                       ifelse(ectEct$Specimen %in% bisMale, "blue", "black"))
text(ectEct$year, ectEct$ectEct, 
     labels = ectEct$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = label_colors)
points(x=seq(from=0, to=38, by=0.1), y=(28.5183)*2.71828^(-2.71828^(-(0.3053)*(seq(from=0, to=38, by=0.1)-(-0.9467)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(ectEct$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(ectEct$year, bis.resid, labels = ectEct$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(ectEct, bis.resid>0)
bis.female<-subset(ectEct, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$ectEct, col="blue", xlab="Age estimate (years)", 
     ylab="Width Between Lateral Margins of Orbits", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(5,45), pch=16)
points(bis.female$year, bis.female$ectEct, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$ectEct, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$ectEct, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$ectEct)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(31.9014)*2.71828^(-2.71828^(-(0.2132)*(seq(from=0, to=40, by=0.1)-(-2.2787)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$ectEct)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(24.8478)*2.71828^(-2.71828^(-(0.5005)*(seq(from=0, to=40, by=0.1)-(-0.5889)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

ectEctEffectSize <- ((31.9014/24.8478)-1)*100
ectEctNum <- length(ectEct$ectEct)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(ectEctEffectSize, 2), 
                     "\nN = ", ectEctNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

ectEctPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 3: Width between bases of frontal horns

set.seed(123)

basFro <- bison[!is.na(bison$basFro),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$basFro <- runif(1, 
                                                                min = 0.75*as.numeric(babyBis$basFro), #min set at 0.75x estimate of skull width between frontal horns
                                                                max = 1.25*as.numeric(babyBis$basFro)) #max set at 1.25x estimate of skull width between frontal horns
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  basFro <- rbind(basFro, baby)
}

plot(basFro$year, basFro$basFro, xlim = c(0,25), ylim = c(5,45))
#fit Gompertz curve to whole data
x<-as.vector(basFro$year)
y<-as.vector(basFro$basFro)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(basFro), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(basFro$Specimen %in% bisFem, "red",
                       ifelse(basFro$Specimen %in% bisMale, "blue", "black"))
text(basFro$year, basFro$basFro, 
     labels = basFro$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = label_colors)
points(x=seq(from=0, to=38, by=0.1), y=(23.98166)*2.71828^(-2.71828^(-(0.47977)*(seq(from=0, to=38, by=0.1)-(-0.06429)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(basFro$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(basFro$year, bis.resid, labels = basFro$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(basFro, bis.resid>0)
bis.female<-subset(basFro, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$basFro, col="blue", xlab="Age estimate (years)", 
     ylab="Width Between Bases of Frontal Horns", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(5,45), pch=16)
points(bis.female$year, bis.female$basFro, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$basFro, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$basFro, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$basFro)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(26.6199)*2.71828^(-2.71828^(-(0.5374)*(seq(from=0, to=40, by=0.1)-(-0.7129)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$basFro)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(20.53662)*2.71828^(-2.71828^(-(0.75115)*(seq(from=0, to=40, by=0.1)-(-0.01767)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

basFroEffectSize <- ((26.6199/20.53662)-1)*100
basFroNum <- length(basFro$basFro)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(basFroEffectSize, 2), 
                     "\nN = ", basFroNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

basFroPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 4: Width between lateral margins of maxillae

set.seed(123)

maxMax <- bison[!is.na(bison$maxMax),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$maxMax <- runif(1, 
                                                                     min = 0.75*as.numeric(babyBis$maxMax), #min set at 0.75x estimate of skull width between lateral margins of maxillae
                                                                     max = 1.25*as.numeric(babyBis$maxMax)) #max set at 1.25x estimate of skull width between lateral margins of maxillae
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  maxMax <- rbind(maxMax, baby)
}

plot(maxMax$year, maxMax$maxMax, xlim = c(0,25), ylim = c(5,20))
#fit Gompertz curve to whole data
x<-as.vector(maxMax$year)
y<-as.vector(maxMax$maxMax)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(maxMax), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(maxMax$Specimen %in% bisFem, "red",
                       ifelse(maxMax$Specimen %in% bisMale, "blue", "black"))
points(x=seq(from=0, to=38, by=0.1), y=(15.1336)*2.71828^(-2.71828^(-(0.3456)*(seq(from=0, to=38, by=0.1)-(-0.7889)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(maxMax$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(maxMax$year, bis.resid, labels = maxMax$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(maxMax, bis.resid>0)
bis.female<-subset(maxMax, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$maxMax, col="blue", xlab="Age estimate (years)", 
     ylab="Width Between Lateral Margins of Maxillae", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(0,25), pch=16)
points(bis.female$year, bis.female$maxMax, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$maxMax, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$maxMax, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$maxMax)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(16.6839)*2.71828^(-2.71828^(-(0.2211)*(seq(from=0, to=40, by=0.1)-(-2.4681)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$maxMax)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(14.0491)*2.71828^(-2.71828^(-(0.2975)*(seq(from=0, to=40, by=0.1)-(-1.2822)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

maxMaxEffectSize <- ((16.5921/13.9509)-1)*100
maxMaxNum <- length(basFro$maxMax)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(maxMaxEffectSize, 2), 
                     "\nN = ", maxMaxNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

maxMaxPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 5: Length from prosthion to nasion

set.seed(123)

proNas <- bison[!is.na(bison$proNas),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$proNas <- runif(1, 
                                                                     min = 0.75*as.numeric(babyBis$proNas), #min set at 0.75x estimate of skull length from prosthion to nasion
                                                                     max = 1.25*as.numeric(babyBis$proNas)) #max set at 1.25x estimate of skull length from prosthion to nasion
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  proNas <- rbind(proNas, baby)
}

plot(proNas$year, proNas$proNas, xlim = c(0,25), ylim = c(5,35))
#fit Gompertz curve to whole data
x<-as.vector(proNas$year)
y<-as.vector(proNas$proNas)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(proNas), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(proNas$Specimen %in% bisFem, "red",
                       ifelse(proNas$Specimen %in% bisMale, "blue", "black"))
points(x=seq(from=0, to=38, by=0.1), y=(26.1531)*2.71828^(-2.71828^(-(0.5431)*(seq(from=0, to=38, by=0.1)-(-0.2745)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(proNas$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(proNas$year, bis.resid, labels = proNas$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(proNas, bis.resid>0)
bis.female<-subset(proNas, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$proNas, col="blue", xlab="Age estimate (years)", 
     ylab="Length from Prosthion to Nasion", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(5,35), pch=16)
points(bis.female$year, bis.female$proNas, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$proNas, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$proNas, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$proNas)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(28.0246)*2.71828^(-2.71828^(-(0.3984)*(seq(from=0, to=40, by=0.1)-(-1.6542)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$proNas)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(24.356021)*2.71828^(-2.71828^(-(0.789233)*(seq(from=0, to=40, by=0.1)-(0.008366)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

proNasEffectSize <- ((28.0246/24.356021)-1)*100
proNasNum <- length(proNas$proNas)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(proNasEffectSize, 2), 
                     "\nN = ", proNasNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

proNasPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 6: Length from Nasion to Nuchal Crest

set.seed(123)

nasNuc <- bison[!is.na(bison$nasNuc),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$nasNuc <- runif(1, 
                                                  min = 0.75*as.numeric(babyBis$nasNuc), #min set at 0.75x estimate of length from nasion to nuchal crest
                                                  max = 1.25*as.numeric(babyBis$nasNuc)) #max set at 1.25x estimate of length from nasion to nuchal crest
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  nasNuc <- rbind(nasNuc, baby)
}

plot(nasNuc$year, nasNuc$nasNuc, xlim = c(0,25), ylim = c(5,30))
#fit Gompertz curve to whole data
x<-as.vector(nasNuc$year)
y<-as.vector(nasNuc$nasNuc)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(nasNuc), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(nasNuc$Specimen %in% bisFem, "red",
                       ifelse(nasNuc$Specimen %in% bisMale, "blue", "black"))
points(x=seq(from=0, to=38, by=0.1), y=(20.0631)*2.71828^(-2.71828^(-(0.2494)*(seq(from=0, to=38, by=0.1)-(-1.2619)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(nasNuc$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(nasNuc$year, bis.resid, labels = nasNuc$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(nasNuc, bis.resid>0)
bis.female<-subset(nasNuc, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$nasNuc, col="blue", xlab="Age estimate (years)", 
     ylab="Length from Nasion to Nuchal Crest", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(5,30), pch=16)
points(bis.female$year, bis.female$nasNuc, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$nasNuc, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$nasNuc, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$nasNuc)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(23.5322)*2.71828^(-2.71828^(-(0.1607)*(seq(from=0, to=40, by=0.1)-(-2.3578)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$nasNuc)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(18.8928)*2.71828^(-2.71828^(-(0.1793)*(seq(from=0, to=40, by=0.1)-(-2.2608)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

nasNucEffectSize <- ((23.5322/18.8928)-1)*100
nasNucNum <- length(nasNuc$nasNuc)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(nasNucEffectSize, 2), 
                     "\nN = ", nasNucNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

proNasPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 7: Length of horn from base to tip

set.seed(123)

basTip <- bison[!is.na(bison$basTip),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$basTip <- runif(1, 
                                                     min = 0.75*as.numeric(babyBis$basTip), #min set at 0.75x estimate of length of frontal horn 
                                                     max = 1.25*as.numeric(babyBis$basTip)) #max set at 1.25x estimate of length of frontal horn
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  basTip <- rbind(basTip, baby)
}

plot(basTip$year, basTip$basTip, xlim = c(0,25), ylim = c(0,30))
#fit Gompertz curve to whole data
x<-as.vector(basTip$year)
y<-as.vector(basTip$basTip)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(basTip), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(basTip$Specimen %in% bisFem, "red",
                       ifelse(basTip$Specimen %in% bisMale, "blue", "black"))
points(x=seq(from=0, to=38, by=0.1), y=(18.0566)*2.71828^(-2.71828^(-(0.6004)*(seq(from=0, to=38, by=0.1)-(2.3701)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(basTip$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(basTip$year, bis.resid, labels = basTip$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(basTip, bis.resid>0)
bis.female<-subset(basTip, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$basTip, col="blue", xlab="Age estimate (years)", 
     ylab="Length of Horn from Base to Tip", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(0,30), pch=16)
points(bis.female$year, bis.female$basTip, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$basTip, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$basTip, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$basTip)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(21.2437)*2.71828^(-2.71828^(-(0.4685)*(seq(from=0, to=40, by=0.1)-(1.8912)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$basTip)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(14.3975)*2.71828^(-2.71828^(-(0.8238)*(seq(from=0, to=40, by=0.1)-(1.5970)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

basTipEffectSize <- ((21.2437/14.3975)-1)*100
basTipNum <- length(basTip$basTip)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(basTipEffectSize, 2), 
                     "\nN = ", basTipNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

basTipPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 8: Distance from posterolateral corner of orbit to prosthion

set.seed(123)

ectPro <- bison[!is.na(bison$ectPro),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$ectPro <- runif(1, 
                                                                                         min = 0.75*as.numeric(babyBis$ectPro), #min set at 0.75x estimate of length of tip to orbit 
                                                                                         max = 1.25*as.numeric(babyBis$ectPro)) #max set at 1.25x estimate of length of tip to obit
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  ectPro <- rbind(ectPro, baby)
}

plot(ectPro$year, ectPro$ectPro, xlim = c(0,25), ylim = c(5,45))
#fit Gompertz curve to whole data
x<-as.vector(ectPro$year)
y<-as.vector(ectPro$ectPro)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(ectPro), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(ectPro$Specimen %in% bisFem, "red",
                       ifelse(ectPro$Specimen %in% bisMale, "blue", "black"))
points(x=seq(from=0, to=38, by=0.1), y=(33.0744)*2.71828^(-2.71828^(-(0.5847)*(seq(from=0, to=38, by=0.1)-(0.0866)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(ectPro$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(ectPro$year, bis.resid, labels = ectPro$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(ectPro, bis.resid>0)
bis.female<-subset(ectPro, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$ectPro, col="blue", xlab="Age estimate (years)", 
     ylab="Length from prosthion to ectorbitale", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(5,45), pch=16)
points(bis.female$year, bis.female$ectPro, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$ectPro, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$ectPro, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$ectPro)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(34.6235)*2.71828^(-2.71828^(-(0.4611)*(seq(from=0, to=40, by=0.1)-(-0.7157)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$ectPro)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(31.4190)*2.71828^(-2.71828^(-(0.7496)*(seq(from=0, to=40, by=0.1)-(0.4177)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

ectProEffectSize <- ((34.6235/31.4190)-1)*100
ectProNum <- length(ectPro$ectPro)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(ectProEffectSize, 2), 
                     "\nN = ", ectProNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

ectTipPlot <- recordPlot()

rm(pred)
gc()

############################

#####METRIC 9: Width of frontal horn at base

set.seed(123)

basWid <- bison[!is.na(bison$basWid),]

for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(bison))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(bison)
  baby$basWid <- runif(1, 
                       min = 0.75*as.numeric(babyBis$basWid), #min set at 0.75x estimate of length of tip to orbit 
                       max = 1.25*as.numeric(babyBis$basWid)) #max set at 1.25x estimate of length of tip to obit
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  basWid <- rbind(basWid, baby)
}

plot(basWid$year, basWid$basWid, xlim = c(0,25), ylim = c(0,12))
#fit Gompertz curve to whole data
x<-as.vector(basWid$year)
y<-as.vector(basWid$basWid)

gompertz.bis <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(basWid), control = nls.lm.control(maxiter = 100))
gompertz.bis
label_colors <- ifelse(basWid$Specimen %in% bisFem, "red",
                       ifelse(basWid$Specimen %in% bisMale, "blue", "black"))
points(x=seq(from=0, to=38, by=0.1), y=(6.3718)*2.71828^(-2.71828^(-(0.3943)*(seq(from=0, to=38, by=0.1)-(1.3780)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis)

# Residuals
bis.resid<-residuals(gompertz.bis)
plot(basWid$year, bis.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(basWid$year, bis.resid, labels = basWid$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
bis.male<-subset(basWid, bis.resid>0)
bis.female<-subset(basWid, bis.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(bis.male$year, bis.male$basWid, col="blue", xlab="Age estimate (years)", 
     ylab="Width of frontal horn at base", 
     cex.lab=1.25, xlim=c(0,25),ylim=c(0,12), pch=16)
points(bis.female$year, bis.female$basWid, col="red", pch=16)
male.label.colors <- ifelse(bis.male$Specimen %in% bisFem, "orange", 
                            ifelse(bis.male$Specimen %in% bisMale, "blue", "black"))
text(bis.male$year, bis.male$basWid, 
     labels = bis.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(bis.female$Specimen %in% bisMale, "purple", 
                              ifelse(bis.female$Specimen %in% bisFem, "red", "black"))
text(bis.female$year, bis.female$basWid, 
     labels = bis.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(bis.male$year)
ym<-as.vector(bis.male$basWid)
gompertz.bis.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(bis.male), control = nls.lm.control(maxiter = 500))
gompertz.bis.m
points(x=seq(from=0, to=40, by=0.1), y=(7.8697)*2.71828^(-2.71828^(-(0.4646)*(seq(from=0, to=40, by=0.1)-(0.7750)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(bis.female$year)
yf<-as.vector(bis.female$basWid)
gompertz.bis.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(bis.female))
gompertz.bis.f
points(x=seq(from=0, to=40, by=0.1), y=(4.6464)*2.71828^(-2.71828^(-(0.9755)*(seq(from=0, to=40, by=0.1)-(0.7798)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.bis.f)

basWidEffectSize <- ((7.8697/4.6464)-1)*100
basWidNum <- length(basWid$basWid)

#CI & PI male
xfit <- with(bis.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(bis.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.bis.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.bis.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(basWidEffectSize, 2), 
                     "\nN = ", basWidNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

ectTipPlot <- recordPlot()

rm(pred)
gc()

#MAKE TABLE WITH METRICS, EFFECT SIZES, AND NUMBERS

metrics <- c("tipFro", "basFro", "ectEct", "proNas",
             "nasNuc", "ectPro", "basTip", "basWid")

metricStats <- data.frame(
  metric = character(8),
  n = numeric(8),
  n. = numeric(8),
  effectSize = numeric(8)
)

metricStats$metric <- metrics

get(paste0(metrics[2],"EffectSize"))

for (i in 1:nrow(metricStats)) {
  metricStats[i,]$n <- get(paste0(metrics[i],"Num"))
  metricStats[i,]$n. <- (get(paste0(metrics[i],"Num")) - 5)
  metricStats[i,]$effectSize <- get(paste0(metrics[i],"EffectSize"))
}

write.csv(metricStats, "bisonStats.csv")
