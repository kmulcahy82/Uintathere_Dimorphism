library(FSA)
library(propagate)
library(minpack.lm)
library(nlstools)

#setwd()

set.seed(123)

uintatheres <- read.csv("Uintatherium_Ages.csv", check.names = F) #load dataframe

uintatheres <- uintatheres[!is.na(uintatheres$`Age`), ]

uintatheres$year <- NA
rownames(uintatheres) <- uintatheres$Specimen
uintatheres <- uintatheres[uintatheres$Specimen != "YPM VP 011039L",] #remove damaged view of YPM 11039
uinFem <- uintatheres$Specimen[uintatheres$sex=="F"]
uinMale <- uintatheres$Specimen[uintatheres$sex=="M"]
uinNB <- uintatheres$Specimen[uintatheres$sex=="?"]
uintatheres <- uintatheres[,-c(2,3)] #remove helper variables
#uintatheres <- uintatheres[,-8] #remove maxCan for now

for (i in 1:nrow(uintatheres)){
  if (uintatheres$Age[i] == 1){
    uintatheres$year[i] <- runif(1, 1, 5)
  }
  else if (uintatheres$Age[i] == 2){
    uintatheres$year[i] <- runif(1, 5, 8)
  }
  else if (uintatheres$Age[i] == 3){
    uintatheres$year[i] <- runif(1, 8, 15)
  }
  else if (uintatheres$Age[i] == 4){
    uintatheres$year[i] <- runif(1, 15, 30)
  }
}

uintatheres["PM 55406A","year"] <- 2 #set age of PM 55406A to 2, necessary for subsequent analyses

babyUin <- matrix(nrow = 1, ncol = ncol(uintatheres))
babyUin <- data.frame(babyUin)
colnames(babyUin) <- colnames(uintatheres)
babyUin <- babyUin[,c(1:10)]
babyUin[1,] <- c("PM 3896 Estimate", 36.8056, 30.304, 35.37337, 11.43149, 16.5132, 
                 12.10645, 12.10645, 0, 4.230811)

#####METRIC 1: Nasals to Occipital Crest

set.seed(123)

nasOct <- uintatheres[!is.na(uintatheres$nasOcc),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$nasOcc <- runif(1, 
                                            min = 0.75*as.numeric(babyUin$nasOcc), 
                                            max = 0.95*as.numeric(uintatheres["PM 8019",]$nasOcc))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  nasOct <- rbind(nasOct, baby)
}

plot(nasOct$year, nasOct$nasOcc, xlim = c(0,30), ylim = c(30,90))
#fit Gompertz curve to whole data
x<-as.vector(nasOct$year)
y<-as.vector(nasOct$nasOcc)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), 
                      data=as.list(nasOct), control = nls.lm.control(maxiter = 100))
gompertz.uin
BIC(gompertz.uin)
logistic.uin <- nlsLM(y~l/(1+2.71828^(q+(k*x))), start=list(q=5,k=0.5,l=75), 
                      data=as.list(nasOct), control = nls.lm.control(maxiter = 100))
logistic.uin
BIC(logistic.uin)
label_colors <- ifelse(nasOct$Specimen %in% uinFem, "red", "black")
points(x=seq(from=0, to=38, by=0.1), y=(71.3438)*2.71828^(-2.71828^(-(0.3635)*(seq(from=0, to=38, by=0.1)-(-0.5555)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)
BIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(nasOct$Specimen %in% uinFem, "red",
                       ifelse(nasOct$Specimen %in% uinMale, "blue", "black"))
plot(nasOct$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(nasOct$year, uin.resid, labels = nasOct$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(nasOct, uin.resid>0)
uin.female<-subset(nasOct, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$nasOcc, col="blue", xlab="Age estimate (years)", 
     ylab="Length of skull from nasals to occipital crest", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(30,90), pch=16)
points(uin.female$year, uin.female$nasOcc, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$nasOcc, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$nasOcc, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$nasOcc)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(73.3228)*2.71828^(-2.71828^(-(0.6380)*(seq(from=0, to=40, by=0.1)-(-0.9228)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best
BIC(gompertz.uin.m)

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$nasOcc)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(66.8876)*2.71828^(-2.71828^(-(0.3874)*(seq(from=0, to=40, by=0.1)-(-0.1702)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

nasOctEffectSize <- ((73.7992/66.8876)-1)*100
nasOctNum <- length(nasOct$nasOcc)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(nasOctEffectSize, 2), 
                     "\nN = ", nasOctNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

nasOccPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 2: Nasals to Occipital Condyles

set.seed(123)

nasCon <- uintatheres[!is.na(uintatheres$nasCon),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$nasCon <- runif(1, 
                       min = 0.75*as.numeric(babyUin$nasCon), 
                       max = 0.95*as.numeric(uintatheres["PM 8019",]$nasCon))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  nasCon <- rbind(nasCon, baby)
}

plot(nasCon$year, nasCon$nasCon, xlim = c(0,30), ylim = c(25,85))
#fit Gompertz curve to whole data
x<-as.vector(nasCon$year)
y<-as.vector(nasCon$nasCon)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(nasCon), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(67.8390)*2.71828^(-2.71828^(-(0.3883)*(seq(from=0, to=38, by=0.1)-(0.1112)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(nasCon$Specimen %in% uinFem, "red",
                       ifelse(nasCon$Specimen %in% uinMale, "blue", "black"))
plot(nasCon$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(nasCon$year, uin.resid, labels = nasCon$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(nasCon, uin.resid>0)
uin.female<-subset(nasCon, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$nasCon, col="blue", xlab="Age estimate (years)", 
     ylab="Length of skull from nasals to occipital condyles", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(25,85), pch=16)
points(uin.female$year, uin.female$nasCon, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$nasCon, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$nasCon, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$nasCon)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(70.4213)*2.71828^(-2.71828^(-(0.4708)*(seq(from=0, to=40, by=0.1)-(-0.3404)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$nasCon)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(63.6550)*2.71828^(-2.71828^(-(0.4151)*(seq(from=0, to=40, by=0.1)-(0.4188)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

nasConEffectSize <- ((70.4213/63.6550)-1)*100
nasConNum <- length(nasCon$nasCon)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(nasConEffectSize, 2), 
                     "\nN = ", nasConNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

nasOccPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 3: Premaxillaries to Occipital Condyles

set.seed(123)

preCon <- uintatheres[!is.na(uintatheres$preCon),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$preCon <- runif(1, 
                       min = 0.75*as.numeric(babyUin$preCon), 
                       max = 0.95*as.numeric(uintatheres["PM 8019",]$preCon))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  preCon <- rbind(preCon, baby)
}

plot(preCon$year, preCon$preCon, xlim = c(0,30), ylim = c(25,85))
#fit Gompertz curve to whole data
x<-as.vector(preCon$year)
y<-as.vector(preCon$preCon)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(preCon), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(63.3582)*2.71828^(-2.71828^(-(0.3684)*(seq(from=0, to=38, by=0.1)-(-0.8570)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(preCon$Specimen %in% uinFem, "red",
                       ifelse(preCon$Specimen %in% uinMale, "blue", "black"))
plot(preCon$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(preCon$year, uin.resid, labels = preCon$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(preCon, uin.resid>0)
uin.female<-subset(preCon, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$preCon, col="blue", xlab="Age estimate (years)", 
     ylab="Length of skull from premaxillaries to occipital condyles", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(25,85), pch=16)
points(uin.female$year, uin.female$preCon, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$preCon, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$preCon, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$preCon)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(66.8721)*2.71828^(-2.71828^(-(0.4139)*(seq(from=0, to=40, by=0.1)-(-1.1864)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$preCon)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(58.8128)*2.71828^(-2.71828^(-(0.4608)*(seq(from=0, to=40, by=0.1)-(-0.3227)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

preConEffectSize <- ((66.8721/58.8128)-1)*100
preConNum <- length(preCon$preCon)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(preConEffectSize, 2), 
                     "\nN = ", preConNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

nasOccPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 4: Occipital Condyles to Occipital Crests

set.seed(123)

occCon <- uintatheres[!is.na(uintatheres$occCon),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$occCon <- runif(1, 
                       min = 0.75*as.numeric(babyUin$occCon), 
                       max = 0.95*as.numeric(uintatheres["PM 8019",]$occCon))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  occCon <- rbind(occCon, baby)
}

plot(occCon$year, occCon$occCon, xlim = c(0,30), ylim = c(5,35))
#fit Gompertz curve to whole data
x<-as.vector(occCon$year)
y<-as.vector(occCon$occCon)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(occCon), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(25.9570)*2.71828^(-2.71828^(-(0.4007)*(seq(from=0, to=38, by=0.1)-(0.1715)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(occCon$Specimen %in% uinFem, "red",
                       ifelse(occCon$Specimen %in% uinMale, "blue", "black"))
plot(occCon$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(occCon$year, uin.resid, labels = occCon$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(occCon, uin.resid>0)
uin.female<-subset(occCon, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$occCon, col="blue", 
     xlab="Age estimate (years)", 
     ylab="Height of skull from occipital condyles to occipital crest", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(5,35), pch=16)
points(uin.female$year, uin.female$occCon, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$occCon, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$occCon, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$occCon)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(27.7404)*2.71828^(-2.71828^(-(0.5397)*(seq(from=0, to=40, by=0.1)-(-0.1624)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$occCon)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(23.3023)*2.71828^(-2.71828^(-(0.4153)*(seq(from=0, to=40, by=0.1)-(0.3255)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

occConEffectSize <- ((27.7404/23.3023)-1)*100
occConNum <- length(occCon$occCon)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(occConEffectSize, 2), 
                     "\nN = ", occConNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

occConPlot <- recordPlot()

rm(pred)
gc()

###########################


#####METRIC 5: PARIETAL PROTUBERANCE

set.seed(123)

parPro <- uintatheres[!is.na(uintatheres$parPro),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$parPro <- runif(1, 
                       min = 0.75*as.numeric(babyUin$parPro), 
                       max = 0.95*as.numeric(uintatheres["PM 8019",]$parPro))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  parPro <- rbind(parPro, baby)
}

plot(parPro$year, parPro$parPro, xlim = c(0,30), ylim = c(10,50))
#fit Gompertz curve to whole data
x<-as.vector(parPro$year)
y<-as.vector(parPro$parPro)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(parPro), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(38.4781)*2.71828^(-2.71828^(-(0.3526)*(seq(from=0, to=38, by=0.1)-(0.2085)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(parPro$Specimen %in% uinFem, "red",
                       ifelse(parPro$Specimen %in% uinMale, "blue", "black"))
plot(parPro$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(parPro$year, uin.resid, labels = parPro$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(parPro, uin.resid>0)
uin.female<-subset(parPro, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$parPro, col="blue", 
     xlab="Age estimate (years)", 
     ylab="Height of parietal protuberance", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(10,50), pch=16)
points(uin.female$year, uin.female$parPro, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$parPro, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$parPro, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$parPro)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(41.1247)*2.71828^(-2.71828^(-(0.5142)*(seq(from=0, to=40, by=0.1)-(-0.1179)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$parPro)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(35.51208)*2.71828^(-2.71828^(-(0.23110)*(seq(from=0, to=40, by=0.1)-(0.08935)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

parProEffectSize <- ((41.1247/35.51208)-1)*100
parProNum <- length(parPro$parPro)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(parProEffectSize, 2), 
                     "\nN = ", parProNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

parProPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 6: MAXILLARY PROTUBERANCE

set.seed(123)

maxAlv <- uintatheres[!is.na(uintatheres$maxAlv),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$maxAlv <- runif(1, 
                       min = 0.75*as.numeric(babyUin$maxAlv), 
                       max = 0.95*as.numeric(uintatheres["PM 55406A",]$maxAlv))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  maxAlv <- rbind(maxAlv, baby)
}

plot(maxAlv$year, maxAlv$maxAlv, xlim = c(0,30), ylim = c(5,35))
#fit Gompertz curve to whole data
x<-as.vector(maxAlv$year)
y<-as.vector(maxAlv$maxAlv)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=75), data=as.list(maxAlv), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(24.22246)*2.71828^(-2.71828^(-(0.28175)*(seq(from=0, to=38, by=0.1)-(0.02104)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(maxAlv$Specimen %in% uinFem, "red",
                       ifelse(maxAlv$Specimen %in% uinMale, "blue", "black"))
plot(maxAlv$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(maxAlv$year, uin.resid, labels = maxAlv$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(maxAlv, uin.resid>0)
uin.female<-subset(maxAlv, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$maxAlv, col="blue", 
     xlab="Age estimate (years)", 
     ylab="Maxillary Protuberance Height", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(5,35), pch=16)
points(uin.female$year, uin.female$maxAlv, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$maxAlv, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$maxAlv, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)


#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$maxAlv)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(27.785185)*2.71828^(-2.71828^(-(0.343833)*(seq(from=0, to=40, by=0.1)-(-0.001614)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$maxAlv)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(21.0328)*2.71828^(-2.71828^(-(0.2253)*(seq(from=0, to=40, by=0.1)-(-0.2824)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

maxAlvEffectSize <- ((27.785185/21.0328)-1)*100
maxAlvNum <- length(maxAlv$maxAlv)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(maxAlvEffectSize, 2), 
                     "\nN = ", maxAlvNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

maxAlvPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 7: DISTANCE FROM MAXILLARY PROTUBERANCE TO TIP OF CANINE

set.seed(123)

maxCan <- uintatheres[!is.na(uintatheres$maxCan),]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$maxCan <- runif(1, 
                       min = 0.75*as.numeric(babyUin$maxCan), 
                       max = 0.95*as.numeric(uintatheres["PM 55406A",]$maxCan))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  maxCan <- rbind(maxCan, baby)
}

plot(maxCan$year, maxCan$maxCan, xlim = c(0,30), ylim = c(5,55))
#fit Gompertz curve to whole data
x<-as.vector(maxCan$year)
y<-as.vector(maxCan$maxCan)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=25), data=as.list(maxCan), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(43.1805)*2.71828^(-2.71828^(-(0.3632)*(seq(from=0, to=38, by=0.1)-(1.3124)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(maxCan$Specimen %in% uinFem, "red",
                       ifelse(maxCan$Specimen %in% uinMale, "blue", "black"))
plot(maxCan$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(maxCan$year, uin.resid, labels = maxCan$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(maxCan, uin.resid>0)
uin.female<-subset(maxCan, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$maxCan, col="blue", 
     xlab="Age estimate (years)", 
     ylab="Height from top of maxillary protuberance to canine tip", 
     cex.lab=1.25, xlim=c(0,40),ylim=c(5,55), pch=16)
points(uin.female$year, uin.female$maxCan, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$maxCan, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$maxCan, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)

#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$maxCan)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(46.8926)*2.71828^(-2.71828^(-(0.4747)*(seq(from=0, to=40, by=0.1)-(0.7092)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$maxCan)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(35.3959)*2.71828^(-2.71828^(-(0.2249)*(seq(from=0, to=40, by=0.1)-(1.6786)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

maxCanEffectSize <- ((46.8926/35.3959)-1)*100
maxCanNum <- length(maxCan$maxCan)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(maxCanEffectSize, 2), 
                     "\nN = ", maxCanNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

maxCanPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 8: Length of Canine

set.seed(123)

alvCan <- uintatheres[!is.na(uintatheres$alvCan),]
#alvCan <- uintatheres[uintatheres$Specimen != "AMNH 1671",]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$alvCan <- runif(1, 
                       min = 0, 
                       max = 0.95*as.numeric(uintatheres["PM 55406A",]$alvCan))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  alvCan <- rbind(alvCan, baby)
}

plot(alvCan$year, alvCan$alvCan, xlim = c(0,30), ylim = c(0,30))
#fit Gompertz curve to whole data
x<-as.vector(alvCan$year)
y<-as.vector(alvCan$alvCan)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=25), data=as.list(alvCan), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(18.895)*2.71828^(-2.71828^(-(0.626)*(seq(from=0, to=38, by=0.1)-(1.959)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(alvCan$Specimen %in% uinFem, "red",
                       ifelse(alvCan$Specimen %in% uinMale, "blue", "black"))
plot(alvCan$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(alvCan$year, uin.resid, labels = alvCan$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(alvCan, uin.resid>0)
uin.female<-subset(alvCan, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$alvCan, col="blue", 
     xlab="Age estimate (years)", 
     ylab="Canine Length", 
     cex.lab=1.25, xlim=c(0,30),ylim=c(0,30), pch=16)
points(uin.female$year, uin.female$alvCan, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$alvCan, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$alvCan, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)

#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$alvCan)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(22.5807)*2.71828^(-2.71828^(-(0.5054)*(seq(from=0, to=40, by=0.1)-(1.7833)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$alvCan)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(15.9033)*2.71828^(-2.71828^(-(0.7882)*(seq(from=0, to=40, by=0.1)-(1.9496)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

alvCanEffectSize <- ((22.5807/15.9033)-1)*100
alvCanNum <- length(alvCan$alvCan)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(alvCanEffectSize, 2), 
                     "\nN = ", alvCanNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

alvCanPlot <- recordPlot()

rm(pred)
gc()

###########################

#####METRIC 9: Diastemal Length

set.seed(123)

diaLen <- uintatheres[!is.na(uintatheres$diaLen),]
#alvCan <- uintatheres[uintatheres$Specimen != "AMNH 1671",]
for (i in 1:5){
  baby <- matrix(nrow = 1, ncol = ncol(uintatheres))
  baby <- data.frame(baby)
  colnames(baby) <- colnames(uintatheres)
  baby$diaLen <- runif(1, 
                       min = 0.75*as.numeric(babyUin$diaLen), 
                       max = 0.95*as.numeric(uintatheres["PM 8019",]$diaLen))
  baby$year <- runif(1, min = 0, max = 1)
  rownames(baby) <- paste("baby", i)
  diaLen <- rbind(diaLen, baby)
}

plot(diaLen$year, diaLen$diaLen, xlim = c(0,30), ylim = c(0,15))
#fit Gompertz curve to whole data
x<-as.vector(diaLen$year)
y<-as.vector(diaLen$diaLen)

gompertz.uin <- nlsLM(y~l*2.71828^(-2.71828^(-k*(x-i))), start=list(i=5,k=0.5,l=25), data=as.list(diaLen), control = nls.lm.control(maxiter = 100))
gompertz.uin
points(x=seq(from=0, to=38, by=0.1), y=(7.4113)*2.71828^(-2.71828^(-(0.4607)*(seq(from=0, to=38, by=0.1)-(-0.6521)))), type="l", col="RED", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin)

# Residuals
uin.resid<-residuals(gompertz.uin)
label_colors <- ifelse(diaLen$Specimen %in% uinFem, "red",
                       ifelse(diaLen$Specimen %in% uinMale, "blue", "black"))
plot(diaLen$year, uin.resid, pch=19,cex=.75,xlab="Age estimate (years)", ylab="Residuals",cex.lab=1.25) # plot residuals
text(diaLen$year, uin.resid, labels = diaLen$Specimen, pos = 4,                   # position: 1=below, 2=left, 3=above, 4=right
     cex = 0.7, col = label_colors) 
abline(h=0)

# Positve residuals are male, negative are female
uin.male<-subset(diaLen, uin.resid>0)
uin.female<-subset(diaLen, uin.resid<0)
#jpeg(filename=paste0("maia.gompertz.jpg"), height=200, width=300, res=300, units="mm") #create empty jpeg object
plot(uin.male$year, uin.male$diaLen, col="blue", 
     xlab="Age estimate (years)", 
     ylab="Diastemal Length", 
     cex.lab=1.25, xlim=c(0,30),ylim=c(0,15), pch=16)
points(uin.female$year, uin.female$diaLen, col="red", pch=16)
male.label.colors <- ifelse(uin.male$Specimen %in% uinFem, "orange", 
                            ifelse(uin.male$Specimen %in% uinMale, "blue", "black"))
text(uin.male$year, uin.male$diaLen, 
     labels = uin.male$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = male.label.colors)
female.label.colors <- ifelse(uin.female$Specimen %in% uinMale, "purple", 
                              ifelse(uin.female$Specimen %in% uinFem, "red", "black"))
text(uin.female$year, uin.female$diaLen, 
     labels = uin.female$Specimen, 
     pos = 3, 
     cex = 0.6, 
     col = female.label.colors)

#Fit Gompertz curve to males
xm<-as.vector(uin.male$year)
ym<-as.vector(uin.male$diaLen)
gompertz.uin.m<-nlsLM(ym~l*2.71828^(-2.71828^(-k*(xm-i))), start=list(i=2,k=0.3,l=38), data=as.list(uin.male), control = nls.lm.control(maxiter = 500))
gompertz.uin.m
points(x=seq(from=0, to=40, by=0.1), y=(8.3605)*2.71828^(-2.71828^(-(0.3313)*(seq(from=0, to=40, by=0.1)-(-1.3026)))), type="l", col="blue", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.m) #Really should be using AIC here to test if Gompertz, Logistic, von Bertalanffy, logarithmic, exponential, linear  etc. fits best

#Fit Gompertz curve to females
xf<-as.vector(uin.female$year)
yf<-as.vector(uin.female$diaLen)
gompertz.uin.f<-nlsLM(yf~l*2.71828^(-2.71828^(-k*(xf-i))), start=list(i=2,k=0.3,l=20), data=as.list(uin.female))
gompertz.uin.f
points(x=seq(from=0, to=40, by=0.1), y=(6.4699)*2.71828^(-2.71828^(-(0.7049)*(seq(from=0, to=40, by=0.1)-(-0.1199)))), type="l", col="red", lwd=3) #Insert the parameters found in the Gompertz model manually
AIC(gompertz.uin.f)

diaLenEffectSize <- ((8.3605/6.4699)-1)*100
diaLenNum <- length(diaLen$diaLen)

#CI & PI male
xfit <- with(uin.male, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=2) 
lines(xfit, hi, col = "blue", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.m,data.frame(xm = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.m,data.frame(xm = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "blue", lwd=3, lty=3) 
lines(xfit, hi, col = "blue", lwd=3, lty=3) # the hi and low confidence

rm(pred)
gc()

#CI & PI female
xfit <- with(uin.female, seq(0,40,by=0.25)) # create x values for predict() and predictNLS() "newdata", 100 values for a plot is usually sufficient
pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "confidence", nsim = 10000) # Get 95% CIs
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=2) 
lines(xfit, hi, col = "red", lwd=3, lty=2) # the hi and low confidence band

rm(pred)
gc()

pred<-predictNLS(gompertz.uin.f,data.frame(xf = xfit), level = 0.95,interval = "prediction", nsim = 10000) # Get 95% Prediction intervals
yfit <- predict(gompertz.uin.f,data.frame(xf = xfit)) # get fitted y values 
low <- pred$summary$`Prop.2.5%` # Nb sub `Prop.2.5%` to `Sim.2.5%` to get monte carlo simulation values
hi <- pred$summary$`Prop.97.5%`
lines(xfit, low, col = "red", lwd=3, lty=3) 
lines(xfit, hi, col = "red", lwd=3, lty=3) # the hi and low confidence

text(x = par("usr")[1] + 1, 
     y = par("usr")[4] - 2, 
     labels = paste0("Es = ", round(diaLenEffectSize, 2), 
                     "\nN = ", diaLenNum),
     adj = c(0, 1),     # left-align horizontally, top-align vertically
     cex = 0.9,         # text size
     font = 2)          # bold

diaLenPlot <- recordPlot()

rm(pred)
gc()

metrics <- c("nasOct", "nasCon", "preCon", "occCon",
             "parPro", "maxAlv", "alvCan", "diaLen")

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
