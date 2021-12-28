library('tibble')
library(evd)
library(ismev)
library(stats)
library(nlme)
library(POT)
library(eva)
library(stringr)



Project2021<-read.table(file="Project2021.txt", sep="")
date2021 <- as.Date(rownames(Project2021),"%Y%m%d")

# We consider negative return, so we consider minus the data in order to now look at positive return

data <- - select(Project2021,Oil,Autos)

## Rolling time window 

rolling_threshold <- function(quant, l, dt){
  thresholds <- dt[-(1:l),]
  for(i in 1:(length(dt[,1])-l)){
    thresholds$Oil[i] <- quantile(dt$Oil[i:(i+l)],quant)[[1]]
    thresholds$Autos[i] <- quantile(dt$Autos[i:(i+l)],quant)[[1]]
  }
  thresholds
}

## Rolling thresholds for different rolling time window

qu = 0.95
th_30 <- rolling_threshold(qu,30,data)
th_90 <- rolling_threshold(qu,90,data)
th_180 <- rolling_threshold(qu,180,data)
th_365 <- rolling_threshold(qu,365,data)
plot(th_180$Oil)
## Exceedances and Exceedances times

X_30 <- data[31:length(data[,1]),]
X_90 <- data[91:length(data[,1]),]
X_180 <- data[181:length(data[,1]),]
X_365 <- data[366:length(data[,1]),]

Y_30 <- data[31:length(data[,1]),]-th_30
Y_90 <- data[91:length(data[,1]),]-th_90
Y_180 <- data[181:length(data[,1]),]-th_180
Y_365 <- data[366:length(data[,1]),]-th_365

year_30 <- as.numeric(str_sub(rownames(X_30),1,4))
year_90 <- as.numeric(str_sub(rownames(X_90),1,4))
year_180 <- as.numeric(str_sub(rownames(X_180),1,4))
year_365 <- as.numeric(str_sub(rownames(X_365),1,4))

time_30_Oil <- which(Y_30$Oil>0)
time_30_Autos <- which(Y_30$Autos>0)
time_90_Oil <- which(Y_90$Oil>0)
time_90_Autos <- which(Y_90$Autos>0)
time_180_Oil <- which(Y_180$Oil>0)
time_180_Autos <- which(Y_180$Autos>0)
time_365_Oil <- which(Y_365$Oil>0)
time_365_Autos <- which(Y_365$Autos>0)

Z_30_Oil <- diff(data[30:length(data[,1]),]$Oil)
ZZ_30_Oil <- diff(data[29:length(data[,1]),]$Oil,differences = 2)
Z_90_Oil <- diff(data[90:length(data[,1]),]$Oil)
ZZ_90_Oil <- diff(data[89:length(data[,1]),]$Oil,differences = 2)
Z_180_Oil <- diff(data[180:length(data[,1]),]$Oil)
ZZ_180_Oil <- diff(data[179:length(data[,1]),]$Oil,differences = 2)
Z_365_Oil<- diff(data[365:length(data[,1]),]$Oil)
ZZ_365_Oil <- diff(data[364:length(data[,1]),]$Oil,differences = 2)

## Non-stationarity fit for the Oil 

fit_30_Oil <- gpd.fit(X_30$Oil[2:length(X_30$Oil)],th_30$Oil[2:length(X_30$Oil)],npy = 252,
                      ydat = as.matrix(scale(cbind(X_30$Oil[1:length(X_30$Oil)-1],th_30$Oil[1:length(X_30$Oil)-1],
                                                   Z_30_Oil[1:length(X_30$Oil)-1],ZZ_30_Oil[1:length(X_30$Oil)-1]))), 
                      sigl = c(1,2,3),siglink = exp)
fit_90_Oil <- gpd.fit(X_90$Oil[2:length(X_90$Oil)],th_90$Oil[2:length(X_90$Oil)],npy = 252, 
                      ydat = as.matrix(scale(cbind(X_90$Oil[1:length(X_90$Oil)-1],th_90$Oil[1:length(X_90$Oil)-1],
                                                   Z_90_Oil[1:length(X_90$Oil)-1],ZZ_90_Oil[1:length(X_90$Oil)-1]))), 
                      sigl = c(1,2,3),siglink = exp)
fit_180_Oil <- gpd.fit(X_180$Oil[2:length(X_180$Oil)],th_180$Oil[2:length(X_180$Oil)],npy = 252, 
                       ydat = as.matrix(scale(cbind(X_180$Oil[1:length(X_180$Oil)-1],th_180$Oil[1:length(X_180$Oil)-1],
                                                    Z_180_Oil[1:length(X_180$Oil)-1],ZZ_180_Oil[1:length(X_180$Oil)-1]))), 
                       sigl = c(1,2,3),siglink = exp)
fit_365_Oil <- gpd.fit(X_365$Oil[2:length(X_365$Oil)],th_365$Oil[2:length(X_365$Oil)],npy = 252, 
                       ydat = as.matrix(scale(cbind(X_365$Oil[1:length(X_365$Oil)-1],th_365$Oil[1:length(X_365$Oil)-1],
                                                    Z_365_Oil[1:length(X_365$Oil)-1],ZZ_365_Oil[1:length(X_365$Oil)-1]))), 
                       sigl = c(1,2,3),siglink = exp)

gpd.diag(fit_30_Oil)
gpd.diag(fit_90_Oil)
gpd.diag(fit_180_Oil)
gpd.diag(fit_365_Oil)

sigma = exp(fit_90_Oil$mle[1] + fit_90_Oil$mle[2]*X_90$Oil[1:length(X_90$Oil)-1] + fit_90_Oil$mle[3]*th_90$Oil[1:length(X_90$Oil)-1])
Y_90_Oil <- Y_90$Oil[2:length(Y_90$Oil)]/sigma
sigma


fpot(Y_90_Oil,0)
