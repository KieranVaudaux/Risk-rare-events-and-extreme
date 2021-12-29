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

## Naïve approach Oil dataset

qu.min <- quantile(data$Oil, 0.5)
qu.max <- quantile(data$Oil,(length(data$Oil)-30)/length(data$Oil))
mrlplot(data$Oil,c(qu.min, qu.max))
par(mfrow=c(1,2))
tcplot(data$Oil,c(qu.min, qu.max))

th <- 2.5

fit_Oil <- fpot(data$Oil,2.5,np=252)
par(mfrow=c(1,1))
plot(fit_Oil)

par(mfrow=c(1,2))
plot(profile(fit_Oil))

fit2<-gpd.fit(data$Oil,threshold=2.5, npy=252)
gpd.diag(fit2)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=100,xlow=12,xup=30,npy=252,conf = 0.95)

rl100 <- th + fit_Oil$est[[1]]/fit_Oil$est[[2]]*((252*fit_Oil$pat)^fit_Oil$est[[2]]-1)
rl100

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

exi(X_30$Oil,th_30$Oil,4)
XX <- data.frame(row.names = rownames(X_30))
XX["obs"] <- X_30$Oil
XX["time"] <- rownames(XX)
exiplot(XX,c(2,5),r=2)

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

model_list <- list(c(),1,2,3,c(1,2),c(1,3),c(2,3),c(1,2,3))
model_name <- list("stationary","1","2","3","1-2","1-3","2-3","1-2-3")
model_selection_30 <- data.frame(row.names = model_name)
model_selection_90 <- data.frame(row.names = model_name)
model_selection_180 <- data.frame(row.names = model_name)
model_selection_365<- data.frame(row.names = model_name)

fit_30_Oil_list <- list()
fit_90_Oil_list <- list()
fit_180_Oil_list <- list()
fit_365_Oil_list <- list()

for(i in 1:8){
  fit_30_Oil <- gpd.fit(X_30$Oil[2:length(X_30$Oil)],th_30$Oil[2:length(X_30$Oil)],npy = 252,
                        ydat = as.matrix(scale(cbind(X_30$Oil[1:length(X_30$Oil)-1],th_30$Oil[1:length(X_30$Oil)-1],
                                                     Z_30_Oil[1:length(X_30$Oil)-1],ZZ_30_Oil[1:length(X_30$Oil)-1]))), 
                        sigl = model_list[[i]],siglink = exp)
  model_selection_30[model_name[[i]],"nllh"] <- fit_30_Oil$nllh[[1]]
  model_selection_30[model_name[[i]],"shape"] <- fit_30_Oil$mle[length(fit_30_Oil$mle)]
  model_selection_30[model_name[[i]],"shape_standard_error"] <- fit_30_Oil$se[length(fit_30_Oil$se)]
  append(fit_30_Oil_list,fit_30_Oil)
  
  fit_90_Oil <- gpd.fit(X_90$Oil[2:length(X_90$Oil)],th_90$Oil[2:length(X_90$Oil)],npy = 252, 
                        ydat = as.matrix(scale(cbind(X_90$Oil[1:length(X_90$Oil)-1],th_90$Oil[1:length(X_90$Oil)-1],
                                                     Z_90_Oil[1:length(X_90$Oil)-1],ZZ_90_Oil[1:length(X_90$Oil)-1]))), 
                        sigl = model_list[[i]],siglink = exp)
  model_selection_90[model_name[[i]],"nllh"] <- fit_90_Oil$nllh
  model_selection_90[model_name[[i]],"shape"] <- fit_90_Oil$mle[length(fit_90_Oil$mle)]
  model_selection_90[model_name[[i]],"shape_standard_error"] <- fit_90_Oil$se[length(fit_90_Oil$se)]
  append(fit_90_Oil_list,fit_90_Oil)
  
  fit_180_Oil <- gpd.fit(X_180$Oil[2:length(X_180$Oil)],th_180$Oil[2:length(X_180$Oil)],npy = 252, 
                         ydat = as.matrix(scale(cbind(X_180$Oil[1:length(X_180$Oil)-1],th_180$Oil[1:length(X_180$Oil)-1],
                                                      Z_180_Oil[1:length(X_180$Oil)-1],ZZ_180_Oil[1:length(X_180$Oil)-1]))), 
                         sigl = model_list[[i]],siglink = exp)
  model_selection_180[model_name[[i]],"nllh"] <- fit_180_Oil$nllh
  model_selection_180[model_name[[i]],"shape"] <- fit_180_Oil$mle[length(fit_180_Oil$mle)]
  model_selection_180[model_name[[i]],"shape_standard_error"] <- fit_180_Oil$se[length(fit_180_Oil$se)]
  append(fit_180_Oil_list,fit_180_Oil)
  
  fit_365_Oil <- gpd.fit(X_365$Oil[2:length(X_365$Oil)],th_365$Oil[2:length(X_365$Oil)],npy = 252, 
                         ydat = as.matrix(scale(cbind(X_365$Oil[1:length(X_365$Oil)-1],th_365$Oil[1:length(X_365$Oil)-1],
                                                      Z_365_Oil[1:length(X_365$Oil)-1],ZZ_365_Oil[1:length(X_365$Oil)-1]))), 
                         sigl = model_list[[i]],siglink = exp)
  model_selection_365[model_name[[i]],"nllh"] <- fit_365_Oil$nllh
  model_selection_365[model_name[[i]],"shape"] <- fit_365_Oil$mle[length(fit_365_Oil$mle)]
  model_selection_365[model_name[[i]],"shape_standard_error"] <- fit_365_Oil$se[length(fit_365_Oil$se)]
  append(fit_365_Oil_list,fit_365_Oil)
  
  png(file = paste("Figures/POT_Oil/fit_30_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_30_Oil)
  dev.off()
  
  png(file = paste("Figures/POT_Oil/fit_90_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_90_Oil)
  dev.off()
  
  png(file = paste("Figures/POT_Oil/fit_180_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_180_Oil)
  dev.off()
  
  png(file = paste("Figures/POT_Oil/fit_365_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_365_Oil)
  dev.off()
}

nb_param <- c(0,1,1,1,2,2,2,3)
model_selection_30["AIC"] <- 2*nb_param + 2*model_selection_30$nllh
model_selection_90["AIC"] <- 2*nb_param + 2*model_selection_90$nllh
model_selection_180["AIC"] <- 2*nb_param + 2*model_selection_180$nllh
model_selection_365["AIC"] <- 2*nb_param + 2*model_selection_365$nllh

fit_90_Oil <- gpd.fit(X_90$Oil[2:length(X_90$Oil)],th_90$Oil[2:length(X_90$Oil)],npy = 252,
                      ydat = as.matrix(scale(cbind(X_90$Oil[1:length(X_90$Oil)-1],th_90$Oil[1:length(X_90$Oil)-1],
                                                   Z_90_Oil[1:length(X_90$Oil)-1],ZZ_90_Oil[1:length(X_90$Oil)-1]))), 
                      sigl = c(1,2,3),siglink = exp)
sigma <- exp(fit_90_Oil$mle[1] + fit_90_Oil$mle[2]*X_90$Oil[1:length(X_90$Oil)-1]  + fit_90_Oil$mle[3]*th_90$Oil[1:length(X_90$Oil)-1]
             + fit_90_Oil$mle[4]*Z_90_Oil[1:length(X_90$Oil)-1])
YY <- Y_90$Oil[2:length(X_90$Oil)]/sigma

fpot(YY,threshold = 0,npp = 252,scale = 1)
model_selection_90
              

