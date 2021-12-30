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

## We consider negative return, so we consider minus the data in order to now look at positive return

data <- - select(Project2021,Oil,Autos)

#### ---- Oil ---- #####


## Naïve approach Oil dataset

qu.min <- quantile(data$Oil, 0.5)
qu.max <- quantile(data$Oil,(length(data$Oil)-30)/length(data$Oil))

# Threshold selection

mrlplot(data$Oil,c(qu.min, qu.max)) # Not so useful
par(mfrow=c(1,2))
tcplot(data$Oil,c(qu.min, qu.max)) # A threshold above 2.5 seems acceptable

th <- 2.5

fit_Oil <- fpot(data$Oil,th,np=252)

par(mfrow=c(2,2))
plot(fit_Oil)

# Normal-based confidence interval for sigma and the shape
low_conf_sigma <- c(fit_Oil$estimate[[1]]-1.96*fit_Oil$std.err[[1]],fit_Oil$estimate[[1]]+1.96*fit_Oil$std.err[[1]])
low_conf_shape <- c(fit_Oil$estimate[[2]]-1.96*fit_Oil$std.err[[2]],fit_Oil$estimate[[2]]+1.96*fit_Oil$std.err[[2]])
low_conf_shape
# The shape is significatively different from zero base on the normal confidence interval

par(mfrow=c(1,2))
plot(profile(fit_Oil))

# We can confirme that the shape parameter is significatively different from zero based on the profile-likelihood

fit2<-gpd.fit(data$Oil,threshold=th, npy=252)
gpd.diag(fit2)

rl1 <- th + fit_Oil$est[[1]]/fit_Oil$est[[2]]*((252*fit_Oil$pat)^fit_Oil$est[[2]]-1)
rl10 <- th + fit_Oil$est[[1]]/fit_Oil$est[[2]]*((2520*fit_Oil$pat)^fit_Oil$est[[2]]-1)
rl100 <- th + fit_Oil$est[[1]]/fit_Oil$est[[2]]*((25200*fit_Oil$pat)^fit_Oil$est[[2]]-1)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=1,xlow=4.1,xup=4.7,npy=252,conf = 0.95)
abline(v=rl1)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=10,xlow=7.5,xup=11,npy=252,conf = 0.95)
abline(v=rl10)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=100,xlow=12.3,xup=28,npy=252,conf = 0.95)
abline(v=rl100)

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

qu = 0.99
th_30 <- rolling_threshold(qu,30,data)
th_90 <- rolling_threshold(qu,90,data)
th_180 <- rolling_threshold(qu,180,data)
th_365 <- rolling_threshold(qu,365,data)

exi(X_30$Oil,th_30$Oil,r=1)
exi(X_90$Oil,th_90$Oil,r=1)
exi(X_180$Oil,th_180$Oil,r=1)
exi(X_365$Oil,th_365$Oil,r=1)

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
ZZ_30_Oil <- data[29:length(data[,1])-2,]$Oil
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
  fit_30_Oil_list$model_name[[i]] <- fit_30_Autos
  
  fit_90_Oil <- gpd.fit(X_90$Oil[2:length(X_90$Oil)],th_90$Oil[2:length(X_90$Oil)],npy = 252, 
                        ydat = as.matrix(scale(cbind(X_90$Oil[1:length(X_90$Oil)-1],th_90$Oil[1:length(X_90$Oil)-1],
                                                     Z_90_Oil[1:length(X_90$Oil)-1],ZZ_90_Oil[1:length(X_90$Oil)-1]))), 
                        sigl = model_list[[i]],siglink = exp)
  model_selection_90[model_name[[i]],"nllh"] <- fit_90_Oil$nllh
  model_selection_90[model_name[[i]],"shape"] <- fit_90_Oil$mle[length(fit_90_Oil$mle)]
  model_selection_90[model_name[[i]],"shape_standard_error"] <- fit_90_Oil$se[length(fit_90_Oil$se)]
  fit_90_Oil_list$model_name[[i]] <- fit_90_Autos
  
  fit_180_Oil <- gpd.fit(X_180$Oil[2:length(X_180$Oil)],th_180$Oil[2:length(X_180$Oil)],npy = 252, 
                         ydat = as.matrix(scale(cbind(X_180$Oil[1:length(X_180$Oil)-1],th_180$Oil[1:length(X_180$Oil)-1],
                                                      Z_180_Oil[1:length(X_180$Oil)-1],ZZ_180_Oil[1:length(X_180$Oil)-1]))), 
                         sigl = model_list[[i]],siglink = exp)
  model_selection_180[model_name[[i]],"nllh"] <- fit_180_Oil$nllh
  model_selection_180[model_name[[i]],"shape"] <- fit_180_Oil$mle[length(fit_180_Oil$mle)]
  model_selection_180[model_name[[i]],"shape_standard_error"] <- fit_180_Oil$se[length(fit_180_Oil$se)]
  fit_180_Oil_list$model_name[[i]] <- fit_180_Autos
  
  fit_365_Oil <- gpd.fit(X_365$Oil[2:length(X_365$Oil)],th_365$Oil[2:length(X_365$Oil)],npy = 252, 
                         ydat = as.matrix(scale(cbind(X_365$Oil[1:length(X_365$Oil)-1],th_365$Oil[1:length(X_365$Oil)-1],
                                                      Z_365_Oil[1:length(X_365$Oil)-1],ZZ_365_Oil[1:length(X_365$Oil)-1]))), 
                         sigl = model_list[[i]],siglink = exp)
  model_selection_365[model_name[[i]],"nllh"] <- fit_365_Oil$nllh
  model_selection_365[model_name[[i]],"shape"] <- fit_365_Oil$mle[length(fit_365_Oil$mle)]
  model_selection_365[model_name[[i]],"shape_standard_error"] <- fit_365_Oil$se[length(fit_365_Oil$se)]
  fit_365_Oil_list$model_name[[i]] <- fit_365_Autos
  
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

# If we use AIC in order to choose which sequences of thresholds use for the rest of the study, we see that the thresholds sequences 
# with a rolling time window of 30 days performe better than the other sequences. Therefore, in the followings of this study we will use 
# this sequence of threshold.

model_selection_30

# The likelihood ratio test does not allow us to reject the null hypothesis which say that the data comes from the 
# model "2" instead of the model "1-2-3", thus as we also have that the AIC is the smallest for the model "2", we choose the 
# parsimonious model "2".
ratio1 <- 2*(model_selection_30$nllh[3]-model_selection_30$nllh[8])
qchisq(0.95,2)
ratio1
p1<-1-pchisq(ratio1,2)
p1


fit_Oil <- gpd.fit(X_30$Oil[2:length(X_30$Oil)],th_30$Oil[2:length(X_30$Oil)],npy = 252,
                  ydat = as.matrix(scale(cbind(th_30$Oil[1:length(X_30$Oil)-1]))), 
                  sigl = c(1),shl=c(1,2,3),siglink = exp,shlink = exp)

gpd.diag(fit_Oil)

# We scale the exceedances in order to have a stationnary sequences of exceedances

Y_rescale <- Y_30$Oil[time_30_Oil]/fit_Oil$vals[,1]
fit_Oil_ <- fpot(Y_rescale,threshold = 0)

# Diagnostics plot of our choosen model
par(mfrow=c(2,2))
plot(fit_Oil_)

# Normal-based confidence interval for the shape parameter
low_conf_shape <- c(fit_Oil_$estimate[[2]]-1.96*fit_Oil_$std.err[[2]],fit_Oil_$estimate[[2]]+1.96*fit_Oil_$std.err[[2]])
low_conf_shape
confint(fit_Oil_)
# The shape is significatively different from zero base on the normal confidence interval

par(mfrow=c(1,2))
plot(profile(fit_Oil_))
# We can confirme that the shape parameter is significatively different from zero based on the profile-likelihood

#### ---- Autos ---- #####

## Naïve approach Autos dataset

qu.min <- quantile(data$Autos, 0.5)
qu.max <- quantile(data$Autos,(length(data$Autos)-30)/length(data$Autos))

# Threshold selection

mrlplot(data$Autos,c(qu.min, qu.max)) # Not so useful
par(mfrow=c(1,2))
tcplot(data$Autos,c(qu.min, qu.max)) # A threshold above 2.5 seems acceptable

quantile(data$Autos, 0.99)
th <- 3.5

fit_Autos <- fpot(data$Autos,th,np=252)
fit_Autos
par(mfrow=c(2,2))
plot(fit_Oil)

# Normal-based confidence interval for sigma and the shape
low_conf_sigma_Autos <- c(fit_Autos$estimate[[1]]-1.96*fit_Autos$std.err[[1]],fit_Autos$estimate[[1]]+1.96*fit_Autos$std.err[[1]])
low_conf_shape_Autos <- c(fit_Autos$estimate[[2]]-1.96*fit_Autos$std.err[[2]],fit_Autos$estimate[[2]]+1.96*fit_Autos$std.err[[2]])
low_conf_shape_Autos
# The shape is not significatively different from zero base on the normal confidence interval

par(mfrow=c(1,2))
plot(profile(fit_Autos))
abline(v=0,col=2,lty=2)

# The profile log-likelihood allow us to reject the hypothesis that the shape parameter is equal to zero

fit2<-gpd.fit(data$Autos,threshold=th, npy=252)
gpd.diag(fit2)

rl1 <- th + fit_Autos$est[[1]]/fit_Autos$est[[2]]*((252*fit_Autos$pat)^fit_Autos$est[[2]]-1)
rl10 <- th + fit_Autos$est[[1]]/fit_Autos$est[[2]]*((2520*fit_Autos$pat)^fit_Autos$est[[2]]-1)
rl100 <- th + fit_Autos$est[[1]]/fit_Autos$est[[2]]*((25200*fit_Autos$pat)^fit_Autos$est[[2]]-1)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=1,xlow=4.6,xup=5.2,npy=252,conf = 0.95)
abline(v=rl1)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=10,xlow=8,xup=11,npy=252,conf = 0.95)
abline(v=rl10)

par(mfrow=c(1,1))
gpd.prof(z=fit2,m=100,xlow=11.8,xup=25,npy=252,conf = 0.95)
abline(v=rl100)

## Rolling time window 

rolling_threshold <- function(quant, l, dt){
  thresholds <- dt[-(1:l),]
  for(i in 1:(length(dt[,1])-l)){
    thresholds$Oil[i] <- quantile(dt$Oil[i:(i+l)],quant)[[1]]
    thresholds$Autos[i] <- quantile(dt$Autos[i:(i+l)],quant)[[1]]
  }
  thresholds
}

##Extremal coefficient

exi(X_30$Oil,th_30$Autos,r=1)
exi(X_90$Oil,th_90$Autos,r=1)
exi(X_180$Oil,th_180$Autos,r=1)
exi(X_365$Oil,th_365$Autos,r=1)

## Covariates 

Z_30_Autos <- diff(data[30:length(data[,1]),]$Autos)
ZZ_30_Autos <- diff(data[29:length(data[,1]),]$Autos,differences = 2)
Z_90_Autos <- diff(data[90:length(data[,1]),]$Autos)
ZZ_90_Autos <- diff(data[89:length(data[,1]),]$Autos,differences = 2)
Z_180_Autos <- diff(data[180:length(data[,1]),]$Autos)
ZZ_180_Autos <- diff(data[179:length(data[,1]),]$Autos,differences = 2)
Z_365_Autos<- diff(data[365:length(data[,1]),]$Autos)
ZZ_365_Autos <- diff(data[364:length(data[,1]),]$Autos,differences = 2)

## Non-stationarity fit for the Oil 

model_list <- list(c(),1,2,3,c(1,2),c(1,3),c(2,3),c(1,2,3))
model_name <- list("stationary","1","2","3","1-2","1-3","2-3","1-2-3")

model_selection_30_Autos <- data.frame(row.names = model_name)
model_selection_90_Autos <- data.frame(row.names = model_name)
model_selection_180_Autos <- data.frame(row.names = model_name)
model_selection_365_Autos <- data.frame(row.names = model_name)

fit_30_Autos_list <- list()
fit_90_Autos_list <- list()
fit_180_Autos_list <- list()
fit_365_Autos_list <- list()

for(i in 1:8){
  fit_30_Autos <- gpd.fit(X_30$Autos[2:length(X_30$Autos)],th_30$Autos[2:length(X_30$Autos)],npy = 252,
                        ydat = as.matrix(scale(cbind(X_30$Autos[1:length(X_30$Autos)-1],th_30$Autos[1:length(X_30$Autos)-1],
                                                     Z_30_Autos[1:length(X_30$Autos)-1],ZZ_30_Autos[1:length(X_30$Autos)-1]))), 
                        sigl = model_list[[i]],siglink = exp)
  model_selection_30_Autos[model_name[[i]],"nllh"] <- fit_30_Autos$nllh[[1]]
  model_selection_30_Autos[model_name[[i]],"shape"] <- fit_30_Autos$mle[length(fit_30_Autos$mle)]
  model_selection_30_Autos[model_name[[i]],"shape_standard_error"] <- fit_30_Autos$se[length(fit_30_Autos$se)]
  fit_30_Autos_list$model_name[[i]] <- fit_30_Autos
  
  fit_90_Autos <- gpd.fit(X_90$Autos[2:length(X_90$Autos)],th_90$Autos[2:length(X_90$Autos)],npy = 252, 
                        ydat = as.matrix(scale(cbind(X_90$Autos[1:length(X_90$Autos)-1],th_90$Autos[1:length(X_90$Autos)-1],
                                                     Z_90_Autos[1:length(X_90$Autos)-1],ZZ_90_Autos[1:length(X_90$Autos)-1]))), 
                        sigl = model_list[[i]],siglink = exp)
  model_selection_90_Autos[model_name[[i]],"nllh"] <- fit_90_Autos$nllh
  model_selection_90_Autos[model_name[[i]],"shape"] <- fit_90_Autos$mle[length(fit_90_Autos$mle)]
  model_selection_90_Autos[model_name[[i]],"shape_standard_error"] <- fit_90_Autos$se[length(fit_90_Autos$se)]
  fit_90_Autos_list$model_name[[i]] <- fit_90_Autos
  
  fit_180_Autos <- gpd.fit(X_180$Autos[2:length(X_180$Autos)],th_180$Autos[2:length(X_180$Autos)],npy = 252, 
                         ydat = as.matrix(scale(cbind(X_180$Autos[1:length(X_180$Autos)-1],th_180$Autos[1:length(X_180$Autos)-1],
                                                      Z_180_Autos[1:length(X_180$Autos)-1],ZZ_180_Autos[1:length(X_180$Autos)-1]))), 
                         sigl = model_list[[i]],siglink = exp)
  model_selection_180_Autos[model_name[[i]],"nllh"] <- fit_180_Autos$nllh
  model_selection_180_Autos[model_name[[i]],"shape"] <- fit_180_Autos$mle[length(fit_180_Autos$mle)]
  model_selection_180_Autos[model_name[[i]],"shape_standard_error"] <- fit_180_Autos$se[length(fit_180_Autos$se)]
  fit_180_Autos_list$model_name[[i]] <- fit_180_Autos
  
  fit_365_Autos <- gpd.fit(X_365$Autos[2:length(X_365$Autos)],th_365$Autos[2:length(X_365$Autos)],npy = 252, 
                         ydat = as.matrix(scale(cbind(X_365$Autos[1:length(X_365$Autos)-1],th_365$Autos[1:length(X_365$Autos)-1],
                                                      Z_365_Autos[1:length(X_365$Autos)-1],ZZ_365_Autos[1:length(X_365$Autos)-1]))), 
                         sigl = model_list[[i]],siglink = exp)
  model_selection_365_Autos[model_name[[i]],"nllh"] <- fit_365_Autos$nllh
  model_selection_365_Autos[model_name[[i]],"shape"] <- fit_365_Autos$mle[length(fit_365_Autos$mle)]
  model_selection_365_Autos[model_name[[i]],"shape_standard_error"] <- fit_365_Autos$se[length(fit_365_Autos$se)]
  fit_365_Autos_list$model_name[[i]] <- fit_365_Autos
  
  png(file = paste("Figures/POT_Autos/fit_30_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_30_Autos)
  dev.off()
  
  png(file = paste("Figures/POT_Autos/fit_90_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_90_Autos)
  dev.off()
  
  png(file = paste("Figures/POT_Autos/fit_180_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_180_Autos)
  dev.off()
  
  png(file = paste("Figures/POT_Autos/fit_365_0",i,".png"), width = 1200, height = 500)
  par(mfrow=c(1,2))
  gpd.diag(fit_365_Autos)
  dev.off()
  
}

nb_param <- c(0,1,1,1,2,2,2,3)
model_selection_30_Autos["AIC"] <- 2*nb_param + 2*model_selection_30$nllh
model_selection_90_Autos["AIC"] <- 2*nb_param + 2*model_selection_90$nllh
model_selection_180_Autos["AIC"] <- 2*nb_param + 2*model_selection_180$nllh
model_selection_365_Autos["AIC"] <- 2*nb_param + 2*model_selection_365$nllh


# If we use AIC in order to choose which sequences of thresholds use for the rest of the study, we see that the thresholds sequences 
# with a rolling time window of 30 days performe better than the other sequences. Therefore, in the followings of this study we will use 
# this sequence of threshold.

model_selection_30_Autos

# The likelihood ratio test does not allow us to reject the null hypothesis which say that the data comes from the 
# model "2" instead of the model "1-2-3", thus as we also have that the AIC is the smallest for the model "2", we choose the 
# parsimonious model "2".
ratio1 <- 2*(model_selection_30_Autos$nllh[3]-model_selection_30_Autos$nllh[8])
qchisq(0.95,2)
ratio1
p1<-1-pchisq(ratio1,2)
p1


fit_Autos <- gpd.fit(X_30$Autos[2:length(X_30$Autos)],th_30$Autos[2:length(X_30$Autos)],npy = 252,
                   ydat = as.matrix(scale(cbind(th_30$Autos[1:length(X_30$Autos)-1]))),
                   sigl = c(1),siglink = exp)

gpd.diag(fit_Autos)

# We scale the exceedances in order to have a stationnary sequences of exceedances

Y_rescale_Autos <- Y_30$Autos[time_30_Autos]/fit_Autos$vals[,1]
fit_Autos_ <- fpot(Y_rescale_Autos,threshold = 0)
fit_Autos_
# Diagnostics plot of our choosen model
par(mfrow=c(2,2))
plot(fit_Autos_)

# Normal-based confidence interval for the shape parameter
low_conf_shape <- c(fit_Autos_$estimate[[2]]-1.96*fit_Autos_$std.err[[2]],fit_Autos_$estimate[[2]]+1.96*fit_Autos_$std.err[[2]])
low_conf_shape
confint(fit_Autos_)
# The shape is significatively different from zero base on the normal confidence interval

par(mfrow=c(1,2))
plot(profile(fit_Autos_))
# We can confirme that the shape parameter is significatively different from zero based on the profile-likelihood



########## ---------- GEV analysis ---------- ##########
max <- data.frame(row.names = seq(1950,2015,1))
max$Oil <- seq(1950,2015,1)
max$Oil["1952"] <- 0
max

dt <- data
max <- data.frame(row.names = seq(1950,2015,1))
max_ <- c()
for(i in 1950:2015){
  max_ <- c(max_,max(dt$Oil[(as.numeric(str_sub(rownames(dt,1,4))==i))]))
}

annual_max <- function(dt){
  max <- data.frame(row.names = seq(1950,2015,1))
  max_ <- c()
  for(i in 1950:2015){
    max_ <- c(max_,max(dt$Oil[(as.numneric(str_sub(rownames(dt,1,4))==i))]))
    max$Oil <- max_
  }
}
#### ---- Oil ---- ####

annual_max <- 
