######################################################################################
# Modified based on the code from Jak & Cheung (2019) 
######################################################################################
# Without covariates
library(metaSEM)
wd = 'working_directory'
dat = read.csv(paste(wd,'data3.csv',sep=''))

# Data preparation 
Ni 		= dat[,3] # primary study sample sizes
Nstudy 	= nrow(dat) # number of primary studies

vR	= as.matrix(dat[,c(4,6,5)])
MFd	= vector('list',Nstudy)
Mat	= diag(1,3)
for(studyi in 1:Nstudy){
	Mat[lower.tri(Mat)] = vR[studyi,]
	Mat[upper.tri(Mat)] = vR[studyi,]
	MFd[[studyi]] = Mat
}

## Create a dataframe with the data and the asymptotic variances and covariances (acov)
my.df <- Cor2DataFrame(MFd, Ni, acov = "weighted")
## Moderator Female proportion (standardized)
my.df$data <- data.frame(my.df$data,covariate=scale(M),check.names=FALSE)
summary(my.df)

# One-Stage MASEM
## Specify the mediation model
model0 <- "Y ~ M + X; M ~ X; X ~~ 1*X"
RAM0 <- lavaan2RAM(model0, obs.variables = c("X","M","Y"))
RAM0

## Create matrices with implicit diagonal constraints
M0 <- create.vechsR(A0=RAM0$A, S0=RAM0$S, F0=RAM0$F)

## Create heterogeneity variances
T0 <- create.Tau2(RAM=RAM0, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)

## Fit the mediation model with One-Stage MASEM
fit0 <- osmasem(model.name="No moderator", Mmatrix=M0, Tmatrix=T0, data=my.df)
summary(fit0, Saturated=TRUE)

## SRMR
osmasemSRMR(fit0)
## Estimated Tau2
sqrt(diag(VarCorr(fit0)))  


#--------------------------------------------------
# With Covariates
library(metaSEM)
#wd = 'working_directory'
wd = 'D:/Research/180307A/Compare3MASEM/'
source(paste(wd,'RFuncsuncs.R',sep=''))
dat = read.csv(paste(wd,'Med/data3.csv',sep=''))

# Data preparation 
# remove studies with missing values on the moderator
na.id	= which(is.na(dat[,"T1DeprR"])==1) 
dat		= dat[-na.id,]
Ni 		= dat[,3] # primary study sample sizes
Nstudy 	= nrow(dat) # number of primary studies

M		= dat[,"T1DeprR"]	# moderator: baseline depression severity
predM	= c(min(M),median(M),max(M)) # low, medorate and high 
predM	= (predM-mean(M))/sd(M)
M		= (M-mean(M))/sd(M)

vR	= as.matrix(dat[,c(4,6,5)])
MFd	= vector('list',Nstudy)
Mat	= diag(1,3)
for(studyi in 1:Nstudy){
	Mat[lower.tri(Mat)] = vR[studyi,]
	Mat[upper.tri(Mat)] = vR[studyi,]
	MFd[[studyi]] = Mat
}

## Create a dataframe with the data and the asymptotic variances and covariances (acov)
my.df <- Cor2DataFrame(MFd, Ni, acov = "weighted")
## Moderator Female proportion (standardized)
my.df$data <- data.frame(my.df$data,covariate=scale(M),check.names=FALSE)
summary(my.df)

# One-Stage MASEM
## Specify the mediation model
model0 <- "Y ~ M + X; M ~ X; X ~~ 1*X"
RAM0 <- lavaan2RAM(model0, obs.variables = c("X","M","Y"))
RAM0

## Create matrices with implicit diagonal constraints
M0 <- create.vechsR(A0=RAM0$A, S0=RAM0$S, F0=RAM0$F)

## Create heterogeneity variances
T0 <- create.Tau2(RAM=RAM0, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)

## Fit the mediation model with One-Stage MASEM
fit0 <- osmasem(model.name="No moderator", Mmatrix=M0, Tmatrix=T0, data=my.df)
summary(fit0, Saturated=TRUE)

## SRMR
osmasemSRMR(fit0)
## Estimated Tau2
sqrt(diag(VarCorr(fit0)))  


## Mediation model  with `covariate` as a moderator on the A matrix
Ax1 <- RAM0$A
Ax1[grep("\\*", Ax1)] <- "0*data.covariate"
#Ax1[3,2] <- "0"
Ax1

## Create matrices with implicit diagonal constraints
M1 <- create.vechsR(A0=RAM0$A, S0=RAM0$S, F0=RAM0$F, Ax=Ax1)

## Fit the mediation model with moderator
fit1 <- osmasem(model.name="Moderation by covariate", Mmatrix=M1, Tmatrix=T0, data=my.df)
summary(fit1, Saturated=TRUE)

## Get the R2
osmasemR2(fit1, fit0)
sqrt(osmasemR2(fit1, fit0)$Tau2.0)

## Compare the models with and without moderator
anova(fit1, fit0)

# mean indirect effect and its SE
# Delta method
Sigma	= vcov(fit0)[c(1,3),c(1,3)]
est		= coef(fit0)[c(1,3)]
# All moderating effects were not significant
ab.est	= est[1]*est[2]
ab.se	= sqrt(t(est[c(2,1)])%*%Sigma%*%est[c(2,1)])
ab.est
ab.se

### Prediction
# Predicted values
est		= coef(fit1)[c(2,5)]
pred	= cbind(rep(1),predM)%*%est

# Calculate prediction SE using Delta method
Sigma	= vcov(fit1)[c(2,5),c(2,5)]
se.pred = rep(NA,3)
for(i in 1:3){
	se.pred[i] = sqrt(t(c(1,predM[i]))%*%Sigma%*%c(1,predM[i]))
}
cbind(pred,se.pred)




