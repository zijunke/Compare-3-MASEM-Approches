######################################################################################
# Modified based on the code from Jak & Cheung (2019) 
######################################################################################
library(metaSEM)
# Data preparation  
## Exclude studies that reported CFA results only
index <- Gnambs18$CorMat==1
Gnambs18 <- lapply(Gnambs18, function(x) x[index])

## Create a dataframe with the data and the asymptotic variances and covariances (acov)
my.df <- Cor2DataFrame(Gnambs18$data, Gnambs18$n, acov = "weighted")

## Add the standardized individualism as the moderator
## Standardization of the moderator improves the convergence.
my.df$data <- data.frame(my.df$data,Individualism=scale(Gnambs18$Individualism),check.names=FALSE)
summary(my.df)

# One-Stage MASEM
## Specify the bifactor model
model0 <- "G =~  g1*I1 + g2*I2 + g3*I3 + g4*I4 + g5*I5 + g6*I6 + g7*I7 + g8*I8 + g9*I9 + g10*I10
          POS =~ p1*I1 + p3*I3 + p4*I4 + p7*I7 + p10*I10
          NEG =~ n2*I2 + n5*I5 + n6*I6 + n8*I8 + n9*I9"
RAM0 <- lavaan2RAM(model0, obs.variables = paste0("I", 1:10),std.lv = TRUE)
RAM0

## Create matrices with implicit diagonal constraints
M0 <- create.vechsR(A0=RAM0$A, S0=RAM0$S, F0=RAM0$F)

## Create heterogeneity variances
T0 <- create.Tau2(RAM=RAM0, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)

## Fit the bifactor model with One-Stage MASEM
fit0 <- osmasem(model.name="No moderator", Mmatrix=M0, Tmatrix=T0, data=my.df)
summary(fit0, Saturated=TRUE)

## SRMR
osmasemSRMR(fit0)
## Estimated Tau2
SD.r = matrix(NA,10,10)
SD.r[lower.tri(SD.r,diag=FALSE)] = sqrt(diag(VarCorr(fit0)))
round(SD.r,2)

#------------------------------------------------------------
#	Studies with missing values on individualism were removed
#------------------------------------------------------------
library(metaSEM)
# Data preparation  
## Exclude studies with missing values on Individualism
index_na <- is.na(Gnambs18$Individualism)
Gnambs18 <- lapply(Gnambs18, function(x) x[!index_na])

## Exclude studies reported CFA results only
index <- Gnambs18$CorMat==1
Gnambs18 <- lapply(Gnambs18, function(x) x[index])

## Create a dataframe with the data and the asymptotic variances and covariances (acov)
my.df <- Cor2DataFrame(Gnambs18$data, Gnambs18$n, acov = "weighted")

## Add the standardized individualism as the moderator
## Standardization of the moderator improves the convergence.
my.df$data <- data.frame(my.df$data,Individualism=scale(Gnambs18$Individualism),       check.names=FALSE)
summary(my.df)

# One-Stage MASEM
## Specify the bifactor model
model0 <- "G =~  g1*I1 + g2*I2 + g3*I3 + g4*I4 + g5*I5 + g6*I6 + g7*I7 + g8*I8 + g9*I9 + g10*I10
          POS =~ p1*I1 + p3*I3 + p4*I4 + p7*I7 + p10*I10
          NEG =~ n2*I2 + n5*I5 + n6*I6 + n8*I8 + n9*I9"
RAM0 <- lavaan2RAM(model0, obs.variables = paste0("I", 1:10))
RAM0

## Create matrices with implicit diagonal constraints
M0 <- create.vechsR(A0=RAM0$A, S0=RAM0$S, F0=RAM0$F)

## Create heterogeneity variances
T0 <- create.Tau2(RAM=RAM0, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)

## Fit the bifactor model with One-Stage MASEM
fit0 <- osmasem(model.name="No moderator", Mmatrix=M0, Tmatrix=T0, data=my.df)
summary(fit0, Saturated=TRUE)

## SRMR
osmasemSRMR(fit0)

## Bifactor model  with `Individualism` as a moderator on the A matrix
## Create the A1 matrix with moderator effects of "Individualism"
Ax1 <- RAM0$A
Ax1[grep("\\*", Ax1)] <- "0*data.Individualism"
#Ax1[10,12] <- "0"
Ax1

## Create matrices with implicit diagonal constraints
M1 <- create.vechsR(A0=RAM0$A, S0=RAM0$S, F0=RAM0$F, Ax=Ax1)

## Fit the bifactor model with moderator
fit1 <- osmasem(model.name="Moderation by individualism", Mmatrix=M1, Tmatrix=T0, data=my.df)
summary(fit1, Saturated=TRUE)

## Get R2 in the metric of bivariate correlations
osmasemR2(fit1, fit0)
vR2 = osmasemR2(fit1, fit0)$R2
MR2 = matrix(NA,10,10)
MR2[lower.tri(MR2,diag=FALSE)] = vR2
round(MR2,2)

## Compare the models with and without moderator
anova(fit1, fit0)

# Compute predicted values and their SEs
# Delta method was used
### SE for prediction error
Sigma = vcov(fit1)[c(1:10,21:30),c(1:10,21:30)]
est = coef(fit1)[c(1:10,21:30)]
Moderator = my.df$data$Individualism[c(6,26)] #(182,-126)

Z = c(1,NA)
SE.pred = matrix(NA,10,2)
pred = matrix(NA,10,2)
for(i in 1:10){
	pred[i,1] <- est[i]+est[i+10]*Moderator[1]
	Z[2] <- Moderator[1]
	SE.pred[i,1] = sqrt(t(Z)%*%Sigma[c(i,i+10),c(i,i+10)]%*%Z)
	
	pred[i,2] <- est[i]+est[i+10]*Moderator[2]
	Z[2] <- Moderator[2]
	SE.pred[i,2] = sqrt(t(Z)%*%Sigma[c(i,i+10),c(i,i+10)]%*%Z)
}
ord = c(1,3,4,7,10,2,5,6,8,9)
# predicted values, SEs, z-stats
round(cbind(pred[ord,1],SE.pred[ord,1],pred[ord,1]/SE.pred[ord,1]),2)
round(cbind(pred[ord,2],SE.pred[ord,2],pred[ord,2]/SE.pred[ord,2]),2)


