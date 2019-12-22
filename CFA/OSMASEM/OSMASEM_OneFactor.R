######################################################################################
# Modified based on the code from Jak & Cheung (2019) 
######################################################################################

library(metaSEM)
# Data preparation  
# Exclude primary studies that did not report bivariate correlations
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
model0 <- "SE =~ p1*I1+p3*I3+p4*I4+p7*I7+p10*I10+n2*I2+n5*I5+n6*I6+n8*I8+n9*I9"
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

## Estimated Tau2
round(diag(VarCorr(fit0)),3) 


