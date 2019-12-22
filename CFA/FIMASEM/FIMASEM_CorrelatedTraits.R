######################################################################################
# Modified based on the supplemental materials of Yu, Downes, Carter, and O'Boyle (2018) 
######################################################################################
library(metaSEM);require('matrixcalc');library(OpenMx);library(Matrix);library(MASS)

mySRMR <- function(oC,mC){
  p = nrow(oC)
  return(sqrt(sum((oC-mC)^2)/p/(p+1)))
}

# Data preparation
#----------------------------------------------------------
## Remove studies that did not report bivariate correlations
index <- Gnambs18$CorMat==1
Gnambs18 <- lapply(Gnambs18, function(x) x[index])
Ni = Gnambs18$n # sample sizes within primary studies
NpS = 1/mean(1/Ni) # harmonic mean of sample sizes within primary studies
#Nps = mean(Ni)
k = length(Gnambs18$data) # number of primary studies
reps <- 10000 # number of replications

#Reformat the data for TSSEM input
vnames <- paste('I',1:10)
cormats <- lapply(Gnambs18$data,function(x) x = x[c(1,3,4,7,10,2,5,6,8,9),c(1,3,4,7,10,2,5,6,8,9)])

# Conduct multivariate FIMASEM
#----------------------------------------------------------
# Stage 1: Run TSSEM to obtain mean correlation vector and its covariance matrix
step.one <- tssem1(cormats,Ni,method="REM",RE.type="Diag")
rho.mult <- diag(1,nrow=length(vnames),ncol=length(vnames))
rho.mult[lower.tri(rho.mult)] <- (coef(step.one,select="fixed"))
temp <- t(rho.mult)
rho.mult[upper.tri(rho.mult)] <- temp[upper.tri(temp)]
sigma.mult <- diag(coef(step.one,select="random")) 
dimnames(rho.mult) <- list(vnames,vnames)

# Stage 2
# Step 1: Simulated a large numbers of correlation vectors
matrices.mult <- rCorPop(rho.mult,sigma.mult,corr=T,k=reps+500,nonPD.pop="nearPD")
matrices.mult <- matrices.mult[which(sapply(matrices.mult,is.positive.definite))][1:reps]

# Step 2: Fit the studied SEM model to each of the simulated correlation vectors
# CFA formulation
# Factor loading matrix
L.values = matrix(c(rep(0.6,5),rep(0,10),rep(0.6,5)),10,2,byrow = F)
L.lbound = matrix(c(0,rep(-1,4),rep(0,11),rep(-1,4)),10,2,byrow = F)
L.ubound = matrix(c(rep(1,5),rep(0,10),rep(1,5)),10,2,byrow = F)
L.labels = matrix(c(paste('L',1:5,sep=''),rep(NA,10),paste('L',6:10,sep='')),10,2,byrow = F)
L.free = L.values!=0
L <- mxMatrix(type = 'Full',free = L.free,values = L.values,labels = L.labels,
              lbound = L.lbound,ubound = L.ubound,name="L") 

# Factor correlation matrix			  
Phi <- mxMatrix(type = 'Symm',nrow = 2,ncol = 2,free = c(FALSE,TRUE,FALSE),values = c(1,.3,1),
	labels = c(NA,'rho',NA),name = 'Phi',lbound = rep(-1,3),ubound = rep(1,3))

# Uniqueness
U.values = diag(0.1,10)
U.lbound = diag(0,10)
U.ubound = diag(1,10)
U.labels = matrix(NA,10,10)
diag(U.labels) = paste('s2e',1:10,sep='')
U.free = U.values!=0
U <- mxMatrix(type = 'Diag',free = U.free,values = U.values,labels = U.labels,
     lbound = U.lbound,ubound = U.ubound, name="U") 

# CFA model-implied covariance matrix
ecov <- mxAlgebra(L%*%Phi%*%t(L) + U , name="expCov")
expectation <- mxExpectationNormal(cov="expCov",dimnames = vnames)

# Run SEM on those random matrices
coefs.fits.multivariate.FIMASEM <- as.data.frame(t(sapply(1:reps,function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.mult[[i]],type="cov",numObs = NpS),
		L,Phi,U,ecov, expectation,funML=mxFitFunctionML());
  openmxfit <- mxRun(openmxmodel,silent=T);
  if (openmxfit$output$status[[1]] == 6) {	openmxfit <- mxRun(openmxfit,silent=T)	}
  modelsummary <- summary(openmxfit);
  coefs <- coef(openmxfit)
  coef.names <- names(coefs)
  mC <- openmxfit$expCov$result
  oC <- matrices.mult[[i]]
  output <- c(coefs,mySRMR(oC,mC),modelsummary$CFI,
	openmxfit$output$status[[1]]); 
  names(output) <- c(names(coefs),'SRMR','CFI','openMxStatus')
  output
}))) #returns a dataframe of SEM parameter estimates (i.e., fit indices and path coefficients)


# Get results
#----------------------------------------------------------
del.id = which(coefs.fits.multivariate.FIMASEM[,24]>0)
m <- sapply(coefs.fits.multivariate.FIMASEM,function(x){mean(x[-del.id])})
sdv <- sapply(coefs.fits.multivariate.FIMASEM,function(x){sd(x[-del.id])}) 
UL <- m + qnorm(.90)*sdv # upper limits
LL <- m - qnorm(.90)*sdv # lower limits
# SEM parameter means, sds, and lower & upper limits of credibility intervals
round(cbind(m,sdv,LL,UL),3)

print('% SRMR < .10')
print(sum(coefs.fits.multivariate.FIMASEM$SRMR[-del.id] < .1)/(reps-length(del.id)))

print('% CFI > .90')
print(sum(coefs.fits.multivariate.FIMASEM$CFI[-del.id] > .90)/(reps-length(del.id)))# 

length(del.id)/reps



