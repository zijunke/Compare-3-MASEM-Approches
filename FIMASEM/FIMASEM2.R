######################################################################################
# Modified based on the supplemental materials of Yu, Downes, Carter, and O'Boyle (2018) 
######################################################################################
library(metaSEM);library(matrixcalc);library(OpenMx);library(Matrix);library(MASS)

# Read in data
#wd = 'D:/Research/180307A/Compare3MASEM/'
wd = 'your_working_directory'
dat = read.csv(paste(wd,'data3.csv',sep=""))

# Data preparation
Ni		= dat[,3] 	# sample sizes for primary studies
Nstudy	= nrow(dat) # number of primary studies
reps	= 10000	# number of replications
NpS		= 1/mean(1/Ni) # harmonic mean of sample size per study

#Reformat the data (bivariate correlations) for TSSEM input
cormats = list(data = vector('list',Nstudy))
vnames	= c('X','M','Y')
vR 	= as.matrix(dat[,c(4,6,5)])
M	= diag(1,3)
for(si in 1:Nstudy){
	M[lower.tri(M)] = vR[si,]
	M[upper.tri(M)] = vR[si,]
	cormats[[si]] = M
	colnames(cormats[[si]]) <- vnames
	rownames(cormats[[si]]) <- vnames
}

# Specify mediation model
#----------------------------------------------------
A <- create.mxMatrix(c(0,0,0,".3*a",0,0,".3*c",".3*b",0),ncol=3,nrow=3,byrow=T,name="A")
dimnames(A) <- list(vnames,vnames)
			  
S <- create.mxMatrix(c(1,0,".2*varM",0,0,".2*varY"),type="Symm",ncol=3,nrow=3,byrow=T,name="S")
dimnames(S) <- list(vnames,vnames)
matrF <- mxMatrix(type="Iden",nrow=3,ncol=3,name="F")
exp <- mxExpectationRAM("A","S","F", dimnames=vnames )

# Conduct multivariate FIMASEM
#----------------------------------------------------
#run TSSEM Step 1
step.one <- tssem1(cormats,Ni,method="REM",RE.type="Diag") 
rho.mult <- diag(1,nrow=length(vnames),ncol=length(vnames))
rho.mult[lower.tri(rho.mult)] <- (coef(step.one,select="fixed"))
temp <- t(rho.mult)
rho.mult[upper.tri(rho.mult)] <- temp[upper.tri(temp)]
sigma.mult <- diag(coef(step.one,select="random")) 
dimnames(rho.mult) <- list(vnames,vnames)
matrices.mult <- rCorPop(rho.mult,sigma.mult,corr=T,k=reps,nonPD.pop="nearPD")

# Run SEM on those random matrices using the same technique we use for V&O FIMASEM
coefs.fits.multivariate.FIMASEM <- as.data.frame(t(sapply(1:reps,function(i) {
  openmxmodel <- mxModel("temp",mxData(matrices.mult[[i]],type="cov",numObs = NpS),
			matrA = A,matrS = S,matrF = matrF,exp=exp,funML=mxFitFunctionML());
  openmxfit <- mxRun(openmxmodel,silent=T);
  if (openmxfit$output$status[[1]] == 6) {	openmxfit <- mxRun(openmxfit,silent=T)	}
  modelsummary <- summary(openmxfit)
  coefs <- coef(openmxfit)
  output <- c(coefs,openmxfit$output$status[[1]])
  names(output) <- c(names(coefs),'openMxStatus')
  output
}))) #returns a dataframe of SEM parameter estimates (i.e., fit indices and path coefficients)

# Get Results
#----------------------------------------------------
del.id = which(coefs.fits.multivariate.FIMASEM[,6]>0)
if(length(del.id)>0){
	m <- sapply(coefs.fits.multivariate.FIMASEM,function(x){mean(x[-del.id])})
	sdv <- sapply(coefs.fits.multivariate.FIMASEM,function(x){sd(x[-del.id])}) 
	UL <- m + qnorm(.90)*sdv # upper limits
	LL <- m - qnorm(.90)*sdv # lower limits
}else{
	m <- sapply(coefs.fits.multivariate.FIMASEM,function(x){mean(x)})
	sdv <- sapply(coefs.fits.multivariate.FIMASEM,function(x){sd(x)}) 
	UL <- m + qnorm(.90)*sdv # upper limits
	LL <- m - qnorm(.90)*sdv # lower limits
}
# SEM parameter means, sds, and lower & upper limits of credibility intervals
round(cbind(m,sdv,LL,UL),3)

# Mean indirect effect
# mean(a)*mean(b) + cov(a,b)
m[1]*m[3]+cov(coefs.fits.multivariate.FIMASEM)[1,3]





