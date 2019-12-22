#-------------------------------------------------------
# CFA Random Factor Loadings: Bifactor Model
# No Covariate
#--------------------------------------------------------
library(matrixcalc);library(MASS);library(Matrix)
library(coda);library(R2OpenBUGS);library(metaSEM)

# set up your working directory
wd = 'working_directory'
setwd(wd)
source(paste(wd,'RFuncs.R',sep=''))

# Data preparation
#--------------------------------------------------------
## Exclude studies that did not report bivariate correlations
index <- Gnambs18$CorMat==1
Gnambs18 <- lapply(Gnambs18, function(x) x[index])

# Convert correlation matrices to correlation vectors
mR = Gnambs18$data
vR = sapply(mR,function(x){	x = x[c(1,3,4,7,10,2,5,6,8,9),c(1,3,4,7,10,2,5,6,8,9)] 
	return(x[lower.tri(x)])	})
vR = t(vR)

N		= Gnambs18$n # sample sizes within primary studies
mu.N	= mean(N) # mean sample size
Nstudy	= length(Gnambs18$data) # the number of primary studies
Ninv	= 1/N # reciprocals of sample sizes

# Coordinates of correlation matrices and vectors
p		= 10	# number of variables
pp 		= p*(p-1)/2	# number of bivariate correlations
vi2jk	= Get.vi2jk(p) # from indices of a vectorized matrix to indices of a matrix
j	= vi2jk[,1]; j10 = j+10
k	= vi2jk[,2]; k10 = k+10
vil	= Get.jk2vi(vi2jk,p)
# Do items load on the same factor? 1=No; 0 = Yes
ind = (j>(p+1)/2)*(k<(p+2)/2)

# Covariance matrices of sample correlation vectors
vR.bar		= apply(vR,2,mean,na.rm = TRUE)
vR.impute	= Mimpute(vR,N,'MCAR')
Stau.vR 	= Vj(vR,vR.impute,vR.bar,N,'correlation',pp,Nstudy,j,k,vil)
tau.vR = Stau.vR$tau.vR; S.vR = Stau.vR$S.vR

# hyperparameter for priors (additional error term)
mu.vR.psi	= rep(0,pp)
df.prelim	= 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)

# Matrices for computing ppp
# Compute the between-study covariance matrix of true study-specific correlation vectors
# Z: First derivative of study-specific correlation vectors with respect to model 
#	 parameters (factor loadings)
# NA: for Openbugs to replace with parameter estimates
# The vi_th element in the vectorized correlation matrix corresponds to the correlation
# 		between the j_th and the k_th items. In the bifactor model, the correlation 
#		between the j_th and the_kth items equals the product of the j_th and the_kth 
#		factor loadings plus the product of the (j+10)_th and the (k+10)_th factor 
#		loadings (the factor loadings of the method factors) if the two items are 
#		loaded on the same method factor. Therefore, the first derivative of the vi_th 
#		correlation equals a nonzero value when the derivative is with respect to the 
#		j(+10)_th or the k(+10)_th factor loading and zero when it is with 
#		respect to other SEM parameters 
Z <- matrix(0,pp,p*2)
for(vi in 1:pp){
	Z[vi,c(j[vi],k[vi])] = NA
	Z[vi,c(j[vi]+10,k[vi]+10)] = NA
}
# Diagonal covariance matrix of study-specific model parameters (factor loadings)
# Random factor loadings are assumed to be uncorrelated
V.theta = matrix(0,20,20)
diag(V.theta) = NA

# Model fitting using openbugs
#--------------------------------------------------------
data<-list("Nstudy","N","Ninv","mu.N",'p',"pp","j","k","j10","k10",'ind',
	"vR","tau.vR","mu.vR.psi",'Z','V.theta','alpha.prior.vE','beta.prior.vE') # data

vR.inits = vR.impute;vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA
initsl <- list(list(mu.L=rep(.6,p*2),sd.L=rep(0.1,20),tau.R=100,vR=vR.inits, 
	vR.psi=matrix(0,Nstudy,pp),vR.rep = vR.impute,thetai.raw = matrix(0.6,Nstudy,p*2)))# Initial values
	
prm = c('mu.L','sd.L','ppp') # Parameters to save
model.fn = paste(wd,'CFA/CFARandom.txt',sep='') #  model file name

# stop every 30000 iterations to check whether convergence is achieved
fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=60000,n.burnin=30000,n.thin = 1,
	debug = TRUE,saveExec = TRUE,working.directory = wd)
for(tryi in 1:20){
	fit.coda = read.openbugs(stem="",thin = 10)
	del.id = na.omit(match(c('ppp'),varnames(fit.coda)))
	tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
	print(paste('tryi = ',tryi,sep=''))
	print(tmp.conv)
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){	print(fit,3);break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=30001,n.burnin=1,n.thin=1,restart = TRUE,saveExec = TRUE,working.directory = wd)
	}
}
tryi # number of iterations
summary(fit.coda) # results summary
HPDinterval(fit.coda,prob = .95) # HPD intervals

