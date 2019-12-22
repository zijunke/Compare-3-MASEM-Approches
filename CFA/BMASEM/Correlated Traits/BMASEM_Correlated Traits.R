##-------------------------------------------------------
# CFA Random Factor Loadings Correlated Traits Model
#--------------------------------------------------------
library(matrixcalc);library(MASS);library(Matrix)
library(coda);library(R2OpenBUGS);library(metaSEM)

# Set up working directory
wd = 'Working_directory' 
setwd(wd)
source(paste(wd,'RealData.R',sep=''))

# Data preparation
#-------------------------------------------------
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
j = vi2jk[,1]; k = vi2jk[,2]; vil = Get.jk2vi(vi2jk,p)



# Covariance matrices of sample correlation vectors
vR.bar		= apply(vR,2,mean,na.rm = TRUE)
vR.impute	= Mimpute(vR,N,'MCAR')
Stau.vR 	= Vj(vR,vR.impute,vR.bar,N,'correlation',pp,Nstudy,j,k,vil)
tau.vR = Stau.vR$tau.vR; S.vR = Stau.vR$S.vR
# Do items load on the same factor? 1=No; 0 = Yes
ind = (j>(p+1)/2)*(k<(p+2)/2)

# Matrices for computing ppp
# Compute the between-study covariance matrix of true study-specific correlation vectors
# Z: First derivative of study-specific correlation vectors with respect to model 
#		parameters (factor loadings)
# NA: for Openbugs to replace with parameter estimates
# The vi_th element in the vectorized correlation matrix corresponds to the correlation
# 		between the j_th and the k_th items. In the correlated traits model, the 
#		correlation between the j_th and the_kth items equals the product of the j_th 
#		and the_kth factor loadings (and the factor correlation if the two items are 
#		not loaded on the same factor). Therefore, the first derivative of the vi_th 
#		correlation equals a nonzero value when the derivative is with respect to the 
#		j_th or the k_th factor loading or the factor correlation, and zero when it is 
#		with respect to other SEM parameters 
Z = matrix(0,pp,p+1)
for(vi in 1:pp){	Z[vi,c(j[vi],k[vi])] = NA	}	
Z[,p+1] = NA
# Diagonal covariance matrix of study-specific model parameters (factor loadings)
# Random factor loadings are assumed to be uncorrelated
V.theta = matrix(0,11,11)
diag(V.theta) = NA

# hyperparameter for priors (additional error term)
mu.vR.psi	= rep(0,pp)
df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)

# Model fitting using openbugs
#--------------------------------------------------------
data<-list("Nstudy","N","Ninv","mu.N",'p',"pp","j","k",'ind','V.theta',
	"vR","tau.vR","mu.vR.psi",'Z','alpha.prior.vE','beta.prior.vE') # data

vR.inits = vR.impute; vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA
initsl <- list(list(mu.L=rep(.6,p),mu.rho = 0,sd.L = rep(0.1,p),sd.rho = 0.1,tau.R=100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,vR.rep = vR.impute)) # Initial values

prm = c('mu.L','sd.L','mu.rho','sd.rho','ppp') # parameters to save;
model.fn = paste(wd,'CFARandom.txt',sep='') # model file name

# stop every 30000 iterations to check whether convergence is achieved
fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=60000,n.burnin=30000,n.thin = 1,
	debug = TRUE,saveExec = TRUE,working.directory = wd)
for(tryi in 1:20){
	fit.coda = read.openbugs(stem="",thin = 10)
	del.id = na.omit(match(c('ppp'),varnames(fit.coda)))
	tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
	print(paste('tryi = ',tryi,sep=''))
	print(tmp.conv)
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){	break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=30001,n.burnin=1,n.thin=1,restart = TRUE,saveExec = TRUE,working.directory = wd)
	}
}
tryi # number of iterations
summary(fit.coda) # results summary
HPDinterval(fit.coda,prob = .95) # HPD intervals



