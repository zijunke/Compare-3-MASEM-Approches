#-------------------------------------------------------
# CFA Random Factor Loadings: Bifactor Model
# Covariate: Individualism
#--------------------------------------------------------
library(matrixcalc);library(MASS);library(coda);library(Matrix);library(R2OpenBUGS);library(metaSEM)

# set up your working directory
wd = 'working_directory'
setwd(wd)
source(paste(wd,'RFuncs.R',sep=''))

# Data preparation
#-------------------------------------------------
## Remove studies that did not report Individualism or bivariate correlations
index_na <- is.na(Gnambs18$Individualism)
Gnambs18 <- lapply(Gnambs18,function(x) x[!index_na])
index 	 <- Gnambs18$CorMat==1
Gnambs18 <- lapply(Gnambs18, function(x) x[index])

# Standardize Individualism
M <- Gnambs18$Individualism
M <- (M-mean(M))/sd(M)

# Convert correlation matrices to correlation vectors
mR = Gnambs18$data
vR = sapply(mR,function(x){	x = x[c(1,3,4,7,10,2,5,6,8,9),c(1,3,4,7,10,2,5,6,8,9)]
	return(x[lower.tri(x)])	})
vR = t(vR)

N		= Gnambs18$n # sample sizes within primary studies
mu.N	= mean(N) # mean sample size
Nstudy	= length(Gnambs18$data) # the number of primary studies

# Coordinates of correlation matrices and vectors
p		= 10	# number of variables
pp 		= p*(p-1)/2	# number of bivariate correlations
vi2jk	= Get.vi2jk(p) # from indices of a vectorized matrix to indices of a matrix
j	= vi2jk[,1]; j10 = j+10
k	= vi2jk[,2]; k10 = k+10
vil	= Get.jk2vi(vi2jk,p)
# Do items load on the same factor? 1=No; 0 = Yes
ind 	= (j>(p+1)/2)*(k<(p+2)/2) 

# Covariance matrices of sample correlation vectors
vR.bar = apply(vR,2,mean,na.rm = TRUE)
vR.impute = Mimpute(vR,N,'MCAR')
Stau.vR <- Vj(vR,vR.impute,vR.bar,N,'correlation',pp,Nstudy,j,k,vil)
tau.vR <- Stau.vR$tau.vR ;S.vR <- Stau.vR$S.vR

mu.vR.psi	= rep(0,pp)	
df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2
beta.prior.vE = alpha.prior.vE*(0.3/mu.N)

# Model fitting using openbugs
data<-list("Nstudy","N","mu.N",'p',"pp","j","k",'j10','k10','ind',
   "vR","tau.vR",'M',"mu.vR.psi",'alpha.prior.vE','beta.prior.vE') # data

vR.inits = vR.impute;vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA
vL.inits = matrix(0.6,Nstudy,p*2);vL.inits[,15] = NA
initsl <- list(list(a=rep(0,p*2),b=c(rep(0,14),NA,rep(0,5)),sd.uL = c(rep(0.1,14),NA,rep(0.1,5)),tau.R = 100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,vR.rep = vR.impute,vL = vL.inits))# initial values

prm =c('a','b','sd.uL','vLH.pred','vLL.pred','vLH.prede','vLL.prede') # Parameters to save
model.fn = paste(wd,'CFACovariate.txt',sep='') # model file name

# stop every 30000 iterations to check whether convergence is achieved
fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=60000,n.burnin=30000,n.thin = 1,
	debug = FALSE,saveExec = TRUE,working.directory = wd)
for(tryi in 1:30){
	fit.coda = read.openbugs(stem="",thin = 10)
	tmp.conv = geweke.diag(fit.coda[,1:61])[[1]]$z
	print(paste('tryi = ',tryi,sep=''))
	print(tmp.conv)

	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){	print(fit,3);break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=30001,n.burnin=1,n.thin = 1,
		restart = TRUE,saveExec = TRUE,working.directory = wd)
	}
}
tryi # number of iterations
summary(fit.coda) # results summary
HPDinterval(fit.coda,prob = .95) # HPD intervals