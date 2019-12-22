#-------------------------------------------------------------------------------
# Mediation: The mediating mechanism underlying mindfulness-based interventions
# Random a,b,cp paths
#--------------------------------------------------------------------------------
library(matrixcalc);library(MASS);library(coda);library(Matrix);library(R2OpenBUGS)

# set working directory
wd = 'workding_directory'
setwd(wd)
source('RealData.R')

# Data preparation
#--------------------------------------------------------
dat = read.csv('data2.csv') # Read in data

vR	= as.matrix(dat[,c(4,6,5)]) # individual study bivariate correlations
N	= dat[,3] # individual study sample sizes
Nstudy	= nrow(dat) # number of studies
mu.N	= mean(N) # mean sample size per study

# Coordinations (matrix <-> vector)
p	= 3 # number of observed variables
pp 	= p*(p-1)/2 # number of bivariate correlations
vi2jk	= Get.vi2jk(p)
j	= vi2jk[,1] # row index from vector to matrix
k	= vi2jk[,2] # column index from vector to matrix
vil	= Get.jk2vi(vi2jk,p) # from matrix to vector

# Compute level-1 error covariance matrix 
# Or covariance matrix of observed correlation vectors 
vR.bar = apply(vR,2,mean,na.rm = TRUE); vR.impute = Mimpute(vR,N,'MCAR')
Stau.vR <- Vj(vR,vR.impute,vR.bar,N,'correlation',pp,Nstudy,j,k,vil)
tau.vR <- Stau.vR$tau.vR; S.vR <- Stau.vR$S.vR

# Hyperparameters for priors (additional error term)
I3 = diag(1,3); u0 = rep(0,3);mu.vR.psi = rep(0,pp)
df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2; beta.prior.vE = alpha.prior.vE*(0.3/mu.N)

# Model fitting using openbugs
#--------------------------------------------------------
data<-list("Nstudy","N","mu.N",'p',"pp","vR","tau.vR","mu.vR.psi",
	'alpha.prior.vE','beta.prior.vE','u0','I3') #data

vR.inits = vR.impute; vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA
initsl <- list(list(mu.a=0,mu.b=0,mu.cp=0,tau.u=diag(100,3),xi=rep(1,3),tau.R = 100,
	vR.psi = matrix(0,Nstudy,pp),vR = vR.inits,vR.rep = vR.impute)) # initial values

prm = c('mu.a','mu.b','mu.cp','mu.ab','sd.a','sd.b','sd.cp',
	'rho.ab','rho.acp','rho.bcp') # Parameters to save; 
model.fn = 'Mediation_Random.txt' # Model file name

# stop every 30000 iterations to check whether convergence is achieved
fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=60000,n.burnin=30000,n.thin=1,
	debug = FALSE,saveExec = TRUE,working.directory = work.d)

for(tryi in 1:20){
	fit.coda <- as.mcmc.list(fit)
	tmp.conv = geweke.diag(fit.coda)[[1]]$z
	print(paste('tryi = ',tryi,sep=''))
	print(tmp.conv)
	if(sum((abs(tmp.conv)>1.96),na.rm = TRUE)==0){	print(fit,3);break
	}else{	
		fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=30001,n.burnin=1,n.thin=1,
		restart = TRUE,saveExec = TRUE,working.directory = work.d)
	}
}

tryi
fit.coda = read.openbugs(stem="",thin = 1)
summary(fit.coda) # results summary
HPDinterval(fit.coda,prob = .95) # HPD intervals


