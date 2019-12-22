#-------------------------------------------------------------------------------
# Mediation: the mediating mechanism underlying mindfulness-based interventions
# Covariate: Baseline depression severity 
#-------------------------------------------------------------------------------
library(matrixcalc);library(MASS);library(coda);library(Matrix);library(R2OpenBUGS)

# set working directory
wd = 'working_directory'
setwd(wd)
source('RealData.R')

# Data preparation
#--------------------------------------------------------
# Read in data
# remove studies with missing values on the moderator
dat = read.csv('data3.csv')
na.id	= which(is.na(dat[,"T1DeprR"])==1) 
dat		= dat[-na.id,]

N		= dat[,3]   # sample sizes within primary studies
mu.N	= mean(N)   # mean sample size
Nstudy	= nrow(dat) # the number of primary studies

vR		= as.matrix(dat[,c(4,6,5)]) # correlation vectors from primary studies
M		= dat[,"T1DeprR"] # moderator: baseline depression severity
predM	= c(min(M),median(M),max(M)) # Low,moderator, and high levels of baseline depression
predM	= (predM-mean(M))/sd(M,na.rm = T) # standardization
M		= (M-mean(M))/sd(M)  # standardization

# Coordinations (matrix <-> vector)
p 	= 3 # number of observed variables
pp 	= p*(p-1)/2 # number of bivariate correlations
vi2jk	= Get.vi2jk(p)
j	= vi2jk[,1] # row index from vector to matrix
k	= vi2jk[,2] # column index from vector to matrix
vil	= Get.jk2vi(vi2jk,p) # from matrix to vector

# Compute level-1 error covariance matrix 
# Or covariance matrix of observed correlation vectors 
vR.bar = apply(vR,2,mean,na.rm = TRUE)
vR.impute = Mimpute(vR,N,'MCAR')
Stau.vR <- Vj(vR,vR.impute,vR.bar,N,'correlation',pp,Nstudy,j,k,vil)
tau.vR <- Stau.vR$tau.vR; S.vR <- Stau.vR$S.vR

# Hyperparameter for priors (additional error term)
mu.vR.psi	= rep(0,pp)
df.prelim = 100*pp/mu.N+pp
alpha.prior.vE = (df.prelim-pp+1)/2; beta.prior.vE = alpha.prior.vE*(0.3/mu.N)

# Model fitting using openbugs
#--------------------------------------------------------
data<-list("Nstudy","N","mu.N","pp","vR","tau.vR",'M','predM',
	"mu.vR.psi",'alpha.prior.vE','beta.prior.vE') # data

vR.inits = vR.impute; vR.inits[which(is.na(vR)==0,arr.ind = TRUE)] = NA
initsl <- list(list(b0.a=0,b0.b=0,b0.cp=0,b1.a=0,b1.cp=0,sd.ua=0.1,sd.ucp=0.1,
	tau.R=100,a = rep(0,Nstudy),cp = rep(0,Nstudy),vR.psi = matrix(0,Nstudy,pp),
	vR=vR.inits,vR.rep = vR.impute))# initial values
	
prm = c(paste('b0.',c('a','b','cp'),sep=''),paste('b1.',c('a','cp'),sep=''),
	paste('sd.u',c('a','cp'),sep=''),'cphat','cphat.e')# Parameters to save
model.fn = paste(wd,'Mediation_Covariate2.txt',sep='')# Model file name

# stop every 30000 iterations to check whether convergence is achieved
fit = bugs(data,initsl,prm,model.fn,n.chains=1,n.iter=60000,n.burnin=30000,n.thin=1,
	debug = TRUE,saveExec = TRUE,working.directory = work.d)
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





