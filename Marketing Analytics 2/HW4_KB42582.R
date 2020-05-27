## Loaing the data
rm(list = ls())

setwd("C:/Users/khyat/Downloads/MA2")
sow =  read.csv('CreditCard_SOW_data2.csv')
attach(sow)

## Question 1
library(censReg)
cr1=censReg(SOW~Promotion+Balance, left=0, right = 1, data=sow)
summary(cr1)

## Question2 


library(truncnorm)
library(mnormt)

rm(list = ls())

#Bayesian estimation for truncated regression
#stage 1. 
##Reading data into R and creating columns for censored data
DataFile = "CreditCard_SOW_data2.csv"
sow.data = read.csv(DataFile, header=T)
sow.data$Cens0 = (sow.data$SOW==0)*1
sow.data$Cens1 = (sow.data$SOW==1)*1

#Extracting right and left censored data
sow.XRC = cbind(1, as.matrix(sow.data[sow.data$Cens0==1, 3:4]))
sow.XLC = cbind(1, as.matrix(sow.data[sow.data$Cens1==1, 3:4])) 
sow.X = cbind(1, as.matrix(sow.data[, 3:4]))
sow.X2 = t(sow.X)%*%sow.X
nRC = dim(sow.XRC)[1]
nLC = dim(sow.XLC)[1]
nObs = dim(sow.data)[1]

#Stage 2. 
##Initial Setup for the algorithm
NIT = 10000       #num of interations
nBurn = 2000      #num of burn-ins  
NIT.eff = NIT - nBurn    #effective sample size
thin.step = 10           #thinning  
NIT.thin = floor(NIT.eff/thin.step)   #effective sample size after thinning

#Stage 3. 
##Record Posterior samples
beta.dim = 3
beta.pos = matrix(0, NIT.thin, beta.dim)
tau.pos = rep(0, NIT.thin)

#stage 4. priors
#for Beta: mNormal(mu.beta, sigma.beta)
mu.beta = rep(0,beta.dim) 
sigma.beta = 1E6 * diag(beta.dim)  
iSigma.beta = 1E-6 * diag(beta.dim)  #inverse prior covariance matrix 

#prior for precision: Gamma(a.tau, b.tau)
a.tau = 1/2
b.tau = 1/2

#stage 5. 
##Gibbs sampler

#initialize the loop
curBeta = c(0.5, 0, 0) 
curTau = 4
g = 1

#main loop
for (m in 1:NIT){
  #step 1. sample latent SOW 
  #step 1.a. sample SOW right-censored at 0
  #Please fill in the blank below the code for sampling the latent SOW when the observed SOW=0
  #Please name your sampled latent SOW curYRC
  curYRC = rtruncnorm(nRC,b=0, mean = sow.XRC%*%curBeta, sd = sqrt(1/curTau))
  
  
  #step 1.b sample SOW left-censored at 1
  #Please fill in the blank below the code for sampling the latent SOW when the observed SOW=1
  #Please name your sampled latent SOW curYLC
  curYLC = rtruncnorm(nLC,a=1, mean = sow.XLC%*%curBeta, sd = sqrt(1/curTau))
  
  #step 2 sample beta
  #step 2.a impute the latent variables
  sow.Y = sow.data$SOW
  sow.Y[sow.data$Cens0==1] = curYRC
  sow.Y[sow.data$Cens1==1] = curYLC
  #step 2.b sample beta
  sigma.hat = solve(curTau*sow.X2 + iSigma.beta)
  betaPos.mn = sigma.hat%*%(curTau*t(sow.X)%*%sow.Y + iSigma.beta%*%mu.beta)
  curBeta = as.vector(rmnorm(1, mean=betaPos.mn, varcov=sigma.hat)) 
  
  #step 3 sample tau (precision = 1/sigma^2)
  sowE.hat = sow.Y-sow.X%*%curBeta
  curTau = rgamma(1, 0.5*nObs+a.tau, 0.5*t(sowE.hat)%*%sowE.hat+b.tau)
  
  #save thinned samples after burn-ins
  if ((m > nBurn) & (m%%thin.step == 0)) {
    beta.pos[g,] = curBeta
    tau.pos[g] = curTau
    g = g+1
  }
}

plot(beta.pos[,1], type = "l",ylab = paste0('Beta - 0'), main = paste0('Beta - 0: Posterior Sampling chain'))

plot(beta.pos[,2], type = "l",ylab = paste0('Beta - 1'), main = paste0('Beta - 1: Posterior Sampling chain'))

plot(beta.pos[,3], type = "l",ylab = paste0('Beta - 2'), main = paste0('Beta - 2: Posterior Sampling chain'))

plot(tau.pos, type = "l",ylab = paste0('Tau'), main = paste0('Tau: Posterior Sampling chain'))

hist(beta.pos[,1],xlab = paste0('Beta - ',0), main = paste0('Beta - 0: Histogram'))
hist(beta.pos[,2],xlab = paste0('Beta - ',1), main = paste0('Beta - 1: Histogram'))
hist(beta.pos[,3],xlab = paste0('Beta - ',2), main = paste0('Beta - 2: Histogram'))


hist(tau.pos,xlab = paste0('Tau '), main = paste0( 'Tau: Histogram'))

quantile(beta.pos[,1],probs=c(0.025, 0.5, 0.975))
quantile(beta.pos[,2],probs=c(0.025, 0.5, 0.975))
quantile(beta.pos[,3],probs=c(0.025, 0.5, 0.975))
quantile(tau.pos,probs=c(0.025, 0.5, 0.975))

mean(beta.pos[,3])
mean(tau.pos)