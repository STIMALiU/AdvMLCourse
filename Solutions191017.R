# Solutions to the 732A96/TDDE15 exam on 19/10-2017.

# Question 1: PGMs

# Learn a BN from the Asia dataset, find a separation statement (e.g. B _|_ E | S, T) and, then, check that
# it corresponds to a statistical independence.

library(bnlearn)
library(gRain)
set.seed(123)
data("asia")
hc3<-hc(asia,restart=10,score="bde",iss=10)
plot(hc3)
hc4<-bn.fit(hc3,asia,method="bayes")
hc5<-as.grain(hc4)
hc6<-compile(hc5)
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("yes","yes","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("yes","yes","no"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("yes","no","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("yes","no","no"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("no","yes","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("no","yes","no"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("no","no","yes"))
querygrain(hc7,c("B"))
hc7<-setFinding(hc6,nodes=c("S","T","E"),states=c("no","no","no"))
querygrain(hc7,c("B"))

# Sample DAGs at random and, then, check which of them coincide with their CPDAGs, i.e. the CPDAGs have no
# undirected edge (recall that a CPDAG has an undirected edge if and only if there are two DAGs in the Markov
# equivalence class that differ in the direction of that edge).

#  The exact ratio according to the literature is 11.2

library(bnlearn)
set.seed(123)
ss<-50000
x<-random.graph(c("A","B","C","D","E"),num=ss,method="melancon",every=50,burn.in=30000)

y<-unique(x)
z<-lapply(y,cpdag)

r=0
for(i in 1:length(y)) {
  if(all.equal(y[[i]],z[[i]])==TRUE)
    r<-r+1
}
length(y)/r

# Question 2: HMMs

# Build the HMM.

library(HMM)
#set.seed(123)

States<-1:100
Symbols<-1:2 # 1=door

transProbs<-matrix(rep(0,length(States)*length(States)), nrow=length(States), ncol=length(States), byrow = TRUE)
for(i in 1:99){
  transProbs[i,i]<-.1
  transProbs[i,i+1]<-.9
}

emissionProbs<-matrix(rep(0,length(States)*length(Symbols)), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
for(i in States){
  if(i %in% c(10,11,12,20,21,22,30,31,32)){
    emissionProbs[i,1]<-.9
    emissionProbs[i,2]<-.1
  }
  else{
    emissionProbs[i,1]<-.1
    emissionProbs[i,2]<-.9
  }
}

startProbs<-rep(1/100,100)
hmm<-initHMM(States,Symbols,startProbs,transProbs,emissionProbs)

# If the robot observes a door, it can be in front of any of the three doors. If it then observes a long
# sequence of non-doors, then it know that it was in front of the third door.

obs<-c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
pt<-prop.table(exp(forward(hmm,obs)),2)

which.maxima<-function(x){ # This function is needed since which.max only returns the first maximum.
  return(which(x==max(x)))
}

apply(pt,2,which.maxima)

# Question 3: GPs

### Question 3(a)
setwd('~/Desktop') # Change to your path
library(kernlab)
library(mvtnorm)

# From KernelCode.R: 
# Squared exponential, k
k <- function(sigmaf = 1, ell = 1)  
{   
  rval <- function(x, y = NULL) 
  {       
    r = sqrt(crossprod(x-y))       
    return(sigmaf^2*exp(-r^2/(2*ell^2)))     
  }   
  class(rval) <- "kernel"   
  return(rval) 
}  


# Simulating from the prior for ell = 0.2
kernel02 <- k(sigmaf = 1, ell = 0.2) # This constructs the covariance function
xGrid = seq(-1,1,by=0.1) 
K = kernelMatrix(kernel = kernel02, xGrid, xGrid)

colors = list("black","red","blue","green","purple")
f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
plot(xGrid,f, type = "l", ylim = c(-3,3), col = colors[[1]])
for (i in 1:4){
  f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
  lines(xGrid,f, col = colors[[i+1]])
}

# Simulating from the prior for ell = 1
kernel1 <- k(sigmaf = 1, ell = 1) # This constructs the covariance function
xGrid = seq(-1,1,by=0.1) 
K = kernelMatrix(kernel = kernel1, xGrid, xGrid)

colors = list("black","red","blue","green","purple")
f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
plot(xGrid,f, type = "l", ylim = c(-3,3), col = colors[[1]])
for (i in 1:4){
  f = rmvnorm(n = 1, mean = rep(0,length(xGrid)), sigma = K)
  lines(xGrid,f, col = colors[[i+1]])
}


# Computing the correlation functions

# ell = 0.2
kernel02(0,0.1) # Note: here correlation=covariance since sigmaf = 1 
kernel02(0,0.5)

# ell = 1
kernel1(0,0.1) # Note: here correlation=covariance since sigmaf = 1 
kernel1(0,0.5)

# The correlation between the function at x = 0 and x=0.1 is much higher than 
# between x = 0 and x=0.5, since the latter points are more distant.
# Thus, the correlation decays with distance in x-space. This decay is much more
# rapid when ell = 0.2 than when ell = 1. This is also visible in the simulations
# where realized functions are much more smooth when ell = 1.


### Question 3(b)

load("GPdata.RData")

### ell = 0.2

sigmaNoise = 0.2

# Set up the kernel function
kernelFunc <- k(sigmaf = 1, ell = 0.2)

# Plot the data and the true 
plot(x, y, main = "", cex = 0.5)

GPfit <- gausspr(x, y, kernel = kernelFunc, var = sigmaNoise^2)
# Alternative: GPfit <- gausspr(y ~ x, kernel = k, kpar = list(sigmaf = 1, ell = 0.2), var = sigmaNoise^2)
xs = seq(min(x),max(x), length.out = 100)
meanPred <- predict(GPfit, data.frame(x = xs)) # Predicting the training data. To plot the fit.
lines(xs, meanPred, col="blue", lwd = 2)

# Compute the covariance matrix Cov(f)
n <- length(x)
Kss <- kernelMatrix(kernel = kernelFunc, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = kernelFunc, x = x, y = x)
Kxs <- kernelMatrix(kernel = kernelFunc, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs)

# Probability intervals for f
lines(xs, meanPred - 1.96*sqrt(diag(Covf)), col = "red")
lines(xs, meanPred + 1.96*sqrt(diag(Covf)), col = "red")

# Prediction intervals for y 
lines(xs, meanPred - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")
lines(xs, meanPred + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")


legend("topright", inset = 0.02, legend = c("data","post mean","95% intervals for f", "95% predictive intervals for y"), 
       col = c("black", "blue", "red", "purple"), 
       pch = c('o',NA,NA,NA), lty = c(NA,1,1,1), lwd = 2, cex = 0.55)



### ell = 1

sigmaNoise = 0.2

# Set up the kernel function
kernelFunc <- k(sigmaf = 1, ell = 1)

# Plot the data and the true 
plot(x,y, main = "", cex = 0.5)
#lines(xGrid,fVals, type = "l", col = "black", lwd = 3) # true mean

GPfit <- gausspr(x, y, kernel = kernelFunc, var = sigmaNoise^2)
xs = seq(min(x), max(x), length.out = 100)
meanPred <- predict(GPfit, data.frame(x = xs)) # Predicting the training data. To plot the fit.
lines(xs, meanPred, col="blue", lwd = 2)

# Compute the covariance matrix Cov(f)
n <- length(x)
Kss <- kernelMatrix(kernel = kernelFunc, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = kernelFunc, x = x, y = x)
Kxs <- kernelMatrix(kernel = kernelFunc, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs)

# Probability intervals for f
lines(xs, meanPred - 1.96*sqrt(diag(Covf)), col = "red")
lines(xs, meanPred + 1.96*sqrt(diag(Covf)), col = "red")

# Prediction intervals for y 
lines(xs, meanPred - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")
lines(xs, meanPred + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "purple")


legend("topright", inset = 0.02, legend = c("data","post mean","95% intervals for f", "95% predictive intervals for y"), 
       col = c("black", "blue", "red", "purple"), 
       pch = c('o',NA,NA,NA), lty = c(NA,1,1,1), lwd = 2, cex = 0.55)


# Question: Explain the difference between the results from ii) and iii). 
# Answer: ii) is about the uncertainty of the function f, which is the MEAN of y
#         iii) is about the uncertainty of individual y values. They are uncertain for 
#              two reasons: you don't know f at the test point, and you don't know
#              the error (epsilon) that will hit this individual observation

# Question: Discuss the differences in results from using the two length scales.
#           Answer: shorter length scale gives less smooth f. We are overfitting the data.
#           Answer: longer length scale gives more smoothness.

# Question: Do you think a GP with a squared exponential kernel is a good model for this data? If not, why?
#           Answer: One would have to experiment with other length scales, or estimate
#           the length scales (see question 3c), but this is not likely to help here.
#           The issue is that the data seems to be have different smoothness for small x
#           than it has for large x (where the function seems much more flat)
#           The solution is probably to have different length scales for different x


### Question 3(c)

# For full points here you should mention EITHER of the following two approaches:
# 1. The marginal likelihood can be used to select optimal hyperparameters, 
# and also the noise variance. We can optimize the log marginal likelihood with
# respect to the hyperparameters. In Gaussian Process Regression the marginal likelihood
# is availble in closed form (a formula).
# 2. We can use sampling methods (e.g. MCMC) to sample from the marginal posterior of the hyperparameters
# We need a prior p(theta) for the hyperparameter and then Bayes rule gives the marginal posterior
# p(theta | data) propto p(data | theta)*p(theta)
# where p(data | theta) is the marginal likelihood (f has been integrated out).

# If the noise variance is unknown, we can treat like any of the kernel hyperparameters and infer the noise variance 
# jointly with the length scale and the prior variance sigma_f



# Question 4: SSMs

# Kalman filter implementation.

set.seed(12345)
start_time <- Sys.time()

T<-10000
mu_0<-50
Sigma_0<-10
R<-1
Q<-5

x<-vector(length=T)
z<-vector(length=T)
err<-vector(length=T)

for(t in 1:T){
  x[t]<-ifelse(t==1,rnorm(1,mu_0,Sigma_0),x[t-1]+1+rnorm(1,0,R))
  z[t]<-x[t]+rnorm(1,0,Q)
}

mu<-mu_0
Sigma<-Sigma_0*Sigma_0 # KF uses covariances
for(t in 2:T){
  pre_mu<-mu+1
  pre_Sigma<-Sigma+R*R # KF uses covariances
  K<-pre_Sigma/(pre_Sigma+Q*Q) # KF uses covariances
  mu<-pre_mu+K*(z[t]-pre_mu)
  Sigma<-(1-K)*pre_Sigma
  
  err[t]<-abs(x[t]-mu)
  
  cat("t: ",t,", x_t: ",x[t],", E[x_t]: ",mu," , error: ",err[t],"\n")
  flush.console()
}

mean(err[2:T])
sd(err[2:T])

end_time <- Sys.time()
end_time - start_time

# Repetition with the particle filter.

set.seed(12345)
start_time <- Sys.time()

T<-10000
n_par<-100
tra_sd<-1
emi_sd<-5
mu_0<-50
Sigma_0<-10

ini_dis<-function(n){
  return (rnorm(n,mu_0,Sigma_0))
}

tra_dis<-function(zt){
  return (rnorm(1,mean=zt+1,sd=tra_sd))
}

emi_dis<-function(zt){
  return (rnorm(1,mean=zt,sd=emi_sd))
}

den_emi_dis<-function(xt,zt){
  return (dnorm(xt,mean=zt,sd=emi_sd))
}

z<-vector(length=T)
x<-vector(length=T)

for(t in 1:T){
  z[t]<-ifelse(t==1,ini_dis(1),tra_dis(z[t-1]))
  x[t]<-emi_dis(z[t])
}

err<-vector(length=T)

bel<-ini_dis(n_par)
w<-rep(1/n_par,n_par)
for(t in 2:T){
  com<-sample(1:n_par,n_par,replace=TRUE,prob=w)
  bel<-sapply(bel[com],tra_dis)
  
  for(i in 1:n_par){
    w[i]<-den_emi_dis(x[t],bel[i])
  }
  w<-w/sum(w)
  
  Ezt<-sum(w * bel)
  err[t]<-abs(z[t]-Ezt)
  
  cat("t: ",t,", z_t: ",z[t],", E[z_t]: ",Ezt," , error: ",err[t],"\n")
  flush.console()
}

mean(err[2:T])
sd(err[2:T])

end_time <- Sys.time()
end_time - start_time

# KF works optimally (i.e. it computes the exact belief function in closed-form) since the SSM sampled is 
# linear-Gaussian. The particle filter on the other hand is approximate. The more particles the closer its 
# performance to the KF's but at the cost of increasing the running time.
