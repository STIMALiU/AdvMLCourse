# PGMs

#source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")

library(bnlearn)
library(gRain)
set.seed(567)
data("asia")
ind <- sample(1:5000, 4000)
tr <- asia[ind,]
te <- asia[-ind,]

nb0<-empty.graph(c("S","A","T","L","B","E","X","D"))
arc.set<-matrix(c("S", "A", "S", "T", "S", "L", "S", "B", "S", "E", "S", "X", "S", "D"), 
                ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))
arcs(nb0)<-arc.set
#nb<-naive.bayes(tr,"S")
plot(nb0)

totalError<-array(dim=6)
k<-1
for(j in c(10,20,50,100,1000,2000)){
  nb<-bn.fit(nb0,tr[1:j,],method="bayes")
  nb2<-as.grain(nb)
  nb2<-compile(nb2)

  error<-matrix(c(0,0,0,0),nrow=2,ncol=2)
  for(i in 1:1000){
    z<-NULL
    for(j in c("A","T","L","B","E","X","D")){
      if(te[i,j]=="no"){
        z<-c(z,"no")
      }
      else{
        z<-c(z,"yes")
      }
    }
  
    nb3<-setFinding(nb2,nodes=c("A","T","L","B","E","X","D"),states=z)
    x<-querygrain(nb3,c("S"))
  
    if(x$S[1]>x$S[2]){
      y<-1
    }
    else{
      y<-2
    }
  
    if(te[i,2]=="no"){
      error[y,1]<-error[y,1]+1
    }
    else{
      error[y,2]<-error[y,2]+1
    }
  }
  totalError[k]<-(error[1,2]+error[2,1])/1000
  k<-k+1
}
totalError

nb0<-empty.graph(c("S","A","T","L","B","E","X","D"))
arc.set<-matrix(c("A", "S", "T", "S", "L", "S", "B", "S", "E", "S", "X", "S", "D", "S"), 
                ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))
arcs(nb0)<-arc.set
#nb<-naive.bayes(tr,"S")
plot(nb0)

totalError<-array(dim=6)
k<-1
for(j in c(10,20,50,100,1000,2000)){
  nb<-bn.fit(nb0,tr[1:j,],method="bayes")
  nb2<-as.grain(nb)
  nb2<-compile(nb2)
  
  error<-matrix(c(0,0,0,0),nrow=2,ncol=2)
  for(i in 1:1000){
    z<-NULL
    for(j in c("A","T","L","B","E","X","D")){
      if(te[i,j]=="no"){
        z<-c(z,"no")
      }
      else{
        z<-c(z,"yes")
      }
    }
    
    nb3<-setFinding(nb2,nodes=c("A","T","L","B","E","X","D"),states=z)
    x<-querygrain(nb3,c("S"))
    
    if(x$S[1]>x$S[2]){
      y<-1
    }
    else{
      y<-2
    }
    
    if(te[i,2]=="no"){
      error[y,1]<-error[y,1]+1
    }
    else{
      error[y,2]<-error[y,2]+1
    }
  }
  totalError[k]<-(error[1,2]+error[2,1])/1000
  k<-k+1
}
totalError

# Discussion
# The NB classifier only needs to estimate the parameters for distributions of the form
# p(C) and P(A_i|C) where C is the class variable and A_i is a predictive attribute. The
# alternative model needs to estimate p(C) and P(C|A_1,...,A_n). Therefore, it requires
# more data to obtain reliable estimates (e.g. many of the parental combinations may not
# appear in the training data and this is actually why you should use method="bayes", to
# avoid zero relative frequency counters).This may hurt performance when little learning
# data is available. This is actually observed in the experiments above. However, when
# the size of the learning data increases, the alternative model should outperform NB, because
# the latter assumes that the attributes are independent given the class whereas the former
# does not. In other words, note that p(C|A_1,...,A_n) is proportional to P(A_1,...,A_n|C) p(C)
# by Bayes theorem. NB assumes that P(A_1,...,A_n|C) factorizes into a product of factors
# p(A_i|C) whereas the alternative model assumes nothing. The NB's assumption may hurt
# performance. This can be observed in the experiments.

# HMMs

library(HMM)
library(entropy)
set.seed(567)
States=1:10
Symbols=1:10
transProbs=matrix(c(.5,.5,0,0,0,0,0,0,0,0,
                    0,.5,.5,0,0,0,0,0,0,0,
                    0,0,.5,.5,0,0,0,0,0,0,
                    0,0,0,.5,.5,0,0,0,0,0,
                    0,0,0,0,.5,.5,0,0,0,0,
                    0,0,0,0,0,.5,.5,0,0,0,
                    0,0,0,0,0,0,.5,.5,0,0,
                    0,0,0,0,0,0,0,.5,.5,0,
                    0,0,0,0,0,0,0,0,.5,.5,
                    .5,0,0,0,0,0,0,0,0,.5), nrow=length(States), ncol=length(States), byrow = TRUE)
emissionProbs=matrix(c(.2,.2,.2,0,0,0,0,0,.2,.2,
                       .2,.2,.2,.2,0,0,0,0,0,.2,
                       .2,.2,.2,.2,.2,0,0,0,0,0,
                       0,.2,.2,.2,.2,.2,0,0,0,0,
                       0,0,.2,.2,.2,.2,.2,0,0,0,
                       0,0,0,.2,.2,.2,.2,.2,0,0,
                       0,0,0,0,.2,.2,.2,.2,.2,0,
                       0,0,0,0,0,.2,.2,.2,.2,.2,
                       .2,0,0,0,0,0,.2,.2,.2,.2,
                       .2,.2,0,0,0,0,0,.2,.2,.2), nrow=length(States), ncol=length(States), byrow = TRUE)
startProbs=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)
sim=simHMM(hmm,100)
logf=forward(hmm,sim$observation[1:100])
ef=exp(logf)
pt=prop.table(ef,2)
maxpt=apply(pt,2,which.max)
table(maxpt==sim$states)
post=posterior(hmm,sim$observation[1:100])
maxpost=apply(post,2,which.max)
table(maxpost==sim$states)

# Forward phase

a<-matrix(NA,nrow=100, ncol=length(States))
for(i in States){
  a[1,i]<-emissionProbs[sim$observation[1],i]*startProbs[i]
}

for(t in 2:100){
  for(i in States){
    a[t,i]<-emissionProbs[i,sim$observation[t]]*sum(a[t-1,]*transProbs[,i])
  }
}

for(t in 1:100){
  a[t,]<-a[t,]/sum(a[t,])
}

maxa=apply(a,1,which.max)
table(maxa==sim$states)

# Backward phase

b<-matrix(NA,nrow=100, ncol=length(States))
for(i in States){
  b[100,i]<-1
}

for(t in 99:1){
  for(i in States){
    b[t,i]<-sum(b[t+1,]*emissionProbs[,sim$observation[t+1]]*transProbs[i,])
  }
}

for(t in 1:100){
  for(i in States){
    b[t,i]<-b[t,i]*a[t,i]
  }
  b[t,]<-b[t,]/sum(b[t,])
}

maxb=apply(b,1,which.max)
table(maxb==sim$states)

# SSMs

set.seed(12345)
start_time <- Sys.time()

T<-100
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

# GPs

source('KernelCode.R') # Reading the Matern32 kernel from file

sigma2f = 1
ell = 0.5
zGrid <- seq(0.01, 1, by = 0.01)
count = 0
covs = rep(0,length(zGrid))
for (z in zGrid){
  count = count + 1
  covs[count] <- sigma2f*k(sigmaf = 1, ell = ell)(0,z)
}
plot(zGrid, covs, type = "l", xlab = "ell")

# The graph plots Cov(f(0),f(z)), the correlation between two FUNCTION VALUES, as a function of the distance between
# two inputs (0 and z)
# As expected the correlation between two points on f decreases as the distance increases. 
# The fact that points of f are dependent, as given by the covariance in the plot, makes the curves smooth. Nearby inputs will
# have nearby outputs when the correlation is large.

sigma2f = 0.5
ell = 0.5
zGrid <- seq(0.01, 1, by = 0.01)
count = 0
covs = rep(0,length(zGrid))
for (z in zGrid){
  count = count + 1
  covs[count] <- sigma2f*k(sigmaf = 1, ell = ell)(0,z)
}
plot(zGrid, covs, type = "l", xlab = "ell")

# Changing sigma2f will have not effect on the relative covariance between points on the curve, i.e. will not affect the
# smoothness. But lowering sigma2f makes the whole covariance curve lower. This means that the variance k(0,0) is lower and has the effect of giving a probability 
# distribution over curves which is tighter (lower variance) around the mean function of the GP. This means that simulated curves
# from the GP will be less variable. 

####################################
### GP inference with the ell = 1
####################################

library(kernlab)
load("lidar.RData") # loading the data
sigmaNoise = 0.05
x = distance
y = logratio

# Set up the kernel function
kernelFunc <- k(sigmaf = 1, ell = 1)

# Plot the data and the true 
plot(x, y, main = "", cex = 0.5)

GPfit <- gausspr(x, y, kernel = kernelFunc, var = sigmaNoise^2)
xs = seq(min(x),max(x), length.out = 100)
meanPred <- predict(GPfit, xs) # Predicting the training data. To plot the fit.
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

####################################
### GP inference with the ell = 5
####################################

load("lidar.RData") # loading the data
sigmaNoise = 0.05
x = distance
y = logratio

# Set up the kernel function
kernelFunc <- k(sigmaf = 1, ell = 5)

# Plot the data and the true 
plot(x, y, main = "", cex = 0.5)

GPfit <- gausspr(x, y, kernel = kernelFunc, var = sigmaNoise^2)
xs = seq(min(x),max(x), length.out = 100)
meanPred <- predict(GPfit, xs) # Predicting the training data. To plot the fit.
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


# Discussion
# The larger length scale gives smoother fits. The smaller length scale seems to generate too jagged fits. 
