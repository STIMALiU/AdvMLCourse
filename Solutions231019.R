# Solutions to the exam 23/10/2019 of 732A96/TDDE15 Advanced Machine Learning
# Author: jose.m.pena@liu.se
# Made for teaching purposes

# Question 1.1 (Monty Hall)

library(bnlearn)
library(gRain)

MH<-model2network("[D][P][M|D:P]") # Structure
cptD = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("D1", "D2", "D3"))) # Parameters
cptP = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("P1", "P2", "P3")))
cptM = c(
  0,.5,.5,
  0,0,1,
  0,1,0,
  0,0,1,
  .5,0,.5,
  1,0,0,
  0,1,0,
  1,0,0,
  .5,.5,0)
dim(cptM) = c(3, 3, 3)
dimnames(cptM) = list("M" = c("M1", "M2", "M3"), "D" =  c("D1", "D2", "D3"), "P" = c("P1", "P2", "P3"))
MHfit<-custom.fit(MH,list(D=cptD,P=cptP,M=cptM))
MHcom<-compile(as.grain(MHfit))

MHfitEv<-setFinding(MHcom,nodes=c(""),states=c("")) # Exact inference
querygrain(MHfitEv,c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","M"),states=c("D1","M2"))
querygrain(MHfitEv,c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","M"),states=c("D1","M3"))
querygrain(MHfitEv,c("P"))

table(cpdist(MHfit,"P",evidence=TRUE)) # Alternatively, one can use approximate inference
table(cpdist(MHfit,"P",(D=="D1" & M=="M2")))
table(cpdist(MHfit,"P",(D=="D1" & M=="M3")))

# Question 1.2 (XOR)

library(bnlearn)
library(gRain)

xor<-model2network("[A][B][C|A:B]") # Structure
cptA = matrix(c(0.5, 0.5), ncol = 2, dimnames = list(NULL, c("0", "1"))) # Parameters
cptB = matrix(c(0.5, 0.5), ncol = 2, dimnames = list(NULL, c("0", "1")))
cptC = c(1,0,0,1,0,1,1,0)
dim(cptC) = c(2, 2, 2)
dimnames(cptC) = list("C" = c("0", "1"), "A" =  c("0", "1"), "B" = c("0", "1"))
xorfit<-custom.fit(xor,list(A=cptA,B=cptB,C=cptC))

for(i in 1:10){
  xordata<-rbn(xorfit,1000) # Sample
  xorhc<-hc(xordata,score="bic") # Learn
  plot(xorhc)
}

# The HC fails because the data comes from a distribution that is not faithful to the graph.
# In other words, the graph has an edge A -> C but A is marginally independent of C in the
# distribution. The same for B and C. And, of course, the same for A and B since there is no 
# edge between them. So, the HC cannot find two dependent variables to link with an edge and,
# thus, it adds no edge to the initial empty graph.

# Question 2 (HMMs)

library(HMM)

States=1:10 # Sectors
Symbols=1:11 # Sectors + malfunctioning
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
emissionProbs=matrix(c(.1,.1,.1,0,0,0,0,0,.1,.1,.5,
                       .1,.1,.1,.1,0,0,0,0,0,.1,.5,
                       .1,.1,.1,.1,.1,0,0,0,0,0,.5,
                       0,.1,.1,.1,.1,.1,0,0,0,0,.5,
                       0,0,.1,.1,.1,.1,.1,0,0,0,.5,
                       0,0,0,.1,.1,.1,.1,.1,0,0,.5,
                       0,0,0,0,.1,.1,.1,.1,.1,0,.5,
                       0,0,0,0,0,.1,.1,.1,.1,.1,.5,
                       .1,0,0,0,0,0,.1,.1,.1,.1,.5,
                       .1,.1,0,0,0,0,0,.1,.1,.1,.5), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
startProbs=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)

obs=c(1,11,11,11)
posterior(hmm,obs)
viterbi(hmm,obs)

# Question 3 (SSMs)

T<-100
mu_0<-50
Sigma_0<-sqrt(100^2/12)
n_par<-100
tra_sd<-.1
emi_sd<-1

ini_dis<-function(n){ # Sampling the initial model
  return (runif(n,min=0,max=100))
}

tra_dis<-function(zt){ # Sampling the transition model
  pi=sample(c(0,1,2),1)
  return (rnorm(1,mean=zt+pi,sd=tra_sd))
}

emi_dis<-function(zt){ # Sampling the emision model
  pi=sample(c(-1,0,1),1)
  return (rnorm(1,mean=zt+pi,sd=emi_sd))
}

den_emi_dis<-function(xt,zt){ # Density of the emission model
  return (((dnorm(xt,mean=zt,sd=emi_sd)
            +dnorm(xt,mean=zt-1,sd=emi_sd)
            +dnorm(xt,mean=zt+1,sd=emi_sd))/3))
}

z<-vector(length=T) # Hidden states
x<-vector(length=T) # Observations

for(t in 1:T){ # Sampling the SSM
  z[t]<-ifelse(t==1,ini_dis(1),tra_dis(z[t-1]))
  x[t]<-emi_dis(z[t])
}

plot(z,col="red",main="True location (red) and observed location (black)")
points(x)

# Particle filter starts

bel<-ini_dis(n_par) # Initial particles
err<-vector(length=T) # Error between expected and true locations
w<-vector(length=n_par) # Importance weights

cat("Initial error: ",abs(z[1]-mean(bel)),"\n")

for(t in 1:T){
  for(i in 1:n_par){
    w[i]<-den_emi_dis(x[t],bel[i])
  }
  w<-w/sum(w)
  
  Ezt<-sum(w * bel) # Expected location
  err[t]<-abs(z[t]-Ezt) # Error between expected and true locations
  
  com<-sample(1:n_par,n_par,replace=TRUE,prob=w) # Sampling new particles according to the importance weights
  bel<-sapply(bel[com],tra_dis) # Project
}

mean(err[(T-10):T])
sd(err[(T-10):T])

# Kalman filter starts

R<-tra_sd
Q<-emi_sd

err<-vector(length=T)

mu<-mu_0
Sigma<-Sigma_0*Sigma_0 # KF uses covariances
for(t in 2:T){
  pre_mu<-mu+1
  pre_Sigma<-Sigma+R*R # KF uses covariances
  K<-pre_Sigma/(pre_Sigma+Q*Q) # KF uses covariances
  mu<-pre_mu+K*(x[t]-pre_mu)
  Sigma<-(1-K)*pre_Sigma
  
  err[t]<-abs(z[t]-mu)
}

mean(err[(T-10):T])
sd(err[(T-10):T])

# PF should outperform KF because the model is not liner-Gaussian, i.e. it is a mixture model.

# Question 4.1 (Hyperparameter selection via log marginal maximization)

tempData <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv', header=TRUE, sep=';')
temp <- tempData$temp
time = 1:length(temp)

# Extract every 5:th observation
subset <- seq(1, length(temp), by = 5)
temp <- temp[subset]
time = time[subset]

polyFit <- lm(scale(temp) ~  scale(time) + I(scale(time)^2))
sigmaNoiseFit = sd(polyFit$residuals)

SEKernel2 <- function(par=c(20,0.2),x1,x2){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- (par[1]^2)*exp(-0.5*( (x1-x2[i])/par[2])^2 )
  }
  return(K)
}

LM <- function(par=c(20,0.2),X,y,k,sigmaNoise){
  n <- length(y)
  L <- t(chol(k(par,X,X)+((sigmaNoise^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  logmar <- -0.5*(t(y)%*%a)-sum(diag(L))-(n/2)*log(2*pi)
  return(logmar)
}

bestLM<-LM(par=c(20,0.2),X=scale(time),y=scale(temp),k=SEKernel2,sigmaNoise=sigmaNoiseFit) # Grid search
bestLM
besti<-20
bestj<-0.2
for(i in seq(1,50,1)){
  for(j in seq(0.1,10,0.1)){
    aux<-LM(par=c(i,j),X=scale(time),y=scale(temp),k=SEKernel2,sigmaNoise=sigmaNoiseFit)
    if(bestLM<aux){
      bestLM<-aux
      besti<-i
      bestj<-j
    }
  }
}
bestLM
besti
bestj

foo<-optim(par = c(1,0.1), fn = LM, X=scale(time),y=scale(temp),k=SEKernel2,sigmaNoise=sigmaNoiseFit, method="L-BFGS-B",
           lower = c(.Machine$double.eps, .Machine$double.eps),control=list(fnscale=-1)) # Alternatively, one can use optim

foo$value
foo$par[1]
foo$par[2]

# Question 4.2 (Hyperparameter selection via validation set)

library(kernlab)

data <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv', header=FALSE, sep=',')
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111); SelectTraining <- sample(1:dim(data)[1], size = 1000, replace = FALSE)
y <- data[,5]
X <- as.matrix(data[,1:4])
yTrain <- y[SelectTraining]
yTest <- y[-SelectTraining]
XTrain <- X[SelectTraining,]
XTest <- X[-SelectTraining,]
SelectVal <- sample(1:1000, size = 200, replace = FALSE) # 800 samples for training, 200 for validation, and the rest for test (optional)
yVal <- yTrain[SelectVal]
XVal <- XTrain[SelectVal,]
yTrain <- yTrain[-SelectVal]
XTrain <- XTrain[-SelectVal,]

acVal <- function(par=c(0.1)){ # Accuracy on the validation set
  gausspr(x = XTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = list(sigma=par[1]))
  predVal <- predict(GPfitFraud,XVal[,selVars])
  table(predVal, yVal)
  accuracyVal <-sum(predVal==yVal)/length(yVal)
  return(accuracyVal)
}

selVars <- c(1,2,3,4)
GPfitFraud <- gausspr(x = XTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = 'automatic')
GPfitFraud
predVal <- predict(GPfitFraud,XVal[,selVars])
table(predVal, yVal)
accuracyVal <-sum(predVal==yVal)/length(yVal)
accuracyVal

bestVal<-accuracyVal # Grid search
for(j in seq(0.1,10,0.1)){
  aux <- acVal(j)
  if(bestVal<aux){
    bestVal<-aux
    bestj<-j
  }
}
bestVal
bestj

gausspr(x = XTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = list(sigma=bestj)) # Accuracy on the test set (optional)
predTest <- predict(GPfitFraud,XTest[,selVars])
table(predTest, yTest)
sum(predTest==yTest)/length(yTest)

foo<-optim(par = c(0.1), fn = acVal, method="L-BFGS-B",
           lower = c(.Machine$double.eps),control=list(fnscale=-1)) # Alternatively, one can use optim

acVal(foo$par)
foo$par

gausspr(x = XTrain[,selVars], y = yTrain, kernel = "rbfdot", kpar = list(sigma=foo$par)) # Accuracy on the test set (optional)
predTest <- predict(GPfitFraud,XTest[,selVars])
table(predTest, yTest)
sum(predTest==yTest)/length(yTest)