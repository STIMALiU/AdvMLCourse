# BNs

library(bnlearn)
library(gRain)

set.seed(123)
data("asia")
tr<-asia
hc<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Asia DAG
hc<-bn.fit(hc,tr[1:10,],method="bayes") # Learn the parameters
hc$D
hc2<-as.grain(hc) # Convert to gRain object
hc2<-compile(hc2) # Compile for inference

for(i in 11:5000){
  z<-NULL
  for(j in c("A","S","T","L","X","D")){
    if(tr[i,j]=="no"){
      z<-c(z,"no")
    }
    else{
      z<-c(z,"yes")
    }
  }
  
  hc3<-setEvidence(hc2,nodes=c("A","S","T","L","X","D"),states=z) # Enter the evidence
  b<-querygrain(hc3,c("B")) # Get posterior distribution
  tr[i,"B"]<-sample(c("no","yes"),size=1,prob=b$B)
  
  hc3<-setEvidence(hc2,nodes=c("A","S","T","L","X","D"),states=z) # Enter the evidence
  e<-querygrain(hc3,c("E")) # Get posterior distribution
  tr[i,"E"]<-sample(c("no","yes"),size=1,prob=e$E)
}

hc<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Asia DAG
hc<-bn.fit(hc,tr,method="bayes") # Learn the parameters
hc$D

hc<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]") # True Asia DAG
hc<-bn.fit(hc,asia,method="bayes") # "True" Asia parameters
hc$D

# The parameters learned from the completed dataset are closer to the "true" ones.

# HMMs

library(HMM)

States=c("1a","1b","2a","2b","2c","3a","3b","4","5a","5b")
Symbols=c("1","2","3","4","5")
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
emissionProbs=matrix(c(1/3,1/3,0,0,1/3,
                       1/3,1/3,0,0,1/3,
                       1/3,1/3,1/3,0,0,
                       1/3,1/3,1/3,0,0,
                       1/3,1/3,1/3,0,0,
                       0,1/3,1/3,1/3,0,
                       0,1/3,1/3,1/3,0,
                       0,0,1/3,1/3,1/3,
                       1/3,0,0,1/3,1/3,
                       1/3,0,0,1/3,1/3), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
startProbs=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)
sim=simHMM(hmm,100)
sim

# RL

# The value of alpha affects the convergence of the algorithm. When alpha=0.1, an optimal
# path is typically found. The q-values corresponding to the actions in the path are close
# to 10, which is the true value as gamma=1. For alpha=0.01, 0.001, an optimal path may still
# be found but the q-values for the actions in it have not converged.

# GPs (I)

X<-seq(0,10,.1)
Yfun<-function(x){
  return (x*(sin(x)+sin(3*x))+rnorm(length(x),0,2))
}
plot(X,Yfun(X),xlim=c(0,10),ylim=c(-15,15))

SEKernel <- function(x1,x2){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- (sigmaF^2)*exp(-0.5*( (x1-x2[i])/l)^2 )
  }
  return(K)
}

posteriorGP <- function(X,y,k,sigmaNoise,xStar){
  n <- length(y)
  L <- t(chol(k(X,X)+((sigmaNoise^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  kStar <- k(X,xStar)
  mu <- t(kStar)%*%a
  v <- solve(L,kStar)
  var <- k(xStar,xStar)-(t(v)%*%v)
  logmar <- -0.5*(t(y)%*%a)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(list("mu"=mu,"var"=var,"logmar"=logmar))
}

sigmaF <- 1.5
l <- .5 # These hyperparameter values work well as the y probability bands cover around 95% of the y values in the dataset.
xData <- X
yData <- Yfun(X)
sigmaN <- 2
xGrid <- seq(0,10,.1)
res<-posteriorGP(X=xData,y=yData,k=SEKernel,sigmaNoise=sigmaN,xStar=xGrid)
plot(xData,yData,xlim=c(0,10),ylim=c(-15,15))
lines(xGrid, res$mu, col="blue", lwd = 2)
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)

# GPs (II)

X<-seq(0,10,2)
Yfun<-function(x){
  return (x*(sin(x)+sin(3*x))+rnorm(length(x),0,.2))
}
plot(X,Yfun(X),xlim=c(0,10),ylim=c(-15,15))

for(i in 1:4){
sigmaF <- 1.5
l <- .5 # These hyperparameter values work well as the y probability bands cover around 95% of the y values in the dataset.
xData <- X
yData <- Yfun(X)
sigmaN <- .2
xGrid <- seq(0,10,.1)
res<-posteriorGP(X=xData,y=yData,k=SEKernel,sigmaNoise=sigmaN,xStar=xGrid)
plot(xData,yData,xlim=c(0,10),ylim=c(-15,15))
lines(xGrid, res$mu, col="blue", lwd = 2)
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)), col = "red")
lines(xGrid, res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)
lines(xGrid, res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2), col = "red", lwd=2)

foo<-.1*which.max(res$mu + 1.96*sqrt(diag(res$var)+sigmaN^2) - (res$mu - 1.96*sqrt(diag(res$var)+sigmaN^2)))
X<-c(X,foo)
}
