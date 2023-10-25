# PGMs

library(bnlearn)
library(gRain)

net <- model2network("[U][C|U][A|C][B|C][D|A:B][Ch|U][Ah|Ch][Bh|Ch][Dh|Ah:Bh]")
plot(net)

cptU <- c(.5,.5)
dim(cptU) <- c(2)
dimnames(cptU) <- list(c("0", "1"))

cptC <- matrix(c(.9,.1,.1,.9), nrow=2, ncol=2)
dim(cptC) <- c(2,2)
dimnames(cptC) <- list("C" = c("0", "1"), "U" =  c("0", "1"))

cptA <- matrix(c(1,0,.2,.8), nrow=2, ncol=2)
dim(cptA) <- c(2,2)
dimnames(cptA) <- list("A" = c("0", "1"), "C" =  c("0", "1"))

cptB <- matrix(c(1,0,.2,.8), nrow=2, ncol=2)
dim(cptB) <- c(2,2)
dimnames(cptB) <- list("B" = c("0", "1"), "C" =  c("0", "1"))

cptD <- matrix(c(.9,.1,0,1,0,1,0,1), nrow=2, ncol=4)
dim(cptD) <- c(2,2,2)
dimnames(cptD) <- list("D" = c("0", "1"), "A" =  c("0", "1"), "B" =  c("0", "1"))

cptCh <- matrix(c(.9,.1,.1,.9), nrow=2, ncol=2)
dim(cptCh) <- c(2,2)
dimnames(cptCh) <- list("Ch" = c("0", "1"), "U" =  c("0", "1"))

cptAh <- matrix(c(1,0,.2,.8), nrow=2, ncol=2)
dim(cptAh) <- c(2,2)
dimnames(cptAh) <- list("Ah" = c("0", "1"), "Ch" =  c("0", "1"))

cptBh <- matrix(c(1,0,.2,.8), nrow=2, ncol=2)
dim(cptBh) <- c(2,2)
dimnames(cptBh) <- list("Bh" = c("0", "1"), "Ch" =  c("0", "1"))

cptDh <- matrix(c(.9,.1,0,1,0,1,0,1), nrow=2, ncol=4)
dim(cptDh) <- c(2,2,2)
dimnames(cptDh) <- list("Dh" = c("0", "1"), "Ah" =  c("0", "1"), "Bh" =  c("0", "1"))

netfit <- custom.fit(net,list(U=cptU, C=cptC, A=cptA, B=cptB, D=cptD, Ch=cptCh, Ah=cptAh, Bh=cptBh, Dh=cptDh))
netcom <- compile(as.grain(netfit))

querygrain(setEvidence(netcom,nodes=c("D","Ah"),states=c("1","0")),c("Dh"))

# HMMs

library(HMM)

States=c("Z=1.C=2","Z=1.C=1","Z=2.C=3","Z=2.C=2","Z=2.C=1","Z=3.C=2","Z=3.C=1","Z=4.C=1","Z=5.C=2","Z=5.C=1")
Symbols=c("1","2","3","4","5")
transProbs=matrix(c(0,1,0,0,0,0,0,0,0,0,
                    0,.5,.5,0,0,0,0,0,0,0,
                    0,0,0,1,0,0,0,0,0,0,
                    0,0,0,0,1,0,0,0,0,0,
                    0,0,0,0,.5,.5,0,0,0,0,
                    0,0,0,0,0,0,1,0,0,0,
                    0,0,0,0,0,0,.5,.5,0,0,
                    0,0,0,0,0,0,0,.5,.5,0,
                    0,0,0,0,0,0,0,0,0,1,
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
startProbs=c(.2,0,.2,0,0,.2,0,.2,.2,0)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)
sim=simHMM(hmm,100)
sim

# \begin{align*}
# p(z^{T-1}|x^{0:T}) &=\sum_{z^T} p(z^{T-1},z^T|x^{0:T})\\
# &=\sum_{z^T} p(z^{T-1}|z^T,x^{0:T}) p(z^T|x^{0:T})\\
# &=\sum_{z^T} p(z^{T-1}|z^T,x^{0:T-1}) p(z^T|x^{0:T})\\
# &=\sum_{z^T} \frac{p(z^T|z^{T-1},x^{0:T-1}) p(z^{T-1}|x^{0:T-1}) p(z^T|x^{0:T})}{p(z^T|x^{0:T-1})}\\
# &=\sum_{z^T} \frac{p(z^T|z^{T-1}) p(z^{T-1}|x^{0:T-1}) p(z^T|x^{0:T})}{p(z^T|x^{0:T-1})}\\
# &=\sum_{z^T} \frac{p(z^T|z^{T-1}) p(z^{T-1}|x^{0:T-1}) p(z^T|x^{0:T})}{\sum_{z^{T-1}} p(z^{T-1},z^T|x^{0:T-1})}\\
# &=\sum_{z^T} \frac{p(z^T|z^{T-1}) p(z^{T-1}|x^{0:T-1}) p(z^T|x^{0:T})}{\sum_{z^{T-1}} p(z^T|z^{T-1},x^{0:T-1}) p(z^{T-1}|x^{0:T-1})}\\
# &=\sum_{z^T} \frac{p(z^T|z^{T-1}) p(z^{T-1}|x^{0:T-1}) p(z^T|x^{0:T})}{\sum_{z^{T-1}} p(z^T|z^{T-1}) p(z^{T-1}|x^{0:T-1})}
# \end{align*}

# RL

theta <- 0.1
gamma <- .95
V <- array(0,dim = 10)
pi <- array(0,dim = 10)

repeat{
  delta <- 0
  for(s in 1:9){
    v <- V[s]
    V[s] <- max(gamma*V[s],(s+1==10)+gamma*V[s+1])
    delta <- max(delta,abs(v-V[s]))
  }
  if(delta<theta) break
}

for(s in 1:9){
  pi[s] <- which.max(c(gamma*V[s],(s+1==10)+gamma*V[s+1]))
}
V
pi

# GPs

tempData <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv', header=TRUE, sep=';')
temp <- tempData$temp
time = 1:length(temp)
subset <- seq(1, length(temp), by = 5)
temp <- temp[subset]
time = time[subset]

LM <- function(par=c(20,0.2),X,y,k,sigmaNoise){
  n <- length(y)
  L <- t(chol(k(par,X,X)+((sigmaNoise^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  logmar <- -0.5*(t(y)%*%a)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(logmar)
}

polyFit <- lm(temp ~  time + I(time^2))
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

SEKernel2(c(20,100),c(1,182,365),c(1,182,365))

# The further apart two points are in x values, the less correlated their f values are.

LM(par=c(20,100),X=time,y=temp,k=SEKernel2,sigmaNoise=sigmaNoiseFit)

LocallyPeriodicSine <- function(par=c(1,1,1,1), x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n1){
    for (j in 1:n2){
      r = sqrt(crossprod(x1[i]-x2[j]))
      K[i,j] <- par[1]^2*exp(-2*(sin(pi*r/par[4])^2)/par[2]^2)*exp(-0.5*r^2/par[3]^2)
    }
  }
  return(K)
} 

LocallyPeriodicSine(c(20,1,1000,365),c(1,182,365),c(1,182,365))

# The further apart two points are in x values MODULE half-a-year, the less correlated their f values are.
# Note that I used l_2=1000 while the exam says l_2=100. My original intention was to have l_2=1000 but
# unfortunately I wrote l_2=100 in the exam. I will correct with l_2=100 then, in which case you do not see
# the "module half-a-year" effect. That is a pity. My bad.

LM(par=c(20,1,1000,365),X=time,y=temp,k=LocallyPeriodicSine,sigmaNoise=sigmaNoiseFit)

# Using a validation set instead of the log marginal is tricky, because the data is not iid.
# You cannot just simply select the last points as validation data, because they do not come
# from the same distribution as the training data (their time values are larger). You cannot
# either randomly select the validation data, because you will be using the future to predict
# the past, which may result in overfitting. Although computationally demanding, one solution
# is to train a GP with the first n points and use it to predict the point n+1, then 
# train a GP with the first n+1 points and use it to predict the point n+2, and so on.
