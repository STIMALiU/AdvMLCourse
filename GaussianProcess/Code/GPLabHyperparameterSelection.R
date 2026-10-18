# Solutions to the lab 4 of 732A96/TDDE15 Advanced Machine Learning
# Author: jose.m.pena@liu.se
# Made for teaching purposes

# Question 2.2: Tullinge temperatures

# Bonus: Hyperparameter selection via log marginal likelihood maximization

LM <- function(par=c(20,100),X,y,k,sigmaNoise){
  
  # Returns the log marginal likelihood. See Algorithm 2.1
  
}

# Grid search

bestLM<-LM(par=c(20,100),X=time,y=temp,k=SEKernel2,sigmaNoise=sigmaNoiseFit)
bestLM
besti<-20
bestj<-100
for(i in seq(1,50,1)){
  for(j in seq(50,150,5)){
    aux<-LM(par=c(i,j),X=time,y=temp,k=SEKernel2,sigmaNoise=sigmaNoiseFit)
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

# optim function

foo<-optim(par = c(1,1), fn = LM, X=time,y=temp,k=SEKernel2,sigmaNoise=sigmaNoiseFit, method="L-BFGS-B",
      lower = c(.Machine$double.eps, .Machine$double.eps),control=list(fnscale=-1))

foo$value
sigmaF <- foo$par[1]
l <- foo$par[2]
sigmaF
l
res<-posteriorGP(X=time,y=temp,k=SEKernel,sigmaNoise=sigmaNoiseFit,xStar=time)
lines(time, res$mu, col="blue", lwd = 2)
lines(time, res$mu - 1.96*sqrt(diag(res$var)), col = "black")
lines(time, res$mu + 1.96*sqrt(diag(res$var)), col = "black")

# Question 2.3: Banknote fraud

# Bonus: Hyperparameter selection via validation set

set.seed(1)
SelectTraining <- sample(1:dim(data)[1], size = 1000)
y <- data[,5]
X <- as.matrix(data[,1:4])
yTrainVal <- y[SelectTraining]
yTest <- y[-SelectTraining]
XTrainVal <- X[SelectTraining,]
XTest <- X[-SelectTraining,]
SelectVal <- sample(1:1000, size = 200, replace = FALSE)
yVal <- yTrainVal[SelectVal]
XVal <- XTrainVal[SelectVal,]
yTrain <- yTrainVal[-SelectVal]
XTrain <- XTrainVal[-SelectVal,]

selVars <- c(1,2,3,4)

acVal <- function(par=c(0.1)){
  GPfitFraud <- # train GP on XTrain using hyperparameter value par
  predVal <- predict(GPfitFraud,XVal[,selVars])
  table(predVal, yVal)
  accuracyVal <-sum(predVal==yVal)/length(yVal)
  return(accuracyVal)
}

# Grid search

bestVal<-0
bestj<-GPfitFraud@kernelf@kpar$sigma
for(j in seq(0.1,10,0.1)){
  aux <- acVal(j)
  if(bestVal<aux){
    bestVal<-aux
    bestj<-j
  }
}
bestVal
bestj

GPfitFraud <- # train GP on XTrainVal using hyperparameter value bestj
predTest <- predict(GPfitFraud,XTest[,selVars])
table(predTest, yTest)
sum(predTest==yTest)/length(yTest)

# optim function

foo<-optim(par = c(0.1), fn = acVal, method="L-BFGS-B",
           lower = c(.Machine$double.eps),control=list(fnscale=-1))

acVal(foo$par)
foo$par

GPfitFraud <- # train GP on XTrainVal using hyperparameter value foo$par
predTest <- predict(GPfitFraud,XTest[,selVars])
table(predTest, yTest)
sum(predTest==yTest)/length(yTest)

        