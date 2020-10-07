##############################################################################################
## Kernlab demo for Gaussian processes regression and classfication.
## Author: Mattias Villani, Linkoping University. http://mattiasvillani.com
##############################################################################################


#############################################################################
###    Prelims: Setting path and installing/loading packages            #####
#############################################################################

#install.packages('kernlab')
#install.packages("AtmRay") # To make 2D grid like in Matlab's meshgrid.
library(kernlab)
library(AtmRay)


###############################################
###    Messin' around with kernels        #####
###############################################

# This is just to test how one evaluates a kernel function
# and how one computes the covariance matrix from a kernel function.
X <- matrix(rnorm(12), 4, 3) # Simulating some data.
Xstar <- matrix(rnorm(15), 5, 3)
ell <- 1
SEkernel <- rbfdot(sigma = 1/(2*ell^2)) # Note how I reparametrize the rbfdot (which is the SE kernel) in kernlab.
SEkernel(1,2) # Just a test - evaluating the kernel in the points x=1 and x'=2.
# Computing the whole covariance matrix K from the kernel. Just a test.
kernelMatrix(kernel = SEkernel, x = X, y = Xstar) # So this is K(X,Xstar).

# Own implementation of Matern with nu = 3/2 (See RW book equation 4.17).
# Note that a call of the form kernelFunc <- Matern32(sigmaf = 1, ell = 0.1) returns a kernel FUNCTION.
# You can now evaluate the kernel at inputs: kernelFunc(x = 3, y = 4).
# Note also that class(kernelFunc) is of class "kernel", which is a class defined by kernlab.
Matern32 <- function(sigmaf = 1, ell = 1) 
{
  rval <- function(x, y = NULL) {
      r = sqrt(crossprod(x-y));
      return(sigmaf^2*(1+sqrt(3)*r/ell)*exp(-sqrt(3)*r/ell))
    }
  class(rval) <- "kernel"
  return(rval)
} 

# Testing our own defined kernel function.
X <- matrix(rnorm(12), 4, 3) # Simulating some data.
Xstar <- matrix(rnorm(15), 5, 3)
MaternFunc = Matern32(sigmaf = 1, ell = 2) # MaternFunc is a kernel FUNCTION.
MaternFunc(c(1,1),c(2,2)) # Evaluating the kernel in x=c(1,1), x'=c(2,2).
# Computing the whole covariance matrix K from the kernel.
kernelMatrix(kernel = MaternFunc, x = X, y = Xstar) # So this is K(X,Xstar).


##########################################
# Author: Jose M. Peña, jose.m.pena@liu.se
# GP regression on the canadian wages data
##########################################

library(kernlab)
CWData <- read.table('https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/CanadianWages.dat', 
                     header = T)
logWage<-CWData$logWage
age<-CWData$age
age<-(age-mean(age))/sd(age) # Standarize the age.

# Estimating the noise variance from a third degree polynomial fit. I() is needed because, otherwise
# age^2 reduces to age in the formula, i.e. age^2 means adding the main effect and the second order
# interaction, which in this case do not exist. See ?I.
polyFit <- lm(logWage ~  age + I(age^2) + I(age^3))
sigmaNoise = sd(polyFit$residuals)
plot(age,logWage)

# Fit the GP with built-in square expontial kernel (called rbfdot in kernlab).
ell <- 0.5
SEkernel <- rbfdot(sigma = 1/(2*ell^2)) # Note the reparametrization.
GPfit <- gausspr(age,logWage, kernel = SEkernel, var = sigmaNoise^2)
meanPred <- predict(GPfit, age) # Predicting the training data.
lines(age, meanPred, col="red", lwd = 2)

# The implementation of kernlab for the probability and prediction intervals seem to have a bug: The intervals
# seem to be too wide, e.g. replace 1.96 with 0.1 to see something.
# GPfit <- gausspr(age,logWage, kernel = SEkernel, var = sigmaNoise^2, variance.model = TRUE)
# meanPred <- predict(GPfit, age)
# lines(age, meanPred, col="red", lwd = 2)
# lines(age, meanPred+1.96*predict(GPfit,age, type="sdeviation"),col="blue")
# lines(age, meanPred-1.96*predict(GPfit,age, type="sdeviation"),col="blue")

# Probability and prediction interval implementation.
x<-age
xs<-age # XStar.
n <- length(x)
Kss <- kernelMatrix(kernel = SEkernel, x = xs, y = xs)
Kxx <- kernelMatrix(kernel = SEkernel, x = x, y = x)
Kxs <- kernelMatrix(kernel = SEkernel, x = x, y = xs)
Covf = Kss-t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), Kxs) # Covariance matrix of fStar.

# Probability intervals for fStar.
lines(xs, meanPred - 1.96*sqrt(diag(Covf)), col = "blue", lwd = 2)
lines(xs, meanPred + 1.96*sqrt(diag(Covf)), col = "blue", lwd = 2)

# Prediction intervals for yStar.
lines(xs, meanPred - 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "blue")
lines(xs, meanPred + 1.96*sqrt((diag(Covf) + sigmaNoise^2)), col = "blue")


###############################################
###      Regression on the LIDAR data       ###
###############################################

lidarData <- read.table('https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/LidarData', 
                        header = T)
LogRatio <- lidarData$LogRatio
Distance <- lidarData$Distance

# Estimating the noise variance from a third degree polynomial fit
polyFit <- lm(LogRatio ~  Distance + I(Distance^2) + I(Distance^3) )
sigmaNoise = sd(polyFit$residuals)

plot(Distance,LogRatio)

# Fit the GP with built in Square expontial kernel (called rbfdot in kernlab)
ell <- 2
GPfit <- gausspr(Distance, LogRatio, kernel = rbfdot, kpar = list(sigma = 1/(2*ell^2)), var = sigmaNoise^2)
meanPred <- predict(GPfit, Distance) # Predicting the training data. To plot the fit.
lines(Distance, meanPred, col="blue", lwd = 2)

# Fit the GP with home made Matern
sigmaf <- 1
ell <- 2
# GPfit <- gausspr(Distance, LogRatio, kernel = Matern32(ell=1)) # NOTE: this also works and is the same as the next line.
GPfit <- gausspr(Distance, LogRatio, kernel = Matern32, kpar = list(sigmaf = sigmaf, ell=ell), var = sigmaNoise^2) 
meanPred <- predict(GPfit, Distance)
lines(Distance, meanPred, col="purple", lwd = 2)

# Trying another length scale
sigmaf <- 1
ell <- 1
GPfit <- gausspr(Distance, LogRatio, kernel = Matern32, kpar = list(sigmaf = sigmaf, ell=ell), var = sigmaNoise^2) 
meanPred <- predict(GPfit, Distance)
lines(Distance, meanPred, col="green", lwd = 2)

# And now with a different sigmaf
sigmaf <- 0.1
ell <- 2
GPfit <- gausspr(Distance, LogRatio, kernel = Matern32, kpar = list(sigmaf = sigmaf, ell=ell), var = sigmaNoise^2) 
meanPred <- predict(GPfit, Distance)
lines(Distance, meanPred, col="black", lwd = 2)


###############################################
####       Classification on Iris data      ###
###############################################

data(iris)
GPfitIris <- gausspr(Species ~  Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data=iris)
GPfitIris

# predict on the training set
predict(GPfitIris,iris[,1:4])
table(predict(GPfitIris,iris[,1:4]), iris[,5]) # confusion matrix

# Now using only Sepal.Length and Sepal.Width to classify
GPfitIris <- gausspr(Species ~  Sepal.Length + Sepal.Width, data=iris)
GPfitIris
# predict on the training set
predict(GPfitIris,iris[,1:2])
table(predict(GPfitIris,iris[,1:2]), iris[,5]) # confusion matrix

# Now using only  Petal.Length + Petal.Width to classify
GPfitIris <- gausspr(Species ~ Petal.Length + Petal.Width, data=iris)
GPfitIris
# predict on the training set
predict(GPfitIris,iris[,3:4])
table(predict(GPfitIris,iris[,3:4]), iris[,5]) # confusion matrix

# class probabilities 
probPreds <- predict(GPfitIris, iris[,3:4], type="probabilities")
x1 <- seq(min(iris[,3]),max(iris[,3]),length=100)
x2 <- seq(min(iris[,4]),max(iris[,4]),length=100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))

gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(iris)[3:4]
probPreds <- predict(GPfitIris, gridPoints, type="probabilities")

# Plotting for Prob(setosa)
contour(x1,x2,matrix(probPreds[,1],100,byrow = TRUE), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Setosa) - Setosa is red')
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green")

# Plotting for Prob(Versicolor)
contour(x1,x2,matrix(probPreds[,2],100,byrow = TRUE), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Versicolor) - Versicolor is green')
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green")


# Plotting for Prob(virginica)
contour(x1,x2,matrix(probPreds[,3],100,byrow = TRUE), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Virginica) - Virginica is blue')
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green")

# Plotting the decision boundaries
meanPred <- matrix(max.col(probPreds),100,byrow = TRUE)
plot(gridPoints,  pch=".", cex=3, col=ifelse(meanPred==1, "red", ifelse(meanPred==2, "green", "blue")))
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red", cex=10, pch=".")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue",  cex=10, pch=".")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green",  cex=10, pch=".")
