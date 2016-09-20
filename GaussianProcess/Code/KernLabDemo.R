# Kernlab demo for Gaussian processes regression and classfication

setwd('/Users/mv/Dropbox/Teaching/AdvML/GaussianProcess/Code')
#install.packages('kernlab')
library(kernlab)
#install.packages("AtmRay") # To make 2D grid like in Matlab's meshgrid
library(AtmRay)

### Regression on the LIDAR data  ###
lidarData <- read.table('LidarData', header = T)
LogRatio <- lidarData$LogRatio
Distance <- lidarData$Distance

plot(Distance,LogRatio)
GPfit <- gausspr(Distance, LogRatio)
ytest <- predict(GPfit, Distance)
lines(Distance,ytest,col="red")


GPfit <- gausspr(Distance, LogRatio, kernel = 'rbfdot', par = list(sigma = 1), variance.model = TRUE)

plot(Distance,LogRatio)
meanPred <- predict(GPfit, Distance)
stdPred <- predict(GPfit, Distance, type="sdeviation")
lines(Distance, meanPred)
lines(Distance, meanPred + 1.96*stdPred, col="red")
lines(Distance, meanPred - 1.96*stdPred, col="red")




#### Regression on the Canadian wages data ####
wagesData <- read.table('CanadianWages.dat', header = T)
logWage <- wagesData$logWage
age <- wagesData$age
plot(age, logWage)
GPfit <- gausspr(logWage ~ age, kernel = 'rbfdot', par = list(sigma = 1))
meanPred <- predict(GPfit, age)
lines(age, meanPred, col="red")


### Regression on the Japan temperature data ####
tempData <- read.table('JapanTemp.dat', header = TRUE)
temp <- tempData$temp
time <- tempData$time
plot(time, temp)
GPfit <- gausspr(temp ~ time, kernel = 'rbfdot', par = list(sigma = 1))
meanPred <- predict(GPfit, time)
lines(time, meanPred, col="red")



#### Classification on Iris data ###
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
contour(x1,x2,matrix(probPreds[,1],100), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Setosa) - Setosa is red')
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green")

# Plotting for Prob(Versicolor)
contour(x1,x2,matrix(probPreds[,2],100), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Versicolor) - Versicolor is green')
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green")


# Plotting for Prob(virginica)
contour(x1,x2,matrix(probPreds[,3],100), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Virginica) - Virginica is blue')
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green")

# Plotting the decision boundaries
meanPred <- matrix(max.col(probPreds),100)
plot(gridPoints,  pch=".", cex=3, col=ifelse(meanPred==1, "red", ifelse(meanPred==2, "green", "blue")))
points(iris[iris[,5]=='setosa',3],iris[iris[,5]=='setosa',4],col="red", cex=10, pch=".")
points(iris[iris[,5]=='virginica',3],iris[iris[,5]=='virginica',4],col="blue",  cex=10, pch=".")
points(iris[iris[,5]=='versicolor',3],iris[iris[,5]=='versicolor',4],col="green",  cex=10, pch=".")
