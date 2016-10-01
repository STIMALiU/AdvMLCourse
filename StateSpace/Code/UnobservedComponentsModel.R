# Unobserved components model
# My own spin of the analysis in the vignette to the rucm package in R
# Mattias Villani, http://mattiasvillani.com


#install.packages('rucm')
library(rucm)

# Plotting the data
# The Nile dataset in is the datasets package
timeNile <- seq(1871, 1970)
nTraining <- length(timeNile) # Total number of data points
nPreds <- 20
par(cex.axis = 0.8, cex.main = 0.8, cex.lab = 0.8)
plot(Nile, ylab = "Flow of Nile", lwd = 2, xlab = 'year', 
     xlim = c(1,nTraining+nPreds), main = "Modeling the Nile flow")
lines(seq(1,nTraining), Nile, lwd = 2)

#### Fitting a local level model ###
# y(t)  = mu(t) + eps(t)   [sigma2_eps is the Irregular_Variance in rucm]
# mu(t) = mu(t-1) + nu(t)  [sigma2_nu is the Level_Variance in rucm]
localLevelModel <- ucm(formula = Nile~0, data = Nile, level = TRUE)
localLevelModel
lines(seq(1,nTraining), localLevelModel$s.level, col = "red", lwd = 2)

#### Fitting a local trend model ###
# y(t)  = mu(t) + eps(t)           [sigma2_eps is the Irregular_Variance in rucm]
# mu(t) = mu(t-1) + beta_t + nu(t) [sigma2_nu is the Level_Variance in rucm]
# beta(t) = beta(t-1) + eta(t)     [sigma2_eta is the Slope_Variance in rucm]
localTrendModel <- ucm(formula = Nile~0, data = Nile, level = TRUE, slope = TRUE)

# Adding the fit from the local trend model
#plot(Nile, ylab = "Flow of Nile", lwd = 2, xlab = 'year')
lines(seq(1,nTraining), localTrendModel$s.level + localTrendModel$s.slope, col = "blue", lwd = 2, xlab = 'year')

# Make predictions from local level model
predsLocalLevel <- predict(localLevelModel$model, n.ahead = nPreds)
predsLocalTrend <- predict(localTrendModel$model, n.ahead = nPreds)

lines(seq(nTraining + 1,nTraining + nPreds), predsLocalLevel, 
      type = "l", lwd = 2, col = "red", lty = 2)
lines(seq(nTraining + 1,nTraining + nPreds), predsLocalTrend, 
      type = "l", lwd = 2, col = "blue", lty = 2)
legend("topright", legend = c("Data","Fit local level", "Fit local trend",
                              "Predictions local level","Predictions local trend"), 
       col = c("black","red","blue","red","blue"), lty = c(1,1,1,2,2), lwd = 2, cex = 0.45)


# Plotting the level and trend components separately. 
plot(timeNile, localTrendModel$s.level, type = "l", col = "red", lwd = 2, 
     main = "Inferred level", xlab = 'year', ylab = "")
plot(timeNile, localTrendModel$s.slope, type = "l", col = "red", lwd = 2, 
     main = "Inferred slope", xlab = 'year', ylab = "")


