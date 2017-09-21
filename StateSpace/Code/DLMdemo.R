# Demo of the dlm package using the nile data
# Author: Mattias Villani, http://mattiasvillani.com

#install.packages("dlm")
library(dlm)

# Plotting the data
# The Nile dataset in is the datasets package
timeNile <- seq(1871, 1970)
nTraining <- length(timeNile) # Total number of data points
nPreds <- 20
par(cex.axis = 0.8, cex.main = 0.8, cex.lab = 0.8)
plot(Nile, ylab = "Flow of Nile", lwd = 2, xlab = 'year', 
     xlim = c(1,nTraining+nPreds), main = "Modeling the Nile flow")
lines(seq(1,nTraining), Nile, lwd = 2)

# Set up the local trend model
# the state is theta_t = (mu_t, beta_)'
m0 = c(0,0)     # Prior mean for the first state theta_0
C0 = 100*diag(2) # Prior covariance for the first state theta_0
FF = matrix(c(1,0),1,2)     # The mapping from the states to the measurements. 
                            # Note FF must be a matrix object even if it is a vector ...
GG = matrix(c(1,1,0,1), 2, 2, byrow = T)
V = 10
W = 0.01*diag(2)

dlmLocalTrend <- dlm(m0 = m0,C0 = C0, FF = FF, GG = GG, V = V, W = W)

# Kalman filtering
dlmKFfilt <- dlmFilter(y = Nile, dlmLocalTrend)

# Plotting the sum of the filtered states
lines(seq(1,nTraining), dlmKFfilt$m[2:(nTraining+1),1] + dlmKFfilt$m[2:(nTraining+1),2], col = "red", lwd = 2)

# Plotting the predicted response, y
lines(seq(1,nTraining), dlmKFfilt$f, col = "green", lwd = 2)

# State smoothing
dlmSmoothed <- dlmSmooth(y = Nile, dlmLocalTrend)

# Plotting the data the sum of the filtered and smoothed series for the two states.
par(cex.axis = 0.8, cex.main = 0.8, cex.lab = 0.8)
plot(Nile, ylab = "Flow of Nile", lwd = 2, xlab = 'year', 
     xlim = c(1,nTraining+nPreds), main = "Modeling the Nile flow")
lines(seq(1,nTraining), Nile, lwd = 2)

# Plotting the sum of the filtered states
lines(seq(1,nTraining), dlmKFfilt$m[2:(nTraining+1),1] + dlmKFfilt$m[2:(nTraining+1),2], col = "red", lwd = 2)

# Plotting the sum of the smoothed states. Note that dlmSmoothed$s[1,] is theta at time t=0
lines(seq(1,nTraining), dlmSmoothed$s[2:(nTraining+1),1] + dlmSmoothed$s[2:(nTraining+1),2], col = "blue", lwd = 2)


# Parameter estimation

# First create so called build function that takes the parameters that you want to estimate
# as input and returns a list (or a dlm object) with the state space matrices as output
buildLocalTrend <- function(x){
  
  V = exp(x[1])   # NOTE THE EXPONENTIAL. x must be unrestricted parameters, otherwise optim will go bananas
  W = diag(exp(x[2:3]),2,2)
  
  return(dlm(
    m0 = rep(0,2),
    C0 = 100*diag(2),   
    FF = matrix(c(1,0),1,2),     
    GG = matrix(c(1,1,0,1), 2, 2, byrow = T),
    V = V,
    W = W))
  
}


initVal <- c(1,1,1) # Initial values for optim on the estimated parameters 
dlmLocalTrend <- buildLocalTrend(initVal) # Just to see that the build function works
MLEs <- dlmMLE(Nile, parm = initVal, build = buildLocalTrend)
dlmWithMLEs <- buildLocalTrend(MLEs$par)

# Let's plot the filtered and smoothed states using the MLEs

# First plot the Nile data
par(cex.axis = 0.8, cex.main = 0.8, cex.lab = 0.8)
plot(Nile, ylab = "Flow of Nile", lwd = 2, xlab = 'year', 
     xlim = c(1,nTraining+nPreds), main = "Modeling the Nile flow")
lines(seq(1,nTraining), Nile, lwd = 2)

# Filtering
dlmKFfilt <- dlmFilter(y = Nile, dlmWithMLEs)

# Plotting the sum of the filtered states
lines(seq(1,nTraining), dlmKFfilt$m[2:(nTraining+1),1] + dlmKFfilt$m[2:(nTraining+1),2], col = "red", lwd = 2)

# State smoothing
dlmSmoothed <- dlmSmooth(y = Nile, dlmWithMLEs)

# Plotting the sum of the smoothed states. Note that dlmSmoothed$s[1,] is theta at time t=0
lines(seq(1,nTraining), dlmSmoothed$s[2:(nTraining+1),1] 
      + dlmSmoothed$s[2:(nTraining+1),2], col = "blue", lwd = 2)


# Bayesian analysis
nDraws <- 1000
nStates <- 2
postDraws <- array(NA, c(nTraining, nDraws, nStates) )
for (i in 1:nDraws){
  postDraws[,i,] <- dlmBSample(dlmKFfilt)[-1,]
}

# Plotting the smoothed series (posterior mean)
# and the first 3 draws from posterior
par(cex.axis = 0.8, cex.main = 0.8, cex.lab = 0.8)
plot(Nile, ylab = "Flow of Nile", lwd = 2, xlab = 'year', 
     xlim = c(1,nTraining+nPreds), main = "Modeling the Nile flow")
lines(seq(1,nTraining), Nile, lwd = 2)
lines(seq(1,nTraining), dlmSmoothed$s[2:(nTraining+1),1] 
      + dlmSmoothed$s[2:(nTraining+1),2], col = "red", lwd = 2)
for (i in 1:3){
  lines(seq(1,nTraining), postDraws[,i,1] + postDraws[,i,2], col = "blue", lwd = 1)
}

# Plotting the marginal posterior of the state at time t
hist(postDraws[50,,1], 40)

