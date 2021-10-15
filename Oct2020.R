##############################################################################
# Questions 1: BNs.

library(bnlearn)
library(gRain)

trials <- 1000
results <- matrix(data = NA, nrow = trials, ncol = 4)

for(i in 1:trials){
  net <- model2network("[D|C][C][A|C][Y|A:C]")
  cptC <- runif(2)
  dim(cptC) <- c(2)
  dimnames(cptC) <- list(c("1", "0"))
  cptC <- cptC/sum(cptC)
  cptD <- runif(4) 
  dim(cptD) <- c(2,2)
  dimnames(cptD) <- list("D" = c("1", "0"), "C" =  c("1", "0"))
  cptD <- prop.table(cptD,2)
  cptA <- runif(4) 
  dim(cptA) <- c(2,2)
  dimnames(cptA) <- list("A" = c("1", "0"), "C" =  c("1", "0"))
  cptA <- prop.table(cptA,2)
  cptY <- runif(8) 
  dim(cptY) <- c(2,2,2)
  dimnames(cptY) <- list("Y" = c("1", "0"), "A" =  c("1", "0"), "C" =  c("1", "0"))
  cptY <- prop.table(cptY,2:3)
  netfit <- custom.fit(net,list(C=cptC, D=cptD, A=cptA, Y=cptY))
  netcom <- compile(as.grain(netfit))
  
  pYAC <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("1","1")),c("Y"))
  pYAc <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("1","0")),c("Y"))
  pYaC <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("0","1")),c("Y"))
  pYac <- querygrain(setEvidence(netcom,nodes=c("A","C"),states=c("0","0")),c("Y"))

  nondecC <- (pYAc$Y[1] <= pYAC$Y[1] & pYac$Y[1] <= pYaC$Y[1])
  nonincC <- (pYAc$Y[1] >= pYAC$Y[1] & pYac$Y[1] >= pYaC$Y[1])
  
  pYAD <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("1","1")),c("Y"))
  pYAd <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("1","0")),c("Y"))
  pYaD <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("0","1")),c("Y"))
  pYad <- querygrain(setEvidence(netcom,nodes=c("A","D"),states=c("0","0")),c("Y"))

  nondecD <- (pYAd$Y[1] <= pYAD$Y[1] & pYad$Y[1] <= pYaD$Y[1])
  nonincD <- (pYAd$Y[1] >= pYAD$Y[1] & pYad$Y[1] >= pYaD$Y[1])
  
  results[i,] <- c(nondecC,nonincC,nondecD,nonincD)
}

colSums(results[which(results[,1]==FALSE & results[,2]==FALSE),])
colSums(results[which(results[,1]==FALSE & results[,2]==TRUE),])
colSums(results[which(results[,1]==TRUE & results[,2]==FALSE),])

##############################################################################
# Question 2: RL.

# install.packages("ggplot2")
# install.packages("vctrs")
library(ggplot2)

# If you do not see four arrows in line 16, then do the following:
# File/Reopen with Encoding/UTF-8

arrows <- c("↑", "→", "↓", "←")
action_deltas <- list(c(1,0), # up
                      c(0,1), # right
                      c(-1,0), # down
                      c(0,-1)) # left

vis_environment <- function(iterations=0, epsilon = 0.5, alpha = 0.1, gamma = 0.95, beta = 0){
  
  # Visualize an environment with rewards. 
  # Q-values for all actions are displayed on the edges of each tile.
  # The (greedy) policy for each state is also displayed.
  # 
  # Args:
  #   iterations, epsilon, alpha, gamma, beta (optional): for the figure title.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  #   H, W (global variables): environment dimensions.
  
  df <- expand.grid(x=1:H,y=1:W)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,1],NA),df$x,df$y)
  df$val1 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,2],NA),df$x,df$y)
  df$val2 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,3],NA),df$x,df$y)
  df$val3 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,q_table[x,y,4],NA),df$x,df$y)
  df$val4 <- as.vector(round(foo, 2))
  foo <- mapply(function(x,y) 
    ifelse(reward_map[x,y] == 0,arrows[GreedyPolicy(x,y)],reward_map[x,y]),df$x,df$y)
  df$val5 <- as.vector(foo)
  foo <- mapply(function(x,y) ifelse(reward_map[x,y] == 0,max(q_table[x,y,]),
                                     ifelse(reward_map[x,y]<0,NA,reward_map[x,y])),df$x,df$y)
  df$val6 <- as.vector(foo)
  
  print(ggplot(df,aes(x = y,y = x)) +
          scale_fill_gradient(low = "white", high = "green", na.value = "red", name = "") +
          geom_tile(aes(fill=val6)) +
          geom_text(aes(label = val1),size = 4,nudge_y = .35,na.rm = TRUE) +
          geom_text(aes(label = val2),size = 4,nudge_x = .35,na.rm = TRUE) +
          geom_text(aes(label = val3),size = 4,nudge_y = -.35,na.rm = TRUE) +
          geom_text(aes(label = val4),size = 4,nudge_x = -.35,na.rm = TRUE) +
          geom_text(aes(label = val5),size = 10) +
          geom_tile(fill = 'transparent', colour = 'black') + 
          ggtitle(paste("Q-table after ",iterations," iterations\n",
                        "(epsilon = ",epsilon,", alpha = ",alpha,"gamma = ",gamma,", beta = ",beta,")")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(breaks = c(1:W),labels = c(1:W)) +
          scale_y_continuous(breaks = c(1:H),labels = c(1:H)))
  
}

GreedyPolicy <- function(x, y){
  
  # Get a greedy action for state (x,y) from q_table.
  #
  # Args:
  #   x, y: state coordinates.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  # 
  # Returns:
  #   An action, i.e. integer in {1,2,3,4}.
  
  # Your code here.
  
  foo <- which(q_table[x,y,] == max(q_table[x,y,]))
  return (ifelse(length(foo)>1,sample(foo, size = 1),foo))
  
}

EpsilonGreedyPolicy <- function(x, y, epsilon){
  
  # Get an epsilon-greedy action for state (x,y) from q_table.
  #
  # Args:
  #   x, y: state coordinates.
  #   epsilon: probability of acting randomly.
  # 
  # Returns:
  #   An action, i.e. integer in {1,2,3,4}.
  
  # Your code here.
  
  foo <- sample(0:1,size = 1,prob = c(epsilon,1-epsilon))
  return (ifelse(foo == 1,GreedyPolicy(x,y),sample(1:4,size = 1)))
  
}

transition_model <- function(x, y, action, beta){
  
  # Computes the new state after given action is taken. The agent will follow the action 
  # with probability (1-beta) and slip to the right or left with probability beta/2 each.
  # 
  # Args:
  #   x, y: state coordinates.
  #   action: which action the agent takes (in {1,2,3,4}).
  #   beta: probability of the agent slipping to the side when trying to move.
  #   H, W (global variables): environment dimensions.
  # 
  # Returns:
  #   The new state after the action has been taken.
  
  delta <- sample(-1:1, size = 1, prob = c(0.5*beta,1-beta,0.5*beta))
  final_action <- ((action + delta + 3) %% 4) + 1
  foo <- c(x,y) + unlist(action_deltas[final_action])
  foo <- pmax(c(1,1),pmin(foo,c(H,W)))
  
  return (foo)
}

q_learning <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                       beta = 0, test = 0){
  
  # Just setting epsilon=0 instead of using the test argument is also OK.
  # But in that case the agent acts greedily while still updating the q-table.
  
  cur_pos <- start_state
  episode_correction <- 0
  ite <- 0
  repeat{
    # Follow policy, execute action, get reward.
    action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon*(1-test))
    new_pos <- transition_model(cur_pos[1], cur_pos[2], action, beta)
    reward <- reward_map[new_pos[1], new_pos[2]]
    
    # Q-table update.
    old_q <- q_table[cur_pos[1], cur_pos[2], action]
    correction <- ifelse(reward==0,-1,reward) + gamma*max(q_table[new_pos[1], new_pos[2], ]) - old_q
    q_table[cur_pos[1], cur_pos[2], action] <<- old_q + alpha*correction*(1-test)
    
    cur_pos <- new_pos
    episode_correction <- episode_correction + correction*(1-test)
    
    if(reward!=0)
      # End episode.
      return (c(reward-ite,episode_correction))
    else
      ite <- ite+1
  }
  
}

SARSA <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                  beta = 0, test = 0){
  
  cur_pos <- start_state
  cur_action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon*(1-test))
  episode_correction <- 0
  ite <- 0
  repeat{
    # Follow policy, execute action, get reward.
    new_pos <- transition_model(cur_pos[1], cur_pos[2], cur_action, beta)
    reward <- reward_map[new_pos[1], new_pos[2]]
    new_action <- EpsilonGreedyPolicy(new_pos[1], new_pos[2], epsilon*(1-test))
    
    # Q-table update.
    old_q <- q_table[cur_pos[1], cur_pos[2], cur_action]
    correction <- ifelse(reward==0,-1,reward) + gamma*q_table[new_pos[1], new_pos[2], new_action] - old_q
    q_table[cur_pos[1], cur_pos[2], cur_action] <<- old_q + alpha*correction*(1-test)
    
    cur_pos <- new_pos
    cur_action <- new_action
    episode_correction <- episode_correction + correction*(1-test)
    
    if(reward!=0)
      # End episode.
      return (c(reward-ite,episode_correction))
    else
      ite <- ite+1
  }
  
}

MovingAverage <- function(x, n){
  
  cx <- c(0,cumsum(x))
  rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  
  return (rsum)
}

H <- 3
W <- 6

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[1,2:5] <- -10
reward_map[1,6] <- 10
# To avoid having to modify vis_environment, I take care of the reward -1 in the functions
# q_learning and SARSA.

q_table <- array(0,dim = c(H,W,4))

rewardqtr <- NULL
for(i in 1:5000){
  foo <- q_learning(start_state = c(1,1), epsilon = 0.5, gamma = 1, beta = 0, test = 0)
  rewardqtr <- c(rewardqtr,foo[1])
}

vis_environment(5000, epsilon = 0.5, gamma = 1, beta = 0)

rewardqte <- NULL
for(i in 1:5000){
  foo <- q_learning(start_state = c(1,1), epsilon = 0.5, gamma = 1, beta = 0, test = 1)
  rewardqte <- c(rewardqte,foo[1])
}

q_table <- array(0,dim = c(H,W,4))

rewardstr <- NULL
for(i in 1:5000){
  foo <- SARSA(start_state = c(1,1), epsilon = 0.5, gamma = 1, beta = 0, test = 0)
  rewardstr <- c(rewardstr,foo[1])
}

vis_environment(5000, epsilon = 0.5, gamma = 1, beta = 0)

rewardste <- NULL
for(i in 1:5000){
  foo <- SARSA(start_state = c(1,1), epsilon = 0.5, gamma = 1, beta = 0, test = 1)
  rewardste <- c(rewardste,foo[1])
}

plot(MovingAverage(rewardqtr,100),type = "l",ylim = c(-15,5))
lines(MovingAverage(rewardqte,100),type = "l",lty=2)
lines(MovingAverage(rewardstr,100),type = "l",col = "blue")
lines(MovingAverage(rewardste,100),type = "l",col = "blue",lty=2)

# During training Q-learning performs worse because it takes the shortest route to the 
# goal state, which means the agent falling off the cliff now and then due to epsilon. 
# Q-learning prefers the shortest path because it assumes in the updating rule that the 
# subsequent moves will be greedy and, thus, the agent will never fall off the cliff. 
# The reality is that the moves are not greedy due to epsilon. During testing Q-learning
# performs better because epsilon is zero and thus the agent never falls off the cliff.

##############################################################################
# Question 3: GPs.

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

sigmaF <- 1
l <- 0.3
xData <- c(-1,-0.6,-0.2,0.4,0.8)
yData <- c(0.768,-0.044,-0.94,0.719,-0.664)
xGrid <- seq(-1,1,0.01)
res<-posteriorGP(X=xData,y=yData,k=SEKernel,sigmaNoise=0,xStar=xGrid)
plot(xData,yData,xlim=c(-1,1),ylim=c(-0.5,0.5))
xGrid[101]
lines(xGrid, res$var[101,], col = "green")
abline(h=0)
abline(v=-1)
abline(v=-0.6)
abline(v=-0.2)
abline(v=0.4)
abline(v=0.8)
abline(v=0)

foo <- SEKernel(xGrid,0)
plot(xGrid,foo)

# The posterior covariance is zero at the training points because the functions must
# go through them, i.e. there is no uncertainty due to this being a noisy-free problem.
# The posterior covariance is not monotone decreasing with the distance because it is
# constrained by the fact of being zero in the training points.

##############################################################################
# Hyperparameter search for sigmaNoise.

tempData <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv', header=TRUE, sep=';')
temp <- tempData$temp
plot(temp, type="l")
time = 1:length(temp)
day = rep(1:365,6)

# Extract every 5:th observation
subset <- seq(1, length(temp), by = 5)
temp <- temp[subset]
time = time[subset]
plot(time,temp, type="l")

sigmaF <- 20
l <- 0.2
polyFit <- lm(scale(temp) ~  scale(time) + I(scale(time)^2))
sigmaNoiseFit = sd(polyFit$residuals)
res<-posteriorGP(X=scale(time),y=scale(temp),k=SEKernel,sigmaNoise=sigmaNoiseFit,xStar=scale(time))
lines(time, res$mu*sd(temp)+mean(temp), col="green", lwd = 2)
lines(time, res$mu*sd(temp)+mean(temp) - 1.96*sd(temp)*sqrt(diag(res$var)), col = "red")
lines(time, res$mu*sd(temp)+mean(temp) + 1.96*sd(temp)*sqrt(diag(res$var)), col = "red")

LM <- function(X,y,k,par){
  n <- length(y)
  L <- t(chol(k(X,X)+((par^2)*diag(n))))
  a <- solve(t(L),solve(L,y))
  logmar <- -0.5*(t(y)%*%a)-sum(log(diag(L)))-(n/2)*log(2*pi)
  return(logmar)
}

# Grid search.

besti<- 0.1
bestLM<-LM(X=scale(time),y=scale(temp),k=SEKernel,par=besti)
bestLM
besti
for(i in seq(0.2,10,0.1)){
    aux<-LM(X=scale(time),y=scale(temp),k=SEKernel,par=i)
    if(bestLM<aux){
      bestLM<-aux
      besti<-i
    }
  }
bestLM
besti

res<-posteriorGP(X=scale(time),y=scale(temp),k=SEKernel,sigmaNoise=besti,xStar=scale(time))
lines(time, res$mu*sd(temp)+mean(temp), col="green", lwd = 2)
lines(time, res$mu*sd(temp)+mean(temp) - 1.96*sd(temp)*sqrt(diag(res$var)), col = "blue")
lines(time, res$mu*sd(temp)+mean(temp) + 1.96*sd(temp)*sqrt(diag(res$var)), col = "blue")

# optim.

foo<-optim(par = 0.1, fn = LM, X=scale(time),y=scale(temp),k=SEKernel, method="L-BFGS-B",
           lower = c(.Machine$double.eps),control=list(fnscale=-1))
foo$value
foo$par

res<-posteriorGP(X=scale(time),y=scale(temp),k=SEKernel,sigmaNoise=besti,xStar=scale(time))
lines(time, res$mu*sd(temp)+mean(temp), col="green", lwd = 2)
lines(time, res$mu*sd(temp)+mean(temp) - 1.96*sd(temp)*sqrt(diag(res$var)), col = "blue")
lines(time, res$mu*sd(temp)+mean(temp) + 1.96*sd(temp)*sqrt(diag(res$var)), col = "blue")
