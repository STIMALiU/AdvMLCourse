# Exercise 1.

# Independencies.

# A _|_ C | B,D
# B _|_ D | A

# Skeleton.

# {A - B, B - C, A - D, D - C}

# Unshielded colliders.

# {B -> C <- D}

# Markov equivalent DAGs.

# {A -> B, B -> C, A <- D, D -> C}
# {A <- B, B -> C, A -> D, D -> C}

# Essential graph.

# {A - B, B -> C, A - D, D -> C}

# Exercise 2.

# HMM library.

library(HMM)

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
                       .2,.2,0,0,0,0,0,.2,.2,.2), nrow=length(States), ncol=length(Symbols), byrow = TRUE)

startProbs=c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)

hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)

viterbi(hmm,c(3,8,9))

# Own implementation.

x <- c(3,8,9)

w <- matrix(NA,length(States),length(x)) 
v <- matrix(NA,length(States),length(x))

w[,1] <- log(startProbs) + log(emissionProbs[x[1],])

for (t in 1:(length(x)-1)) {
  for (i in 1:length(States)) {
    w[i,t+1] <- log(emissionProbs[x[t+1],i]) + max(log(transProbs[,i]) + w[,t])
    v[i,t+1] <- which.max(log(transProbs[,i]) + w[,t])
  }
}

z <- vector("integer",length(x))
z[length(x)] <- which.max(w[,length(x)])

for (t in (length(x)-1):1) {
  z[t] <- v[z[t+1],t+1]
}

z

# Exercise 3.

library(ggplot2)

arrows <- c("^", ">", "v", "<")
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
                       beta = 0){
  
  # Perform one episode of Q-learning. The agent should move around in the 
  # environment using the given transition model and update the Q-table.
  # The episode ends when the agent reaches a terminal state.
  # 
  # Args:
  #   start_state: array with two entries, describing the starting position of the agent.
  #   epsilon (optional): probability of acting randomly.
  #   alpha (optional): learning rate.
  #   gamma (optional): discount factor.
  #   beta (optional): slipping factor.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  # 
  # Returns:
  #   reward: reward received in the episode.
  #   correction: sum of the temporal difference correction terms over the episode.
  #   q_table (global variable): Recall that R passes arguments by value. So, q_table being
  #   a global variable can be modified with the superassigment operator <<-.
  
  cur_pos <- start_state
  episode_correction <- 0
  
  repeat{
    # Follow policy, execute action, get reward.
    action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon)
    new_pos <- transition_model(cur_pos[1], cur_pos[2], action, beta)
    reward <- reward_map[new_pos[1], new_pos[2]]
    
    # Q-table update.
    old_q <- q_table[cur_pos[1], cur_pos[2], action]
    correction <- reward + gamma*max(q_table[new_pos[1], new_pos[2], ]) - old_q
    q_table[cur_pos[1], cur_pos[2], action] <<- old_q + alpha*correction
    
    cur_pos <- new_pos
    episode_correction <- episode_correction + correction
    
    if(reward!=0)
      # End episode.
      return (c(reward,episode_correction))
  }
  
}

expected_SARSA <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, 
                                    beta = 0){
  
  # Perform one episode of Q-learning. The agent should move around in the 
  # environment using the given transition model and update the Q-table.
  # The episode ends when the agent reaches a terminal state.
  # 
  # Args:
  #   start_state: array with two entries, describing the starting position of the agent.
  #   epsilon (optional): probability of acting randomly.
  #   alpha (optional): learning rate.
  #   gamma (optional): discount factor.
  #   beta (optional): slipping factor.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  #   q_table (global variable): a HxWx4 array containing Q-values for each state-action pair.
  # 
  # Returns:
  #   reward: reward received in the episode.
  #   correction: sum of the temporal difference correction terms over the episode.
  #   q_table (global variable): Recall that R passes arguments by value. So, q_table being
  #   a global variable can be modified with the superassigment operator <<-.
  
  cur_pos <- start_state
  episode_correction <- 0
  
  repeat{
    # Follow policy, execute action, get reward.
    action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon)
    new_pos <- transition_model(cur_pos[1], cur_pos[2], action, beta)
    reward <- reward_map[new_pos[1], new_pos[2]]
    
    # Q-table update.
    old_q <- q_table[cur_pos[1], cur_pos[2], action]
    exp <- max(q_table[new_pos[1], new_pos[2], ])*(1-epsilon) + sum(q_table[new_pos[1], new_pos[2], ])*epsilon/4
    correction <- reward + gamma*exp - old_q
    q_table[cur_pos[1], cur_pos[2], action] <<- old_q + alpha*correction
    
    cur_pos <- new_pos
    episode_correction <- episode_correction + correction
    
    if(reward!=0)
      # End episode.
      return (c(reward,episode_correction))
  }
  
}

# Environment A.

H <- 5
W <- 7

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[3,6] <- 10
reward_map[2:4,3] <- -1

MovingAverage <- function(x, n){
  
  cx <- c(0,cumsum(x))
  rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  
  return (rsum)
}

q_table <- array(0,dim = c(H,W,4))

reward <- NULL

for(i in 1:10000){
  foo <- q_learning(start_state = c(3,1))
  reward <- c(reward,foo[1])
}

plot(MovingAverage(reward,100),type = "l")

q_table <- array(0,dim = c(H,W,4))

reward <- NULL

for(i in 1:10000){
  foo <- expected_SARSA(start_state = c(3,1))
  reward <- c(reward,foo[1])
}

plot(MovingAverage(reward,100),type = "l")

# Expected SARSA outperforms Q-learning because it uses the true next action (in expectation) in the 
# updating rule, as opposed to Q-learning which always uses the best action in the updating rule, 
# although that action is not always chosen as the next action due to epsilon.

# Exercise 4.

library(kernlab)

tempData <- read.csv('https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv', header=TRUE, sep=';')
temp <- tempData$temp
time = 1:length(temp)

# Extract every 5:th observation
subset <- seq(1, length(temp), by = 5)
temp <- temp[subset]
time = time[subset]
plot(time,temp, type="l")

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

sigmaF <- 20
l <- 100
polyFit <- lm(temp ~  time + I(time^2))
sigmaNoiseFit = sd(polyFit$residuals)
res<-posteriorGP(X=time,y=temp,k=SEKernel,sigmaNoise=sigmaNoiseFit,xStar=time)
plot(time,temp, type="l")
lines(time, res$mu, col="blue", lwd = 4)
lines(time, res$mu - 1.96*sqrt(diag(res$var)), col = "red", lwd = 4)
lines(time, res$mu + 1.96*sqrt(diag(res$var)), col = "red", lwd = 4)

# Nyström approximation.

x <- time
xs <- time
n <- length(x)
z <- time[sample(1:length(time),100)]

sigmaF <- 20
l <- 100
polyFit <- lm(temp[z] ~  z + I(z^2))
sigmaNoiseFit = sd(polyFit$residuals)

Kxs <- SEKernel(x,xs)
Kss <- SEKernel(xs,xs)
Kxz <- SEKernel(x,z)
Kzz <- SEKernel(z,z)

Nystrom <- diag(n)/sigmaNoiseFit^2 - Kxz%*%solve(Kzz+t(Kxz)%*%Kxz/sigmaNoiseFit^2,t(Kxz)/sigmaNoiseFit^4, tol = 0)

Meanf <- t(Kxs)%*%Nystrom%*%temp
Covf = Kss-t(Kxs)%*%Nystrom%*%Kxs
lines(time, Meanf, col="black", lwd = 2)
lines(time, Meanf - 1.96*sqrt(diag(Covf)), col = "black", lwd = 2)
lines(time, Meanf + 1.96*sqrt(diag(Covf)), col = "black", lwd = 2)

# The Nyström approximation is suitable for large datasets, since it only requires the inversion of a
# matrix of dimension m times m where m is the size of z (i.e. the number of inducing points) which is
# typically much smaller than the total number of points in the dataset.

