# BNs.

library(bnlearn)
library(gRain)

data("asia")
net<-model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
net<-bn.fit(net,asia,method="bayes")

mydata<-matrix(NA,1000,8)
for(i in 1:1000){
  a<-sample(1:2,1,prob=net$A$prob)
  s<-sample(1:2,1,prob=net$S$prob)
  t<-sample(1:2,1,prob=net$T$prob[,a])
  
  l<-sample(1:2,1,prob=net$L$prob[,s])
  b<-sample(1:2,1,prob=net$B$prob[,s])
  
  e<-sample(1:2,1,prob=net$E$prob[,l,t])
  d<-sample(1:2,1,prob=net$D$prob[,b,e])
  
  x<-sample(1:2,1,prob=net$X$prob[,e])
  
  mydata[i,]<-c(a,s,t,l,b,d,e,x)
}

foo<-mydata[which(mydata[,6]==2),2]
table(foo)/length(foo)

net<-as.grain(net)
net<-compile(net)
net<-setEvidence(net,nodes=c("D"),states=c("yes"))
querygrain(net,c("S"))

# HMMs.

library(HMM)
library(entropy)
set.seed(567)
States=1:3 # healthy, sick1, sick2
Symbols=1:2 #healthy, sick
transProbs=matrix(c(.9,.1,0,
                    0,0,1,
                    .2,0,.8), nrow=length(States), ncol=length(States), byrow = TRUE)
emissionProbs=matrix(c(.6,.4,
                       .3,.7,
                       .3,.7), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
startProbs=c(.5,.5,0)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)
simHMM(hmm,100)

# RL.

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
                       beta = 0, tr = 1){
  
  # Perform one episode of Q-learning. The agent should move around in the 
  # environment using the given transition model and update the Q-table.
  # The episode ends when the agent reaches a terminal state.
  # 
  # Args:
  #   start_state: array with two entries, describing the starting position of the agent.
  #   epsilon (optional): probability of acting greedily.
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
    action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon*tr)
    new_pos <- transition_model(cur_pos[1], cur_pos[2], action, beta)
    reward <- reward_map[new_pos[1], new_pos[2]]
    
    # Q-table update.
    old_q <- q_table[cur_pos[1], cur_pos[2], action]
    correction <- reward + gamma*max(q_table[new_pos[1], new_pos[2], ]) - old_q
    q_table[cur_pos[1], cur_pos[2], action] <<- old_q + alpha*correction*tr
    
    cur_pos <- new_pos
    episode_correction <- episode_correction + correction*tr
    
    if(reward!=0)
      # End episode.
      return (c(reward,episode_correction))
  }
  
}

#####################################################################################################
# Q-Learning Environments
#####################################################################################################

# Environment B (the effect of epsilon and gamma)

H <- 7
W <- 8

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[1,] <- -1
reward_map[7,] <- -1
reward_map[4,5] <- 5
reward_map[4,8] <- 10

q_table <- array(0,dim = c(H,W,4))

vis_environment()

mreward <- NULL

for(i in c(0.1,0.25,0.5)){
  for(j in c(0.5,0.75,0.95)){
    q_table <- array(0,dim = c(H,W,4))
    reward <- NULL
    
    for(k in 1:30000){
      q_learning(epsilon = i, gamma = j, start_state = c(4,1), tr = 1) # training
    }
    
    for(k in 1:1000){
      foo <- q_learning(epsilon = i, gamma = j, start_state = c(4,1), tr = 0) # validation
      reward <- c(reward,foo[1])
    }
    
    vis_environment(k, epsilon = i, gamma = j)
    mreward <- c(mreward,mean(reward))
  }
}

mreward

# GPs.

# ell=1 is too smooth, and sigmaNoise=1 implies that the points are not trustworthy as
# function values. So, the GP does not really try to go through them and the probability
# bands are wide. Setting sigmaNoise=.1 reduces the probability bands drastically, as the
# points are now trustworthy. However, the GP is still too smooth to fit the data well.

# ell=0.3 is quite flexible and, thus, the GP goes through the points with sigmaNoise=0.1.
# However, with sigmaNoise=1, the GP does not see any reason to go through the points, since
# they are very noisy versions of the function values. Hence the probability bands are wide too.

# See also p. 21 in the book by Rasmussen and Williams.