##############################################################################
# Question 1: PGMs.

library(bnlearn)

data("lizards")

lizardsnet<-model2network("[Species][Diameter|Species][Height|Species]") # True DAG
plot(lizardsnet)
plot(cpdag(lizardsnet)) # Plot the true pattern

ci.test("Diameter", "Species", test = "x2", data = lizards) # Keep edge D-S.
ci.test("Height", "Species", test = "x2", data = lizards) # Keep edge H-S.
ci.test("Height", "Diameter", test = "x2", data = lizards) # Remove edge D-H.
ci.test("Diameter", "Species", "Height", test = "x2", data = lizards) # Keep edge D-S.
ci.test("Height", "Species", "Diameter", test = "x2", data = lizards) # Keep edge H-S.

# Orientate D->S<-H. Wrong model !

##############################################################################
# Question 2: HMMs.

library(bnlearn)
library(gRain)

net <- model2network("[Z0][X0|Z0][Z1|Z0][X1|Z1][Z2|Z1][X2|Z2][Z3|Z2]")
plot(net)

cptZ0 <- c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
dim(cptZ0) <- c(10)
dimnames(cptZ0) <- list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

cptZ1 <- matrix(c(.5,.5,0,0,0,0,0,0,0,0,
                  0,.5,.5,0,0,0,0,0,0,0,
                  0,0,.5,.5,0,0,0,0,0,0,
                  0,0,0,.5,.5,0,0,0,0,0,
                  0,0,0,0,.5,.5,0,0,0,0,
                  0,0,0,0,0,.5,.5,0,0,0,
                  0,0,0,0,0,0,.5,.5,0,0,
                  0,0,0,0,0,0,0,.5,.5,0,
                  0,0,0,0,0,0,0,0,.5,.5,
                  .5,0,0,0,0,0,0,0,0,.5), nrow=10, ncol=10)
dim(cptZ1) <- c(10,10)
dimnames(cptZ1) <- list("Z1" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), "Z0" =  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

cptZ2 <- matrix(c(.5,.5,0,0,0,0,0,0,0,0,
                  0,.5,.5,0,0,0,0,0,0,0,
                  0,0,.5,.5,0,0,0,0,0,0,
                  0,0,0,.5,.5,0,0,0,0,0,
                  0,0,0,0,.5,.5,0,0,0,0,
                  0,0,0,0,0,.5,.5,0,0,0,
                  0,0,0,0,0,0,.5,.5,0,0,
                  0,0,0,0,0,0,0,.5,.5,0,
                  0,0,0,0,0,0,0,0,.5,.5,
                  .5,0,0,0,0,0,0,0,0,.5), nrow=10, ncol=10)
dim(cptZ2) <- c(10,10)
dimnames(cptZ2) <- list("Z2" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), "Z1" =  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

cptZ3 <- matrix(c(.5,.5,0,0,0,0,0,0,0,0,
                  0,.5,.5,0,0,0,0,0,0,0,
                  0,0,.5,.5,0,0,0,0,0,0,
                  0,0,0,.5,.5,0,0,0,0,0,
                  0,0,0,0,.5,.5,0,0,0,0,
                  0,0,0,0,0,.5,.5,0,0,0,
                  0,0,0,0,0,0,.5,.5,0,0,
                  0,0,0,0,0,0,0,.5,.5,0,
                  0,0,0,0,0,0,0,0,.5,.5,
                  .5,0,0,0,0,0,0,0,0,.5), nrow=10, ncol=10)
dim(cptZ3) <- c(10,10)
dimnames(cptZ3) <- list("Z3" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), "Z2" =  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

cptX0 <- matrix(c(.2,.2,.2,0,0,0,0,0,.2,.2,
                  .2,.2,.2,.2,0,0,0,0,0,.2,
                  .2,.2,.2,.2,.2,0,0,0,0,0,
                  0,.2,.2,.2,.2,.2,0,0,0,0,
                  0,0,.2,.2,.2,.2,.2,0,0,0,
                  0,0,0,.2,.2,.2,.2,.2,0,0,
                  0,0,0,0,.2,.2,.2,.2,.2,0,
                  0,0,0,0,0,.2,.2,.2,.2,.2,
                  .2,0,0,0,0,0,.2,.2,.2,.2,
                  .2,.2,0,0,0,0,0,.2,.2,.2), nrow=10, ncol=10)
dim(cptX0) <- c(10,10)
dimnames(cptX0) <- list("X0" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), "Z0" =  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

cptX1 <- matrix(c(.2,.2,.2,0,0,0,0,0,.2,.2,
                  .2,.2,.2,.2,0,0,0,0,0,.2,
                  .2,.2,.2,.2,.2,0,0,0,0,0,
                  0,.2,.2,.2,.2,.2,0,0,0,0,
                  0,0,.2,.2,.2,.2,.2,0,0,0,
                  0,0,0,.2,.2,.2,.2,.2,0,0,
                  0,0,0,0,.2,.2,.2,.2,.2,0,
                  0,0,0,0,0,.2,.2,.2,.2,.2,
                  .2,0,0,0,0,0,.2,.2,.2,.2,
                  .2,.2,0,0,0,0,0,.2,.2,.2), nrow=10, ncol=10)
dim(cptX1) <- c(10,10)
dimnames(cptX1) <- list("X1" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), "Z1" =  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

cptX2 <- matrix(c(.2,.2,.2,0,0,0,0,0,.2,.2,
                  .2,.2,.2,.2,0,0,0,0,0,.2,
                  .2,.2,.2,.2,.2,0,0,0,0,0,
                  0,.2,.2,.2,.2,.2,0,0,0,0,
                  0,0,.2,.2,.2,.2,.2,0,0,0,
                  0,0,0,.2,.2,.2,.2,.2,0,0,
                  0,0,0,0,.2,.2,.2,.2,.2,0,
                  0,0,0,0,0,.2,.2,.2,.2,.2,
                  .2,0,0,0,0,0,.2,.2,.2,.2,
                  .2,.2,0,0,0,0,0,.2,.2,.2), nrow=10, ncol=10)
dim(cptX2) <- c(10,10)
dimnames(cptX2) <- list("X2" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), "Z2" =  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

netfit <- custom.fit(net,list(Z0=cptZ0, Z1=cptZ1, Z2=cptZ2, Z3=cptZ3, X0=cptX0, X1=cptX1, X2=cptX2))
netcom <- compile(as.grain(netfit))

querygrain(setEvidence(netcom,nodes=c("X0","X2"),states=c("1","3")),c("Z0"))
querygrain(setEvidence(netcom,nodes=c("X0","X2"),states=c("1","3")),c("Z1"))
querygrain(setEvidence(netcom,nodes=c("X0","X2"),states=c("1","3")),c("Z2"))
querygrain(setEvidence(netcom,nodes=c("X0","X2"),states=c("1","3")),c("Z3"))

##############################################################################
# Question 3: RL.

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

SARSA <- function(start_state, epsilon = 0.5, alpha = 0.1, gamma = 0.95, beta = 0){
  
  cur_pos <- start_state
  cur_action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon)
  
  # Follow policy, execute action, get reward.
  new_pos <- transition_model(cur_pos[1], cur_pos[2], cur_action, beta)
  reward <- reward_map[new_pos[1], new_pos[2]]
  new_action <- EpsilonGreedyPolicy(new_pos[1], new_pos[2], epsilon)
  
  repeat{
    if(reward!=0){
      # End episode.
      old_q <- q_table[cur_pos[1], cur_pos[2], cur_action]
      q_table[cur_pos[1], cur_pos[2], cur_action] <<- old_q + alpha*(reward - old_q)
      break
    }
    else{
      # Follow policy, execute action, get reward.
      new_pos2 <- transition_model(new_pos[1], new_pos[2], new_action, beta)
      reward2 <- reward_map[new_pos2[1], new_pos2[2]]
      new_action2 <- EpsilonGreedyPolicy(new_pos2[1], new_pos2[2], epsilon)

      if(reward2!=0){
        # End episode.
        old_q <- q_table[cur_pos[1], cur_pos[2], cur_action]
        q_table[cur_pos[1], cur_pos[2], cur_action] <<- old_q + alpha*(-1+gamma*reward2 - old_q)
        old_q <- q_table[new_pos[1], new_pos[2], new_action]
        q_table[new_pos[1], new_pos[2], new_action] <<- old_q + alpha*(reward2 - old_q)
        break
      }
      else{
        old_q <- q_table[cur_pos[1], cur_pos[2], cur_action]
        q_table[cur_pos[1], cur_pos[2], cur_action] <<- old_q + alpha*(-1+gamma*(-1)+gamma*gamma*q_table[new_pos2[1], new_pos2[2], new_action2] - old_q)
        cur_pos <- new_pos
        cur_action <- new_action
        new_pos <- new_pos2
        reward <- reward2
        new_action <- new_action2
      }
    }
  }
  
}

H <- 3
W <- 6

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[1,2:5] <- -10
reward_map[1,6] <- 10
# To avoid having to modify vis_environment, I take care of the reward -1 in the function SARSA.

q_table <- array(0,dim = c(H,W,4))

for(i in 1:5000)
  SARSA(start_state = c(1,1), epsilon = 0.5, gamma = 1, beta = 0)

vis_environment(5000, epsilon = 0.5, gamma = 1, beta = 0)