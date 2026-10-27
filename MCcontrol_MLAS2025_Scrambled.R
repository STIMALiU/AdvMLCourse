# By Jose M. Peña.
# For teaching purposes.
# jose.m.pena@liu.se.

#####################################################################################################
# MC Control
#####################################################################################################

# install.packages("ggplot2")
# install.packages("vctrs")
library(ggplot2)

arrows <- c("^", ">", "v", "<")
action_deltas <- list(c(1,0), # up
                      c(0,1), # right
                      c(-1,0), # down
                      c(0,-1)) # left

vis_environment <- function(iterations=0, epsilon = 0.5, gamma = 0.95, beta = 0){
  
  # Visualize an environment with rewards. 
  # Q-values for all actions are displayed on the edges of each tile.
  # The (greedy) policy for each state is also displayed.
  # 
  # Args:
  #   episodes, epsilon, gamma, beta (optional): for the figure title.
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
          ggtitle(paste("Q-table after ",iterations," episodes\n",
                        "(epsilon = ",epsilon,"gamma = ",gamma,", beta = ",beta,")")) +
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

MC_episode <- function(start_state, epsilon = 0.5, beta = 0){
  
  # Generate an episode. The agent should move around in the environment using the given
  # transition model and store the episode in the corresponding global variables. The episode
  # ends when the agent reaches a terminal state.
  # 
  # Args:
  #   start_state: array with two entries, describing the starting position of the agent.
  #   epsilon (optional): probability of acting randomly.
  #   beta (optional): slipping factor.
  #   reward_map (global variable): a HxW array containing the reward given at each state.
  # 
  # Returns:
  #   episode_pos1, episode_pos2, episode_action, episode_reward (global variables): Recall that R passes
  #   arguments by value. So, a global variable should be modified with the superassigment operator <<-.
  
# Scrambled code.
 
} 
break
repeat{ 
if(reward!=0)
cur_pos <- new_pos
cur_pos <- start_state
reward <- reward_map[new_pos[1], new_pos[2]]
episode_pos1 <<- c(episode_pos1, cur_pos[1])
episode_pos2 <<- c(episode_pos2, cur_pos[2])
episode_action <<- c(episode_action, action)
episode_reward <<- c(episode_reward, reward)
action <- EpsilonGreedyPolicy(cur_pos[1], cur_pos[2], epsilon)   
new_pos <- transition_model(cur_pos[1], cur_pos[2], action, beta)
    
# End scrambled code.
  
}

MC_control <- function(n_episode, gamma = 0.95){
  
  # Update the Q-values for the state-action pairs in an episode.
  # 
  # Args:
  #   n_episode: episode number.
  #   gamma (optional): discount factor.
  #   episode_pos1, episode_pos2, episode_action, episode_reward (global variables): lists containing
  #   respectively the X and Y coordinates of the states, the actions and rewards in the episode. 
  # 
  # Returns:
  #   q_table (global variable): Recall that R passes arguments by value. So, q_table being
  #   a global variable can be modified with the superassigment operator <<-.

# Scrambled code.
  
}
}
}
G <- 0
if(t>1){
if(! repeated){
repeated <- FALSE
G <- gamma * G + episode_reward[t]
repeated <- any(foo1 & foo2 & foo3)
for(t in length(episode_reward):1){
foo1 <- episode_pos1[1:t-1]==episode_pos1[t]
foo2 <- episode_pos2[1:t-1]==episode_pos2[t]
foo3 <- episode_action[1:t-1]==episode_action[t]
q_table[episode_pos1[t], episode_pos2[t], episode_action[t]] <<- 
  (q_table[episode_pos1[t], episode_pos2[t], episode_action[t]] * (n_episode - 1) + G) / n_episode
    
# End scrambled code.
  
}

#####################################################################################################
# MC Control Environments
#####################################################################################################

# Environment A (learning)

H <- 5
W <- 7

reward_map <- matrix(0, nrow = H, ncol = W)
reward_map[3,6] <- 10
reward_map[2:4,3] <- -1

q_table <- array(0,dim = c(H,W,4))

vis_environment()

for(i in 1:1000){
  episode_pos1 <- NULL
  episode_pos2 <- NULL
  episode_action <- NULL
  episode_reward <- NULL
  
  MC_episode(start_state = c(3,1))
  MC_control(n_episode = i)
  
  if(any(i==c(10,100,1000,1000)))
    vis_environment(i)
}
