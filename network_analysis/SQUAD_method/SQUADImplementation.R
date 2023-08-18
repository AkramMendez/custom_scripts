# This script exemplifies the implementation of regulatory networks using SQUAD.

library(BoolNet)
library(deSolve)

# Simulation of the regulatory network defined in figure 1 as a discrete model using BoolNet R package

# load the network
net<-loadNetwork("cartoonNetwork.R")

# define a initial state (A,B,C) = (0,1,1)
initialState<-generateState(net,specs=c("A"=0,"B"=1,"C"=1))

# simulate a discrete trajectory 
plotSequence(net,initialState)

#########################################################################################################

# Transformation of the regulatory network model as defined in figure 1 as a continuous model using SQUAD

# define network parameters
parameters <- c(h = 50,gamma = 1)

# define SQUAD generic function
SQUAD<-function(w,gamma,h){
  val<- ((-exp(0.5*h) + exp(-h*(interaction(x)-0.5))) / ((1-exp(0.5*h)) * (1+exp(-h*(interaction(x)-0.5))))) - (gamma*x)
  return(val)
}

# define ODE SQUAD system using deSolve R package
squadInteractions<-function(times,state,parameters) {
  with(as.list(c(state,parameters)),{
    w_A <- A
    w_B <- min(A,1-C)
    w_C <- 1 - max(A,B)
    A <- SQUAD(w_A,gamma,h)
    B <- SQUAD(w_B,gamma,h)
    C <- SQUAD(w_C,gamma,h)
    return(list(c(A,B,C)))
  })
}

# define interval for numeric integration
times<-seq(0,10,0.01)

# simulate the same trajectory but in a continuous manner
result<-ode(y=initialState,times=times,func=squadInteractions,parms = parameters,atol=10e-6, rtol=10e-6)
result



