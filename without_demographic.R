# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)

# MODEL INPUTS:

# Vector storing the initial number of animals in each compartment
# (at timestep 0)
initial_state_values <- c(S = 999999
                          ,  # the whole population we're modelling
                          # is susceptible to infection
                          E = 1,       
                          I = 0)       
# population

# Vector storing the parameters describing the transition rates in
# units of years^-1
parameters <- c(beta = 0.0276*365,      # the infection rate
                delta = 0.0164*365)   # the rate of recovery, which acts on 
# those infected

# TIMESTEPS:

# Vector storing the sequence of timesteps to solve the model at
times <- seq(from = 0, to = 20, by = 1)   
# from 0 to 60 days in daily intervals

# SIR MODEL FUNCTION: 

# The model function takes as input arguments (in the following order):
# time, state and parameters
sie_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {  # tell R to unpack variable names 
    # from the state and parameters inputs    
    
    # New: calculating the total population size N        
    N <- S+E+I # (the sum of the number of bovines in each compartment)
    
    # New: defining lambda as a function of beta and I:
    lambda <- beta * I/N
    # Another option is simply replacing lambda with this
    # expression in the differential equations below
    
    # The differential equations
    dS <- -lambda * S               
    dE <- lambda * S - delta * E    
    dI <- delta * E               
    
    return(list(c(dS, dE, dI))) 
  })
  
}

# MODEL OUTPUT (solving the differential equations):

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sie_model,
                            parms = parameters))

# Plotting the output
output_long <- melt(as.data.frame(output), id = "time")                  # turn output dataset into long format

# Adding a column for the prevalence proportion to the long-format output
output_long$prevalence <- output_long$value/sum(initial_state_values)

# Plot the prevalence proportion
ggplot(data = output_long,                                               # specify object containing data to plot
       aes(x = time, y = prevalence, colour = variable, group = variable)) +  # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                   # add label for x axis
  ylab("Prevalence") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title
       title = "Prevalence of susceptibility, exposed and infected over time")   # add plot title