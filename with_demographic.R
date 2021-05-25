# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)

# MODEL INPUTS:

# Vector storing the initial number of animals in each compartment (at timestep 0)
initial_state_values <- c(S = 999999,  # the whole population we are modelling is susceptible to infection
                          E = 1, 
                          I = 0)       

# Vector storing the parameters describing the transition rates in units of days^-1
parameters <- c(beta = 0.4,     # the infection rate in units of years^-1
                delta = 0.2,     # the latency period in units of years^-1
                u = 0.4,#death rate
                a = 0.3, #cull due to infection
                b = 0.3) # birth rate

# TIMESTEPS:

# Vector storing the sequence of timesteps to solve the model at
times <- seq(from = 0, to = 100, by = 1/365)   # from 0 to 100 years in daily intervals

# SEID MODEL FUNCTION: 

# The model function takes as input arguments (in the following order): time, state and parameters
seid_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {  # tell R to unpack variable names from the state and parameters inputs
    
    N <- S + E + I
    lambda <- beta * I/N
    
    # The differential equations
    dS <- -lambda * S - u * S + b * N            
    dE <- lambda * S - delta * E - u * E 
    dI <- delta * E - a * I - u * I
    
    # Return the number of animals in the S, E and I compartments at each timestep 
    # (in the same order as the input state variables)
    return(list(c(dS, dE, dI))) 
  })
  
}

# MODEL OUTPUT (solving the differential equations):

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = seid_model,
                            parms = parameters))

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

#R0= beta/(a+u) * (delta+u)/u