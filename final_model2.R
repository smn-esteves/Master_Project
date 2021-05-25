# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)

# MODEL INPUTS:

initial_state_values <- c(S = 1000000-3,
                          I = 1,        
                          E = 2,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
parameters <- c(beta = 0.61,     # the infection rate in units of years^-1
                delta = 0.09,     # the latency period in units of years^-1
                c_s = 0.4,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.4,# the reduction in the infectivity of vaccinated infected people  
                u = 0.1,#death rate in units of years^-1
                a = 0.7, #cull due to infection in units of years^-1
                b = 0.3, #birth rate in units of years^-1
                vc = 0.8) # vaccine coverage    

# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 20, by =1/365)#from 0 to 20 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv
    lambda <- beta * E/N + c_i * beta * Ev/N 
    # the Ev compartment gets c_i times less infected than the E compartment
   
    
    # The differential equations
    dS <- -lambda * S - u * S + b * N - vc * S           
    dE <- lambda * S - delta * E - u * E 
    dI <- delta * E - a * I - u * I  
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S              
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev
    dIv <- delta * Ev - a * Iv - u * Iv
    
    return(list(c(dS, dI, dE, dSv, dIv,dEv))) 
  })
  
}

# MODEL OUTPUT:

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = vaccine_model,
                            parms = parameters))

# PLOT THE OUTPUT

# turn output dataset into long format
output_long <- melt(as.data.frame(output), id = "time")

# Adding a column for the prevalence proportion to the long-format output
output_long$prevalence <- output_long$value/sum(initial_state_values)

# Plot the number in each compartment over time
ggplot(data = output_long,                                               
       aes(x = time, y = prevalence, colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("Proportion of the population") +
  labs(title = paste("Leaky vaccine with coverage of", 80, "%"), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2")