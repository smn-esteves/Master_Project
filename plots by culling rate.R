# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.77, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g1<- ggplot(data = output_long,                                               
       aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("a) Test and slaughter,", alpha, "= 0.77 and vc = 0")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))


############

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.77, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.5,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g2<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("b) Test and slaughter and vaccination,", alpha, "= 0.77 and vc = 0.5")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

###########

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.77, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.7,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g3<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("c) Test and slaughter and vaccination,", alpha, "= 0.77 and vc = 0.7")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

###############
# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.77, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.9,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g4<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("d) Test and slaughter and vaccination,", alpha, "= 0.77 and vc = 0.9")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))


# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.77, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 1,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g5<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("e) Test and slaughter and vaccination,", alpha, "= 0.77 and vc = 1")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

library("gridExtra")
library("egg")
grid.arrange(g1, g2, g3, g4, g5, nrow = 3)

############################################################################################################################################################################################################################################

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.9, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g1<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("a) Test and slaughter,", alpha, "= 0.9 and vc = 0")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))


############

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.9, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.5,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g2<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("b) Test and slaughter and vaccination,", alpha, "= 0.9 and vc = 0.5")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

###########

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.9, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.7,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g3<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("c) Test and slaughter and vaccination,", alpha, "= 0.9 and vc = 0.7")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

###############
# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.9, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.9,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g4<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("d) Test and slaughter and vaccination,", alpha, "= 0.9 and vc = 0.9")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))


# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 0.9, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 1,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g5<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("e) Test and slaughter and vaccination,", alpha, "= 0.9 and vc = 1")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

library("gridExtra")
library("egg")
grid.arrange(g1, g2, g3, g4, g5, nrow = 3)

############################################################################################################################################################################################################################################

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 1, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g1<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("a) Test and slaughter,", alpha, "= 1 and vc = 0")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))


############

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 1, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.5,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g2<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("b) Test and slaughter and vaccination,", alpha, "= 1 and vc = 0.5")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

###########

# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 1, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.7,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g3<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("c) Test and slaughter and vaccination,", alpha, "= 1 and vc = 0.7")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2")+
  theme(plot.title = element_text(vjust = 3))

###############
# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 1, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 0.9,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g4<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("d) Test and slaughter and vaccination,", alpha, "= 1 and vc = 0.9")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))


# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)


# MODEL INPUTS:

initial_state_values <- c(S = 1224000- 4672,
                          E = 1000,        
                          I = 3672,        
                          Sv = 0,      
                          Ev = 0,
                          Iv = 0)      

# Parameters
#R0=0.42

parameters <- c(beta = 0.01,     # the infection rate in units of years^-1  5.2
                delta = 0.01*365,     # the latency period in units of years^-1 
                c_s = 0.39,       # the reduction in the force of infection
                # acting on those vaccinated
                c_i = 0.39,# the reduction in the infectivity of vaccinated infected bovines  
                u = 1/5,#death rate in units of years^-1 
                a = 1, #testing rate in units of years^-1
                b = 1/5, #birth rate in units of years^-1
                vc = 1,  # vaccination rate
                w = 0.0009 ) #wildife infection rate 0.131



# TIMESTEPS:

# Sequence of timesteps to solve the model at
times <- seq(from = 0, to = 15, by =0.1)#from 0 to 15 years, daily intervalS
# MODEL FUNCTION: 

vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv 
    lambda <- beta * I/N + c_i * beta * Iv/N 
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -lambda * S - u * S - S * w - vc * S + (b * N * (1-vc))+ a * I + a * Iv
    dE <- lambda * S - delta * E - u * E - vc * E 
    dI <- delta * E - a * I - u * I + S * w 
    dSv <- -c_s * lambda * Sv - u * Sv + vc * S  + b * N * vc - Sv * w            
    dEv <- c_s * lambda * Sv - delta * Ev - u * Ev + vc * E 
    dIv <- delta * Ev - a * Iv - u * Iv + Sv * w 
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
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
g5<- ggplot(data = output_long,                                               
            aes(x = time, y = log(prevalence), colour = variable, group = variable)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                   
  ylab("log (Proportion of the population)") +
  labs(title = expression(paste("e)Test and slaughter and vaccination,", alpha, "= 1 and vc = 1")), 
       colour = "Compartment") +
  scale_colour_brewer(palette = "Set2") +
  theme(plot.title = element_text(vjust = 3))

library("gridExtra")
library("egg")
grid.arrange(g1, g2, g3, g4, g5, nrow = 3)
