# Poisson distribution
require(deSolve)
require(ggplot2)

db_bTB_PS

# INPUT
initial_state_values <- c(S = 100000,  
                          E = 1,       
                          I = 0,
                          Sv=0,
                          Ev=0,
                          Iv=0)

times <- seq(from = 0, to = 50, by = 1)

# INPUT
vaccine_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {    
    
    # Defining lambda as a function of beta and E:
    N <- S + E + I + Sv + Ev + Iv
    lambda <- beta * I/N + beta + Iv/N
  
    
    # the Ev compartment gets c_i times less infected than the E compartment
    
    
    # The differential equations
    dS <- -(beta * E/N) * S           
    dE <- (beta * E/N) * S - delta * E 
    dI <- delta * E
    dSv <- - (beta * E/N) * Sv              
    dEv <-  (beta * E/N) * Sv - delta * Ev 
    dIv <- delta * Ev  
    
    return(list(c(dS, dE, dI, dSv, dEv, dIv))) 
  })
  
}
# DISTANCE FUNCTION

loglik_function <- function(parameters, dat) {   # takes as inputs the parameter values and dataset
  
  beta <- parameters[1]    # extract and save the first value in the "parameters" input argument as beta
  delta <- parameters[2]   # extract and save the second value in the "parameters" input argument as gamma
  
  # Simulate the model with initial conditions and timesteps defined above, and parameter values from function call
  output <- as.data.frame(ode(y = initial_state_values, 
                              times = times, 
                              func = vaccine_model,
                              parms = c(beta = beta,       # ode() takes the values for beta and gamma extracted from
                                        delta = delta)))   # the "parameters" input argument of the loglik_function()
  
  # Calculate log-likelihood using code block 4 from the previous etivity, accounting for the reporting rate of 60%:
  LL <- sum(dpois(x = dat$I, lambda = 0.01 * output$I[output$time %in% dat$time], log = TRUE))
  
  return(LL) 
}
# OPTIMISATION

optim(par = c(0.01, 0.01),           # starting values for beta and gamma - you should get the same result no matter 
      # which values you choose here
      fn = loglik_function,        # the distance function to optimise
      dat = db_bTB_PS,         # the dataset to fit to ("dat" argument is passed to the function specified in fn)
      control = list(fnscale=-1))  # tells optim() to look for the maximum number instead of the minimum (the default)


# Simulate the model with the estimated best-fitting parameter values
parameters <- c(beta = ,
                delta = )

output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = vaccine_model,
                            parms = parameters))

# PLOT OF THE MODEL FIT

ggplot() +
  geom_line(data = output, aes(x = time, y = I)) +                              
  geom_point(data = db_bTB_PS, aes(x = time, y = I, colour="I")) + 
  xlab("Time (days)")+                                              
  ylab("Number of bovines") +                                 
  labs(title = paste("Model fit to the epidemic curve with beta =", parameters["beta"], 
                     "and delta =", parameters["delta"]), colour = "")