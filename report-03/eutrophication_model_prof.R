# Model 1: Basic P model - analogous to logistic growth
# dP/dt = growth - linear_loss = r*P*(1-P/K) - m*P

# Parameters
r <- 0.5    # growth rate
K <- 100     # carrying capacity (max biomass)
m <- 0.3    # linear loss rate (mortality, washout)
dt <- 0.1   # smaller dt for better stability
time <- seq(0, 100, by = dt)

# Initialize
P <- numeric(length(time))
P[1] <- 5   # initial phytoplankton biomass

# Explicit Euler integration loop
for (i in 2:length(time)) {
  # Logistic growth minus linear loss
  growth <- r * P[i-1] * (1 - P[i-1]/K)
  loss <- m * P[i-1]
  dP_dt <- growth - loss
  
  P[i] <- P[i-1] + dt * dP_dt
}
dev.new()
# Create dataframe and plot
results1 <- data.frame(time = time, P = P)
plot(results1$time, results1$P, type = "l", lwd = 2,
     xlab = "Time", ylab = "Phytoplankton Biomass (P)",
     main = "Basic Model: Single Stable State",
     ylim = c(-20, 50))
abline(h = K * (1 - m/r), col = "red", lty = 2, lwd = 2)
#legend("topright", legend = c("Dynamics", "Stable Equilibrium"),
       #col = c("black", "red"), lty = c(1, 2), lwd = 2)








# Model 2: Non-linear grazing creates alternative stable states
# dP/dt = r*P*(1-P/K) - (g*P^2)/(h^2 + P^2)
# where the second term is a type II/sigmoidal grazing function

# Parameters
r <- 0.3
K <- 15
g <- 1.2   # maximum grazing rate
h <- 2     # half-saturation constant
disturbance <- -0.0 # a number between -1 and 1.

dt <- 0.1
time <- seq(0, 400, by = dt)

# Test multiple starting points
initial_Ps <- c(0.5, 6, 12)  # clear, intermediate, turbid
colors <- c("blue", "green", "red")
dev.new()
# Set up plot
plot(NA, xlim = range(time), ylim = c(0, 15),
     xlab = "Time", ylab = "Phytoplankton (P)",
     main = "Alternative Stable States: Priority Effect in Lakes")

# Loop over different starting conditions
for (run in 1:length(initial_Ps)) {
  P <- numeric(length(time))
  P[1] <- initial_Ps[run]
  
  # Euler integration
  for (i in 2:length(time)) {
    growth <- r * P[i-1] * (1 - P[i-1]/K)
    # Non-linear grazing loss
    grazing <- (g * P[i-1]^2) / (h^2 + P[i-1]^2)
    dP_dt <- growth - grazing
    
    
    P[i] <- P[i-1] + dt * dP_dt 
    if(i==700)P[i]<-P[i]+P[i]*disturbance
    
  }
  print(P[700])
  lines(time, P, col = colors[run], lwd = 2)
}








# Create function for dP/dt
dP_dt_func <- function(P_val) { # creates a function to predict how the speed of growth will behave in function of P
  growth <- r * P_val * (1 - P_val/K)
  grazing <- (g * P_val^2) / (h^2 + P_val^2)
  return(growth - grazing)
}
r=0.5
h=1.5
g=1.5
# Calculate over P range
P_vals <- seq(0, 15, length.out = 1000) 
dP_vals <- sapply(P_vals, dP_dt_func)

#abs_dP<-abs(dP_vals)

#closest_index<-order(abs_dP)[1:5]
  
#equilibria<-P_vals[closest_index]
dev.new()
# Plot the 1D phase portrait
plot(P_vals, dP_vals, type = "l", lwd = 2,
     xlab = "Phytoplankton (P)", ylab = "dP/dt",
     main = "1D Phase Portrait: Find Equilibria")
abline(h = 0, col = "gray", lty = 2)

# Mark equilibria (visually or by calculation)
# Zeros occur where dP/dt crosses 0: ~0.8, ~4.2, ~11.5
equilibria <- c(1.29838, 6, 7.900719)
points(equilibria, rep(0, 3), col = c("green","green", "red", "green","green"), 
       pch = 19, cex = 1.5)
#text(equilibria, c(0.2, 0.2, -0.2), 
 #    labels = c("Stable Clear", "Unstable", "Stable Turbid"),
  #   pos = c(4, 4, 2))

# Model 3: Hysteresis simulation - slowly changing nutrient input
# K = scaling * N_input

simulate_lake <- function(N_input, P_initial, time_horizon = 200) {
  # Parameters
  r <- 0.1
  scaling <- 0.065  # how nutrients translate to carrying capacity
  g <- 0.085
  h <- 0.035
  
  # Calculate effective carrying capacity
  K_eff <- scaling * N_input
  
  # Integration
  dt <- 0.1
  time <- seq(0, time_horizon, by = dt)
  P <- numeric(length(time))
  P[1] <- P_initial
  
  for (i in 2:length(time)) {
    growth <- r * P[i-1] * (1 - P[i-1]/K_eff)
    grazing <- (g * P[i-1]^2) / (h^2 + P[i-1]^2)
    dP_dt <- growth - grazing
    
    P[i] <- P[i-1] + dt * dP_dt
  }
  
  # Return final (equilibrium) state
  return(tail(P, 1))
}

# Ramp nutrient input UP (starting from clear state)
N_up <- seq(5, 200, by = 0.25)
P_final_up <- numeric(length(N_up))

for (j in 1:length(N_up)) {
  P_final_up[j] <- simulate_lake(N_input = N_up[j], P_initial = 1)
}

P_final_up
# Ramp nutrient input DOWN (starting from turbid state)
N_down <- seq(200, 5, by = -0.25)
P_final_down <- numeric(length(N_down))

for (j in 1:length(N_down)) {
  P_final_down[j] <- simulate_lake(N_input = N_down[j], P_initial = 14)
}

dev.new()
# Plot the hysteresis loop
plot(N_up, P_final_up, type = "l", lwd = 3, col = "blue",
    xlab = "Nutrient Input (N)", ylab = "Equilibrium Phytoplankton (P)",
    main = "Hysteresis in Lake Eutrophication",
    ylim = c(0, 16), xlim = c(5,200))

lines(N_down, P_final_down, lwd = 3, col = "red")



# Mark tipping points (approximate)
# abline(v = 12, lty = 3, col = "darkgray")  # Collapse point
# abline(v = 8, lty = 3, col = "darkgray")   # Recovery point
# text(c(12, 8), c(15, 15), 
#      labels = c("Collapse", "Recovery"), 
#      pos = c(4, 2))

legend("topleft", 
      legend = c("Increasing Nutrients", "Decreasing Nutrients"),
      col = c("blue", "red"), lwd = 3)
