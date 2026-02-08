dt <- 0.05
t <- seq(0, 50, by = dt)
r1 <- 0.5
r2 <- 0.5
n0_1 <- 50
n0_2 <- 50

# Carrying capacities
k1 <- 400
k2 <- 250

# Competition coefficients for each scenario
alpha_1a <- 0.5
alpha_2a <- 2

alpha_1b <- 2
alpha_2b <- 0.5

alpha_1c <- 2.5
alpha_2c <- 1.5

alpha_1d <- 0.5
alpha_2d <- 0.5

# Run simulations for each scenario
df_a <- N_asym(dt, t, r1, r2, k1, k2, n0_1, n0_2, alpha_1a, alpha_2a)
df_b <- N_asym(dt, t, r1, r2, k1, k2, n0_1, n0_2, alpha_1b, alpha_2b)
df_c <- N_asym(dt, t, r1, r2, k1, k2, n0_1, n0_2, alpha_1c, alpha_2c)
df_d <- N_asym(dt, t, r1, r2, k1, k2, n0_1, n0_2, alpha_1d, alpha_2d)

# Create 2x2 plot layout
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))

# Plot (a): Species 1 wins
plot(df_a$Time, df_a$N1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(k1, k2)),
  xlab = "Time", ylab = "Population Size", main = "(a) Species 1 wins")
lines(df_a$Time, df_a$N2, col = "red", lwd = 2)
legend("right", legend = c("Species 1", "Species 2"), col = c("blue", "red"), 
  lwd = 2, cex = 0.8)

# Plot (b): Species 2 wins
plot(df_b$Time, df_b$N1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(k1, k2)),
  xlab = "Time", ylab = "Population Size", main = "(b) Species 2 wins")
lines(df_b$Time, df_b$N2, col = "red", lwd = 2)
legend("right", legend = c("Species 1", "Species 2"), col = c("blue", "red"), 
  lwd = 2, cex = 0.8)

# Plot (c): Stable coexistence
plot(df_c$Time, df_c$N1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(k1, k2)),
  xlab = "Time", ylab = "Population Size", main = "(c) Stable coexistence")
lines(df_c$Time, df_c$N2, col = "red", lwd = 2)
legend("right", legend = c("Species 1", "Species 2"), col = c("blue", "red"), 
  lwd = 2, cex = 0.8)

# Plot (d): Unstable equilibrium
plot(df_d$Time, df_d$N1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(k1, k2)),
  xlab = "Time", ylab = "Population Size", main = "(d) Unstable equilibrium")
lines(df_d$Time, df_d$N2, col = "red", lwd = 2)
legend("right", legend = c("Species 1", "Species 2"), col = c("blue", "red"), 
  lwd = 2, cex = 0.8)

mtext("Lotka-Volterra Competition Dynamics Over Time", outer = TRUE, cex = 1.2)
