


# true parameters
mu1    <- 1
mu2    <- -1
delta1 <- 0.5
delta2 <- 9
sigma1 <- 1.5
sigma2 <- 0.5
xi1    <- 0.5
xi2    <- 1
dep    <- 0.5   

pars_true <- c(mu1, mu2, delta1, delta2,
               sigma1, sigma2, xi1, xi2, dep)
set.seed(1)
# simulate data
Y <- rmbev(
  n = 1000,
  mu1 = mu1, mu2 = mu2,
  delta1 = delta1, delta2 = delta2,
  sigma1 = sigma1, sigma2 = sigma2,
  xi1 = xi1, xi2 = xi2,
  dep = dep
)

# objective: estimate (mu1, mu2, dep), others fixed
f <- function(pars) {
  val = -mbev_log_likelihood(
    Y,
    pars = pars
  )
  if (!is.finite(val)) return(1e99)
  val
}

grid_size = 10
# bounds
par_lower = c(-rep(grid_size,2), -0.99, -0.99, 0.01,0.01, -rep(grid_size,2), 0.1)
par_upper = c(rep(grid_size,8), 0.9)


# get good starting values via DEoptim
pars_estimated <- DEoptim::DEoptim(
  fn = f,
  lower = par_lower,
  upper = par_upper,
  control = DEoptim::DEoptim.control(
    itermax = 100,
    trace = FALSE
  )
)$optim$bestmem
cbind(pars_true,round(pars_estimated,1))


plot_profile <- function(j, grid, theta_hat, Y, theta_true, obs = "") {
  ll <- sapply(grid, function(x) {
    theta <- theta_hat
    theta[j] <- x
    -mbev_log_likelihood(Y, theta)
  })
  plot(grid, ll, type = "l",
       main = paste("Profile likelihood for par", j, "\n",obs),
       ylab = "log-likelihood", xlab = "")
  abline(v = theta_true, col = "red")
}
pars_new = pars_estimated
pars_new[6] = pars_true[6]
pars_new[2] = pars_true[2]
par(mfrow = c(3,3))
for(j in 1:9)
  plot_profile(
    j = j,
    grid = seq(par_lower[j], par_upper[j], length.out = 50),
    theta_hat = pars_estimated,
    Y = Y,
    theta_true = pars_true[j],
    ""
  )
#(mu1, mu2, delta1, delta2, sigma1, sigma2, xi1, xi2, dep)


# par 8 is maybe correlated with pars 6 and 2 (fixing those two it becomes identifiable)

