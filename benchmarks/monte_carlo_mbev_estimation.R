


source("R/mbev_dist.R")
source("R/mbev_estimation.R")












# create design
library(tibble)
library(dplyr)
theta_tbl <- tibble(
  theta_id = paste0("Theta", 1:6),
  dep  = c(1,   1,   1/2,   1,   1,   1/2.5),
  mu1    = c(0,   0,   0,   0,   0,   0),
  sigma1 = c(1,   1,   1,   1,   1,   1),
  xi1    = c(0.5, 0.5, 1,  -1,  0.5, 0),
  delta1 = c(0,   0,   0,   1,   1,   0.8),
  mu2    = c(0,   0,   0,   0,   0,   0),
  sigma2 = c(1,   1,   1,   1,   1,   1),
  xi2    = c(0.5, 0.5, 0,  -1,  0.5, 0),
  delta2 = c(0,   1,   1,   2,   1,   9)
)
n_vals <- c(50, 100, 500)
Design <- tidyr::crossing(
  n = n_vals,
  theta_tbl
) %>%
  select(
    n, dep,
    mu1, sigma1, xi1, delta1,
    mu2, sigma2, xi2, delta2
  )




Generate <- function(condition, fixed_objects) {
  dat <- with(condition, rmbev(n = n, mu1 = mu1, mu2 = mu2, delta1 = delta1, delta2 = delta2, sigma1 = sigma1, sigma2 = sigma2, xi1 = xi1, xi2 = xi2, dep = dep)  )
  dat
} 
Analyse <- function(condition, dat, fixed_objects) {
  est = mbev_estimation(dat)
  if(est$value == 0)
    return(rep(NA,9))
  ret = as.vector(est$par)
  names(ret) = c("mu1", "mu2", "delta1", "delta2", "sigma1", "sigma2", "xi1", "xi2", "dep") 
  return(ret)
}


analysis <- function(condition, dat, fixed_objects) {
  est <- mbev_estimation(dat)
  
  # detect degenerate / penalty solution
 # if (!is.finite(est$value) || est$value <= 0) {
  #  return(rep(NA_real_, 9))  # one NA per estimand
  #}
  
  c(
    mu1 = est$par["mu1"],
    mu2 = est$par["mu2"],
    delta1 = est$par["delta1"],
    delta2 = est$par["delta2"],
    sigma1 = est$par["sigma1"],
    sigma2 = est$par["sigma2"],
    xi1 = est$par["xi1"],
    xi2 = est$par["xi2"],
    dep = est$par["dep"]
  )
}






Summarise <- function(condition, results, fixed_objects) {
  # assuming your Design object columns match these names)
  true_mu1 <- condition$mu1
  true_mu2 <- condition$mu2
  true_delta1 <- condition$delta1
  true_delta2 <- condition$delta2
  true_sigma1 <- condition$sigma1
  true_sigma2 <- condition$sigma2
  true_xi1 <- condition$xi1
  true_xi2 <- condition$xi2
  true_dep <- condition$dep
  
  # Return a named vector of the summary statistics (bias and RMSE)
  ret <- c(
    bias_mu1 = bias(results[, "mu1"], parameter = true_mu1),
    bias_mu2 = bias(results[, "mu2"], parameter = true_mu2),
    bias_delta1 = bias(results[, "delta1"], parameter = true_delta1),
    bias_delta2 = bias(results[, "delta2"], parameter = true_delta2),
    bias_sigma1 = bias(results[, "sigma1"], parameter = true_sigma1),
    bias_sigma2 = bias(results[, "sigma2"], parameter = true_sigma2),  
    bias_xi1 = bias(results[, "xi1"], parameter = true_xi1),
    bias_xi2 = bias(results[, "xi2"], parameter = true_xi2),
    bias_dep = bias(results[, "dep"], parameter = true_dep),
    RMSE_mu1 = RMSE(results[, "mu1"], parameter = true_mu1),
    RMSE_mu2 = RMSE(results[, "mu2"], parameter = true_mu2),
    RMSE_delta1 = RMSE(results[, "delta1"], parameter = true_delta1),
    RMSE_delta2 = RMSE(results[, "delta2"], parameter = true_delta2),
    RMSE_sigma1 = RMSE(results[, "sigma1"], parameter = true_sigma1),
    RMSE_sigma2 = RMSE(results[, "sigma2"], parameter = true_sigma2),  
    RMSE_xi1 = RMSE(results[, "xi1"], parameter = true_xi1),
    RMSE_xi2 = RMSE(results[, "xi2"], parameter = true_xi2),
    RMSE_dep = RMSE(results[, "dep"], parameter = true_dep) 
  )
  
  return(ret)
}

Design2 = Design[c(6,12,18),]
Final <- SimDesign::runSimulation(design=Design2, replications=10,
                                  generate=Generate, analyse=Analyse, summarise=Summarise,
                                  progress = FALSE, verbose = FALSE)

t(Final)  # see results


# SAVE RESULTS
saveRDS(
  Final,
  file = paste0(
    "benchmarks/MC_delta_9_rep_500_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".rds"
  )
)


# does bgev works to estimate mu = 0, sigma = 1, xi = 0, delta = 9
library(bgev)
source("benchmarks/bgevmle_new.R")
for(i in 1:10){
  print("\n\n=====================")
  x = bgev::rbgev(50, mu = 0, sigma = 1, xi = 0, delta = 9)
  bgev.mle(x)
  #bgev.mle.new(x, lower = c(-20, 0.001, -20, -0.999), upper = c(20, 10, 20, 20), verbose = TRUE)
}

# one run of mbev_estimation
for(i in 1:10){
  dat = rmbev(n =50, mu1 = 0, mu2 = 0, delta1 = 0.8, delta2 = 9, sigma1 = 1, sigma2 = 1, xi1 = 0, xi2 = 0, dep = 0.4)
  est = mbev_estimation(dat)
  print("=====================")
  print(est$par)
  print(est$value)
}

# bgev estimate of second coordinate of mbev data, handle cases with errors in initial estimates
for(i in 1:10){
  print("=====================")
  dat = rmbev(n = 500, mu1 = 0, mu2 = 0, delta1 = 0.8, delta2 = 9, sigma1 = 1, sigma2 = 1, xi1 = 0, xi2 = 0, dep = 0.4)
  bgev.mle.new(as.vector(dat[,2]), lower = c(-20, 0.001, -20, -0.999), upper = c(20, 10, 20, 20), verbose = TRUE)
}

