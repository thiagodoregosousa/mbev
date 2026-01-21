require(bgev)

test_that("likelihood behaviour at invalid parameters", {
  
  # invalid delta1
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,-2,0,1,1,0.5,0.5,1)), 1e99)
  
  # invalid delta2
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,2,-1,1,1,0.5,0.5,1)), 1e99)
  
  # invalid dep
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,2,-1,1,1,0.5,0.5,-0.1)), 1e99)
  
  # invalid dep (depends on my parametrization)
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,2,-1,1,1,0.5,0.5,1)), 1e99) 
  
})


test_that("likelihood finite at truth", {
  
  Y <- rmbev(n = 200, mu1=0, mu2=0, delta1=0.5, delta2=0.2,
             sigma1=1, sigma2=2, xi1=0.1, xi2=0.2, dep=0.5)
  
 pars_true <- c(0,0,0.5,0.2,1,2,0.1,0.2,0.5)
  
 log_likelihood_Y = mbev_log_likelihood(Y, pars_true)
  
 expect_true(is.finite(log_likelihood_Y))
  
 expect_true(log_likelihood_Y != 1e99)
  
})


test_that("likelihood worsen when parameters are perturbed", {
  
  Y <- rmbev(n = 200, mu1=0, mu2=0, delta1=0.5, delta2=0.2,
             sigma1=1, sigma2=2, xi1=0.1, xi2=0.2, dep=0.5)
  pars_true <- c(0,0,0.5,0.2,1,2,0.1,0.2,0.5)
  log_likelihood_pars_true = mbev_log_likelihood(Y, pars_true)
  
  # pertub each of the 9 parameters by 0.4
  for( i in 1:9){
    pars_true_perturbed = pars_true
    pars_true_perturbed[i] = pars_true_perturbed[i] + 0.4
    log_likelihood_pars_true_perturbed = mbev_log_likelihood(Y, pars_true_perturbed)
    
    expect_gt(log_likelihood_pars_true, log_likelihood_pars_true_perturbed )
  }
  
})



test_that("Estimating ONLY dependency parameter using log-likelihood works with rnd starting values in optim", {
  
  mc_dep <- function(n = 200, dep_true = 0.5) {
    Y <- rmbev(n = n, dep = dep_true)
    
    f <- function(d)
      -mbev_log_likelihood(
        Y,
        pars = c(0,0,0,0,1,1,0.5,0.5,d)
      )
    
    
    opt = try(optim(runif(1,0,1), f, method = "L-BFGS-B", lower = 0.01, upper = 0.99), silent = TRUE)
    if(class(opt) == "try-error")
      return(c(1e99, 1))
    
    c(dep_hat = opt$par, conv = opt$convergence)
  }
  
  mc_dep_results = replicate(10, mc_dep(), simplify = FALSE)
  
  conv    <- sapply(mc_dep_results, function(x) x["conv"])
  converged_indexes = which(conv == 0)
  dep_hat <- sapply(mc_dep_results[converged_indexes], function(x) x["dep_hat"])
  
  # 1) at least 80% converged
  conv_rate <- length(converged_indexes)/length(conv)
  expect_gte(conv_rate, 0.8)
  
  # 2) mean close to true value 0.5
  expect_equal(mean(dep_hat), 0.5, tolerance = 0.05)
  
  # 3) std close to expected dispersion (0.2)
  expect_lt(sd(dep_hat), 0.1)
  
})





test_that("Estimating ONLY dependency and means using log-likelihood works with rnd starting values in optim", {
  
 
   mc <- function(n = 200, mu1_true = 1, mu2_true = -2, dep_true = 0.5) {
    Y <- rmbev(n = n, mu1 = mu1_true, mu2 = mu2_true, dep = dep_true)
    
    f <- function(mus_and_d)
      -mbev_log_likelihood(
        Y,
        pars = c(mus_and_d[1], mus_and_d[2], 0,0, 1,1, 0.5,0.5, mus_and_d[3])
      )
    
    pars.lower = c(-5,-5,0.01)
    pars.upper = c(5,5,0.99)
    
    # get good starting values
    pars.start = DEoptim::DEoptim(fn = f, lower = pars.lower, upper = pars.upper, control = DEoptim::DEoptim.control(itermax = 5, 
                                                                                                      trace = FALSE))$optim$bestmem
    
    opt = try(suppressWarnings(Rsolnp::solnp(pars = pars.start, fun = f, LB = pars.lower, UB = pars.upper, control = list(trace = FALSE))), silent = TRUE)
    if(class(opt) == "try-error")
      return(c(1e99, 1))
    
    c(estimate = opt$par, conv = opt$convergence)
  }
  
  mc_results = replicate(10, mc(), simplify = FALSE)
  
  conv    <- sapply(mc_results, function(x) x["conv"])
  converged_indexes = which(conv == 0)
  est_mat <- t(sapply(mc_results, function(x)
    x[grep("^estimate\\.", names(x))]
  ))
  est_mat_mean = as.vector(apply(est_mat, mean,MARGIN = 2))
  est_mat_sd = as.vector(apply(est_mat, sd,MARGIN = 2))
  
  # 1) at least 80% converged
  conv_rate <- length(converged_indexes)/length(conv)
  expect_gte(conv_rate, 0.8)
  
  # 2) mean close to true vector (1,-2,0.5)
  expect_equal(est_mat_mean[1], 1, tolerance = 0.2)
  expect_equal(est_mat_mean[2], -2, tolerance = 0.2)
  expect_equal(est_mat_mean[3], 0.5, tolerance = 0.1)
  
  # 3) std small
  expect_lt(est_mat_sd[1], 0.1)
  expect_lt(est_mat_sd[2], 0.1)
  expect_lt(est_mat_sd[3], 0.1)
  
})




test_that("mbev_estimation runs when specifying only dataset", {
  expect_no_error(mbev_estimation(Y = rmbev(n = 500)))
  
})




  