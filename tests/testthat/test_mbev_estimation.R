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


test_that("mbev_estimation runs when specifying only dataset", {
  set.seed(1)
  expect_no_error(mbev_estimation2(Y = rmbev(n = 500)))
})




  