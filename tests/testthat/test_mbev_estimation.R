test_that("Force likelihood to return BIG VALUE", {
  
  # invalid delta1
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,-2,0,1,1,0.5,0.5,1)), 1e99)
  
  # invalid delta2
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,2,-1,1,1,0.5,0.5,1)), 1e99)
  
  # invalid dep
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,2,-1,1,1,0.5,0.5,-0.1)), 1e99)
  
  # invalid dep (depends on my parametrization)
  expect_equal(mbev_log_likelihood(rmbev(), pars = c(0,0,2,-1,1,1,0.5,0.5,1)), 1e99) 
  
})


