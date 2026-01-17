test_that("Passing invalid to my function", {
  # passing invalid sample size
  expect_error(rmbev(n = 1))
  
  # invalid delta
  expect_error(rmbev(delta1 = -1))
  expect_error(rmbev(delta2 = -1))
  
  # invalid sigma
  expect_error(rmbev(sigma1 = -1))
  expect_error(rmbev(sigma2 = -1))
  expect_error(rmbev(sigma1 = 0))
  expect_error(rmbev(sigma2 = 0))
  
  # invalid model (currently only log allowed)
  expect_error(rmbev(modelo = "alog"))
  
  # invalid dep
  expect_error(rmbev(dep = 0))
  expect_error(rmbev(dep = 0.9999))
  expect_error(rmbev(dep = -100))
  
})


test_that("Output format", {
  
  
  # (nx2) matrix of output 
  expect_equal( dim(rmbev(n = 2)), c(2,2))
  
  
  # (nx2) matrix of output 
  expect_equal( dim(rmbev(n = 43)), c(43,2))
})


test_that("Checking correct behaviour for knonw values", {
  
  # set.seed fixes generated value
  set.seed(1)
  xy = rmbev(n = 2)
  expect_equal(as.vector(xy[2,2]), 2.788448, tolerance = 0.0001)
})


test_that("Edge cases behaviour", {
  
  # xi close to zero gives valid values
  expect_silent(rmbev(xi1 = 1e-30))
  
  # Inf values appear when xi is big
  x1 = rmbev(n = 10, xi1 = 1e30)[,1]
  expect_gt(sum(x1 == Inf), 1)
  x2 = rmbev(n = 10, xi2 = 1e30)[,2]
  expect_gt(sum(x2 == Inf), 1)
  x = rmbev(n = 10, xi1 = 1e30, xi2 = 1e30)
  expect_gt(sum(x[,1] == Inf), 1)
  expect_gt(sum(x[,2] == Inf), 1)
  
  # abs(Observations) get close to 1 if delta is big 
  expect_equal(abs(as.vector(rmbev(n = 10, delta1 = 100000000)[10,1])), 1, tolerance = 0.0001)
  
})






