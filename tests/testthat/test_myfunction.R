test_that("Passing invalid to my function", {
  # passing two input values instead of one
  expect_error(my_function(3,3))
  
})


test_that("Checking correct behaviour for knonw values", {
  # sd cannot be <= 0
  expect_equal(my_function(3), 3.001, tolerance = 0.01)
})



