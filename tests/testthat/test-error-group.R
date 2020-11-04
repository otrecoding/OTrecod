test_that("error_group works", {

  Z1 = as.factor(sample(1:3,50,replace = TRUE))
  Z2 = sample(c("A","B","C","D"),50, replace = TRUE)
  Z3 = as.factor(sample(1:2,50,replace = TRUE))
  Z4 = as.factor(sample(c("A","B","C","D","E"),50, replace = TRUE))

  expect_equal(is.data.frame(error_group(Z1,Z4)),TRUE)
  expect_error(error_group(Z4,Z1))
  expect_error(error_group(Z4,Z2))
  expect_error(error_group(Z3,Z3))
  expect_equal(ncol(error_group(Z1,Z4)),5)
  expect_lt(nrow(error_group(Z1,Z4)),nrow(error_group(Z1,Z4,FALSE)))


})
