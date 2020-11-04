test_that("power_set works", {

  expect_equal(is.list(power_set(5)),TRUE)
  expect_equal(length(power_set(5)),2^5-1)

  expect_error(power_set(0))
  expect_lt(length(power_set(4,ordinal = TRUE)),2^4-1)

})
