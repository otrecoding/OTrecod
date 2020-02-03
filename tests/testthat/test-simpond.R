context("test-simpond")

test_that("simpond function works", {
  x = c(0.5,0,0.5,0.5,3.5)
  simplex_new = x
  expect_equal(simplex_new , c(0.5,0,0.5,0.5,3.5))
})
