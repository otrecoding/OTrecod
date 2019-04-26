context("test-simpond")

test_that("simpond function works", {
  x = c(0.5,0,0.5,0.5,3.5)
  simplex_new = simpond(x)
  expect_equal(simplex_new , c(1.0, 0.0, 1.0, 0.0, 3.0))
})
