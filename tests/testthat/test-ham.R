test_that("ham works", {

  # Complete Case
  aaa = sample(c(0,1),12,replace = TRUE)
  bbb = sample(c(0,1),21,replace = TRUE)
  A = matrix(aaa, ncol = 3)
  B = matrix(bbb, ncol = 3)
  C = c(1.74,6.23,5.80)

  # Incomplete case
  A_NA        = A
  A_NA[3,1]   = NA
  A_NA[2,2:3] = rep(NA,2)

  B_NA = B
  B_NA[3,1:3] = NA

  expect_equal(dim(ham(A,B)),c(nrow(A),nrow(B)))
  expect_equal(dim(ham(A,A)),c(nrow(A),nrow(A)))
  expect_identical(is.vector(ham(C,A)),TRUE)
  expect_that(ham(A,C),is_a("matrix"))
  expect_equal(dim(ham(C,C)),rep(length(C),2))



})
