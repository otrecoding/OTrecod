test_that("the dimensions of the output objects are validated", {

  # complete case: global and "by output object" tests
  Y      = rnorm(150,50,10)
  test1  = transfo_target(Y)

  expect_that(test1,is_a("list"))
  expect_equal(length(test1),2)
  expect_equal(length(test1[[1]]),length(Y))
  expect_equal(length(test1),2)

  # with NAs: "by output object" tests

  Z      = Y; Z[c(5,32,26,112)] = NA
  test2  = transfo_target(Z)

  expect_equal(length(test2),2)
  expect_equal(length(test2[[1]]),length(Z))

})


test_that("the formats are respected", {

  # continuous variable
  Y      = rnorm(150,50,10)
  Z      = c(rep(1:5,10))
  test1  = transfo_target(Y)
  lev1   = c("5","3","2","4","1")
  test1b = transfo_target(Z, levels_order = lev1)

  expect_equal(is.factor(test1[[1]]),TRUE)
  expect_equal(is.ordered(test1b[[1]]),TRUE)
  expect_identical(test1b[[2]],lev1)

  # character variable
  YY     = c(rep("A",27),rep("B",42),rep("C",36))
  test2  = transfo_target(YY)

  expect_equal(is.factor(test2[[1]]),TRUE)
  expect_equal(levels(test2[[1]]),c("A","B","C"))

  test3  = transfo_target(YY,levels_order = c("C","B","A"))
  expect_equal(is.ordered(test3[[1]]),TRUE)
  expect_equal(levels(test3[[1]]),c("C","B","A"))


  # factor
  ZZ    = as.factor(YY)
  test4 = transfo_target(ZZ)
  expect_equal(is.factor(test4[[1]]),TRUE)

  test5  = transfo_target(ZZ,levels_order = c("C","B","A"))
  expect_equal(is.ordered(test5[[1]]),TRUE)

  test6  = transfo_target(ZZ,levels_order = c("A","B","D","C"))
  expect_identical(test6[[2]],levels(ZZ))

  # ordered factor
  UU    = as.ordered(YY)
  test7 = transfo_target(UU)
  expect_equal(is.ordered(test7[[1]]),TRUE)


})
