test_that("OT_joint works", {

  data(tab_test)
  tb_test     = tab_test[c(1:80,5001:5090),]
  tb_test[,2] = factor(tb_test[,2])
  tb_test[,3] = ordered(tb_test[,3])

  test1 = OT_joint(tb_test, dist.choice = "M", percent.knn = 0.30, which.DB = "A", ordinal = 1:6, prox.X = 0.10, maxrelax = 0.4, lambda.reg = 0.1)

  expect_that(test1,is_a("otres"))
  expect_equal(length(test1),9)

  expect_equal(sum(test1[[2]]),1)
  # expect_equal(sum(test1[[3]]),1)

  # expect_equal(dim(test1[[2]]),dim(test1[[3]]))
  expect_equal(nrow(test1[[2]]),length(levels(tb_test[,2])))
  # expect_equal(ncol(test1[[3]]),length(levels(tb_test[,3])))

  expect_equal(ncol(test1[[4]]),ncol(tab_test)-2)

  expect_that(test1[[5]],is_a("list"))

  expect_identical(levels(test1[[8]][,ncol(test1[[8]])]),levels(tb_test[,3]))
  # expect_identical(levels(test1[[9]][,ncol(test1[[9]])]),levels(tb_test[,2]))
  # expect_equal(ncol(test1[[8]]),ncol(test1[[9]]))


  # test2 = OT_joint(tb_test, dist.choice = "G", percent.knn = 0.40, which.DB = "B", ordinal = 1:6, prox.X = 0.20)

  # expect_null(test2[[2]])
  # expect_equal(sum(test1[[3]]),1)
  # expect_null(test2[[6]])
  # expect_equal(ncol(test2[[8]])+1,ncol(test2[[9]]))
  # expect_that(test2[[7]],is_a("array"))
  # expect_equal(dim(test2[[7]][,,1]),c(nrow(test2[[4]]),length(levels(tb_test[,3]))))
  # expect_identical(levels(test2[[9]][,ncol(test2[[9]])]),levels(tb_test[,2]))


  # Using the maxrelax and lambda.reg arguments

  # test3 = OT_joint(tb_test, dist.choice = "E", percent.knn = 0.40, which.DB = "BOTH", ordinal = 1:6, prox.X = 0.10, maxrelax = 0.4, lambda.reg = 0.1)

  # expect_equal(sum(test3[[2]]),1)
  # expect_equal(sum(test3[[3]]),1)
  # expect_equal(ncol(test3[[8]]),ncol(test3[[9]]))



})
