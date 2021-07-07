test_that("OT_outcome works", {

  data(simu_data)
  simu_dat  = simu_data[c(1:200,301:520),]

  ## individual method: optimal with 2 databases

  # Global
  test1 = OT_outcome(simu_dat, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
                     dist.choice = "M", percent.knn = 0.90, maxrelax = 0,
                     indiv.method = "sequential")  # Ex - optimal

  expect_that(test1,is_a("otres"))
  expect_equal(length(test1),9)

  expect_equal(dim(test1[[2]]),dim(test1[[3]]))
  expect_equal(nrow(test1[[2]]),length(levels(simu_data[,2])))
  expect_equal(ncol(test1[[2]]),length(levels(as.factor(simu_data[,3]))))

  # expect_that(test1[[5]],is_a("list"))
  expect_equal(nrow(test1[[5]][[5]]),nrow(test1[[4]]))

  expect_that(test1[[6]],is_a("array"))
  expect_that(test1[[7]],is_a("array"))

  expect_identical(levels(test1[[8]][,10]),levels(as.factor(simu_data[,3])))
  expect_identical(levels(test1[[9]][,10]),levels(simu_data[,2]))

  expect_equal(nrow(test1[[8]]),nrow(simu_dat[simu_dat[,1]=="A",]))
  expect_equal(nrow(test1[[9]]),nrow(simu_dat[simu_dat[,1]=="B",]))

  expect_equal(sum(is.na(test1[[8]][,10])),0)
  expect_equal(sum(is.na(test1[[9]][,10])),0)


  ## individual method: optimal with 1 database only A

  test2 = OT_outcome(simu_dat, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
                     dist.choice = "M", percent.knn = 0.60, maxrelax = 0.40,
                     indiv.method = "optimal", which.DB = "A")

  expect_equal(length(test2),9)
  expect_equal(dim(test2[[2]]),dim(test2[[3]]))
  expect_null(test2[[7]])

  expect_equal(sum(is.na(test2[[8]][,10])),0)
  expect_lt(ncol(test2[[9]]),ncol(test2[[8]]))


  ## individual method: optimal with 1 database only B

  # test3 = OT_outcome(simu_dat, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
  #                    dist.choice = "G", percent.knn = 0.60, maxrelax = 0.40,
  #                   indiv.method = "optimal", which.DB = "B")

  # expect_equal(length(test3),9)
  # expect_equal(dim(test3[[2]]),dim(test3[[3]]))
  # expect_null(test3[[6]])

  # expect_equal(sum(is.na(test3[[9]][,9])),0)
  # expect_lt(ncol(test3[[8]]),ncol(test3[[9]]))

  # expect_identical(levels(test3[[9]][,9]),levels(simu_data[,2]))


})
