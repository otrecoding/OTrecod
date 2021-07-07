test_that("verif_OT works", {


  data(simu_data)

  # with two imputed databases from OT_outcome
  test1 = OT_outcome(simu_data[c(1:150,301:450),], quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
                    dist.choice = "G",percent.knn = 0.90, maxrelax = 0,
                    convert.num = 8, convert.clss = 3,
                    indiv.method = "sequential", which.DB = "BOTH",prox.dist = 0.30)

  ver1 = verif_OT(test1)

  # Structure of the expected output object

  # expect_that(ver1,is_a("list"))
  expect_equal(length(ver1),7)

  expect_equal(is.numeric(ver1[[1]]),TRUE)

  expect_equal(nrow(ver1[[2]]),length(levels(simu_data[,2]))+1)
  expect_equal(ncol(ver1[[2]]),length(levels(as.factor(simu_data[,3])))+1)

  expect_equal(is.matrix(ver1[[3]]),TRUE)
  expect_equal(dim(ver1[[3]]),c(3,3))

  expect_equal(is.null(ver1[[4]]),TRUE)
  expect_equal(dim(ver1[[5]]),c(1,2))
  expect_equal(is.null(ver1[[6]]),TRUE)
  expect_equal(is.null(ver1[[7]]),TRUE)


  # two overlayed databases, all options
  ver2 = verif_OT(test1, group.clss = TRUE, stab.prob = TRUE)

  expect_equal(ncol(ver2[[4]]),5)
  expect_equal(!is.null(ver2[[6]]),TRUE)
  expect_equal(dim(ver2[[7]]),c(3,4))

  # One database (A), all options
  test2 = OT_outcome(simu_data[c(1:150,301:450),], quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
                     dist.choice = "G",percent.knn = 0.90, maxrelax = 0,
                     convert.num = 8, convert.clss = 3,
                     indiv.method = "sequential", which.DB = "A",prox.dist = 0.30)

  ver3 = verif_OT(test2, group.clss = TRUE, stab.prob = TRUE)

  # expect_equal(length(ver3[[3]]),3)
  expect_equal(dim(ver3[[5]]),c(1,2))
  expect_equal(is.na(ver3[[5]][,1]),TRUE)
  expect_equal(dim(ver3[[7]]),c(1,4))


  # One database (B), all options

  # test3 = OT_outcome(simu_data, quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
  #                   dist.choice = "G",percent.knn = 0.90, maxrelax = 0,
  #                   convert.num = 8, convert.clss = 3,
  #                   indiv.method = "sequential",which.DB = "B",prox.dist = 0.30)

  # ver4 = verif_OT(test3, group.clss = TRUE, stab.prob = TRUE)

  # expect_equal(length(ver4[[4]]),3)
  # expect_equal(dim(ver4[[6]]),c(1,2))
  # expect_equal(is.na(ver4[[6]][,2]),TRUE)
  # expect_equal(dim(ver4[[8]]),c(1,5))


})
