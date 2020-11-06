test_that("the dimensions of the output objects are validated", {

  # Build situation

  data(simu_data)
  data_A     = simu_data[1:200  ,c(2,4:8)]
  data_B     = simu_data[301:500,c(7,4,6,5,8,3)]
  ident      = paste("L",1:nrow(data_B),sep="_")
  data_B$ident = ident
  data_B     = data_B[,c(1:4,7,5:6)]
  data_B$Yb2 = as.factor(data_B$Yb2)

  # Complete outcomes and MICE

  test1   = merge_dbs(data_A,data_B, row_ID2 = 5,
                      NAME_Y = "Yb1",NAME_Z = "Yb2",
                      ordinal_DB1 = c(1,4), ordinal_DB2 = c(3,7),
                      impute = "MICE",R_MICE = 2, seed_func = 2067)

  expect_equal(length(test1),13)
  expect_equal(nrow(test1[[1]]),400)
  expect_equal(ncol(test1[[1]]),ncol(simu_data))
  expect_equal(sum(is.na(test1[[1]][,2])),nrow(data_B))
  expect_equal(sum(is.na(test1[[1]][,3])),nrow(data_A))


  expect_equal(length(test1[[2]]),0)
  expect_equal(length(test1[[3]]),0)

  expect_identical(test1[[4]],levels(data_A$Yb1))
  expect_equal(test1[[13]],2067)

  expect_equal(is.null(test1[[6]]),TRUE)
  expect_equal(is.null(test1[[7]]),TRUE)

  expect_that(test1[[10]],is_a("list"))
  expect_equal(length(test1[[10]]),2)

  expect_equal(length(test1[[10]][[1]]),2)
  expect_equal(nrow(test1[[10]][[1]][[2]][[1]]),200)

  # Incomplete outcomes and MICE

  data_A$Yb1[c(5,10,125,42)] = NA
  data_B$Yb2[c(21,151,182)]  = NA

  test2   = merge_dbs(data_A,data_B, row_ID2 = 5,
                      NAME_Y = "Yb1",NAME_Z = "Yb2",
                      ordinal_DB1 = c(1,4), ordinal_DB2 = c(3,7),
                      impute = "MICE",R_MICE = 2, seed_func = 2067)


  expect_equal(dim(test2[[1]]),c(400-7,8))
  expect_equal(length(test2[[2]]),4)
  expect_equal(length(test2[[3]]),3)
  expect_equal(dim(test2[[11]]),c(200-4,6))
  expect_equal(dim(test2[[12]]),c(200-3,7))

  expect_equal(nrow(test2[[10]][[1]][[2]][[1]]),200-4)


  # Incomplete outcomes and FAMD

  test3   = merge_dbs(data_A,data_B, row_ID2 = 5,
                      NAME_Y = "Yb1",NAME_Z = "Yb2",
                      ordinal_DB1 = c(1,4), ordinal_DB2 = c(3,7),
                      impute = "FAMD", NCP_FAMD = 4, seed_func = 3078)

  expect_equal(length(test3),12)
  expect_equal(dim(test3[[1]]),c(400-7,8))


  # Incomplete outcomes and NO

  test4   = merge_dbs(data_A, data_B, row_ID2 = 5,
                      NAME_Y = "Yb1",NAME_Z = "Yb2",
                      ordinal_DB1 = c(1,4), ordinal_DB2 = c(3,7),
                      impute = "CC", seed_func = 3078)

  expect_equal(length(test4),12)
  expect_equal(dim(test4[[1]]),c(304,8))

  # MICE and REMOVE_VAR

  data_A$Weight1   = stats::rnorm(200,70,5)
  data_B$Weight2   = stats::rnorm(200,65,5)
  data_B$Treatment = as.factor(ifelse(data_B$Treatment == "Trt A","trt A",data_B$Treatment))
  data_A$Treatment = as.factor(data_A$Treatment)
  data_A$Smoking = as.character(data_A$Smoking)

  test5   = merge_dbs(data_A, data_B, row_ID2 = 5,
                      NAME_Y = "Yb1",NAME_Z = "Yb2",
                      ordinal_DB1 = c(1,4), ordinal_DB2 = c(3,7),
                      impute = "CC", seed_func = 3078)

  expect_equal(length(test5$REMOVE1),1)
  expect_equal(length(test5$REMOVE2),1)
  expect_equal(nrow(test5[[1]]),304)


})
