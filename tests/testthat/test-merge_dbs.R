test_that("the dimensions of the output objects are validated", {

  # From two distinct databases

  data(simu_data)
  data_A     = simu_data[1:200  ,c(2,4:8)]
  data_B     = simu_data[301:500,c(7,4,6,5,8,3)]
  data_B$Yb2 = as.factor(data_B$Yb2)

  test1   = merge_dbs(data_A,data_B,
                      NAME_Y = "Yb1",NAME_Z = "Yb2",
                      ordinal_DB1 = c(1,4), ordinal_DB2 = c(3,6),
                      impute = "MICE",R_MICE = 2, seed_func = 2067)

  expect_equal(length(test1),13)
  expect_equal(is.null(test1[[2]]),TRUE)
  expect_identical(test1[[4]],levels(data_A$Yb1))
  expect_equal(test1[[13]],2067)

})
