test_that("transfo_dist works", {

  data(simu_data)

  sim_data     = simu_data
  sim_data$Yb2 = as.ordered(sim_data$Yb2)

  test1 = transfo_dist(sim_data, quanti = c(3,8), nominal = c(1,4:5,7),
                       ordinal = c(2,6), logic = NULL, prep_choice = "E")

  # Euclidean distance ok
  expect_that(test1,is_a("data.frame"))
  expect_equal(nrow(test1),nrow(sim_data))
  expect_equal(is.numeric(test1[,8]),TRUE)
  expect_equal(sum(is.na(test1[,ncol(test1)])),sum(is.na(sim_data[,ncol(sim_data)])))


  # Manhattan distance ok
  test2 = transfo_dist(sim_data, quanti = c(3,8), nominal = c(1,4:5,7),
                       ordinal = c(2,6), logic = NULL, prep_choice = "M")

  expect_identical(test1,test2)


  # Gower distance ok
  test3 = transfo_dist(sim_data, quanti = c(3,8), nominal = c(1,4:5,7),
                       ordinal = c(2,6), logic = NULL, prep_choice = "G")
  #
  # expect_equal(dim(test3),dim(sim_data))

  # FAMD ok
  # test4 = transfo_dist(sim_data, quanti = c(3,8), nominal = c(1,4:5,7),
  #                     ordinal = c(2,6), logic = NULL, prep_choice = "FAMD")

  test5 = transfo_dist(sim_data, quanti = c(3,8), nominal = c(1,4:5,7),
                       ordinal = c(2,6), logic = NULL, prep_choice = "FAMD",info = 0.5)

  expect_equal(dim(test3),dim(sim_data))
  expect_lt(ncol(test5),ncol(test2))

  # Boolean ok

  sim_data2      = sim_data
  sim_data2$logi = sample(c(T,F),nrow(sim_data),replace = TRUE)

  test6 = transfo_dist(sim_data2, quanti = c(3,8), nominal = c(1,4:5,7),
                       ordinal = c(2,6), logic = 9, prep_choice = "E")


  # Hamming distance and conversion options ok

  test7 = transfo_dist(sim_data,quanti = c(3,8), nominal = c(1,4:5,7),ordinal = c(2,6),
                       convert_num = 8, convert_clss = 3, prep_choice = "H")

  expect_equal(nrow(test7),nrow(sim_data))

  sim_data3 = sim_data
  sim_data3$nw_var = stats::rnorm(nrow(sim_data),10,2)

  # test8 = transfo_dist(sim_data3,quanti = c(3,8,9), nominal = c(1,4:5,7),ordinal = c(2,6),
  #                      convert_num = c(8,9), convert_clss = c(3,4), prep_choice = "H")
  #
  # expect_equal(nrow(test8),nrow(sim_data))



})
