test_that("select_pred works", {

  data(simu_data)
  sim_data = simu_data[c(1:150,301:520),]

  ## overlayed without RF

  test1 = select_pred(sim_data,Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Y",
                      quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
                      thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
                      RF = FALSE)


  # expect_that(test1,is_a("list"))
  expect_equal(length(test1),11)

  expect_true(is.numeric(test1[[1]]))
  expect_length(test1[[1]],1)
  expect_identical(test1[[2]],"Yb1")
  expect_length(test1[[4]],0)
  expect_equal(dim(test1[[5]]),dim(sim_data[sim_data[,1]=="A",]))
  expect_that(test1[[6]],is_a("data.frame"))
  expect_equal(ncol(test1[[6]]),5)
  expect_that(test1[[7]],is_a("data.frame"))
  expect_equal(ncol(test1[[7]]),5)
  expect_that(test1[[8]],is_a("data.frame"))
  expect_equal(ncol(test1[[8]]),5)
  expect_that(test1[[9]],is_a("data.frame"))
  expect_equal(ncol(test1[[9]]),5)
  expect_length(test1[[10]],3)
  # expect_that(test1[[11]],is_a("list"))
  expect_true(all(test1[[11]][[1]][,3]>0.30))


  ## overlayed with RF

  test2 = select_pred(sim_data,Y = "Yb1", Z = "Yb2", ID = 1, OUT = "Z",
                      quanti = c(3,8), nominal = c(1,4:5,7), ordinal = c(2,6),
                      thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
                      RF = TRUE, RF_SEED = 448571)

  expect_equal(test2[[1]],448571)
  # expect_that(test2,is_a("list"))
  expect_equal(length(test2),14)


  expect_identical(test2[[2]],"Yb2")
  expect_length(test2[[4]],0)
  expect_equal(dim(test2[[5]]),dim(sim_data[sim_data[,1]=="B",]))
  # expect_true(is.vector(test2[[13]]))
  expect_identical(names(test2[[13]])[test2[[13]]> 20],test2[[14]])


  ## convert + ident + logical

  sim_data2      = sim_data[,c(2:4,1,5:8)]
  sim_data2$logi = rep(c(TRUE,FALSE),nrow(sim_data2)/2)

  # test3 = select_pred(sim_data2,Y = "Yb1", Z = "Yb2", ID = 4, OUT = "Y",
  #                    quanti = c(2,8), nominal = c(4,3,5,7), ordinal = c(1,6),
  #                    logic = 9, convert_num = 8, convert_clss = 4,
  #                    thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
  #                    RF = TRUE)

  # expect_equal(dim(test3[[5]]),dim(sim_data2[sim_data2[,4]=="A",]))


  ## single database

  sim_data3 = sim_data2[sim_data2[,4]=="B",]
  sim_data3 = sim_data3[,-1]

  test4 = select_pred(sim_data3,Y = "Yb2", Z = NULL, ID = 3, OUT = "Y",
                      quanti = c(1,7), nominal = c(3,2,4,6), ordinal = 5,
                      logic = 8, convert_num = 7, convert_clss = 4,
                      thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
                      RF = TRUE)

  expect_equal(dim(test4[[5]]),dim(sim_data3))

  # RF condi

  # ident_NA  = unique(unlist(lapply(sim_data2[,3:9],function(x)which(is.na(x)))))
  # sim_data4 = sim_data2[-ident_NA,]

  # test5 = select_pred(sim_data4,Y = "Yb1", Z = "Yb2", ID = 4, OUT = "Y",
    #                  quanti = c(2,8), nominal = c(4,3,5,7), ordinal = c(1,6),
    #                  logic = 9, convert_num = 8, convert_clss = 4,
    #                  thresh_cat = 0.30, thresh_num = 0.70, thresh_Y = 0.20,
    #                  RF = TRUE, RF_condi = TRUE, RF_SEED = 2036)

  # expect_equal(length(test5),14)
  # expect_equal(test5[[1]],2036)
  # expect_equal(dim(test5[[5]]),dim(sim_data4[sim_data4[,4]=="A",]))


})
