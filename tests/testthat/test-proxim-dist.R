test_that("proxim_dist works", {

  data(simu_data)
  tab1  = transfo_dist(simu_data[c(1:150,301:500),],quanti = c(3,8), nominal = c(1,4:5,7),
                       ordinal = c(2,6), logic = NULL, prep_choice = "M")


  ## test1: Manhattan dist. with NAs, standard order for DB, Y and Z

  test1 = proxim_dist(tab1, norm = "M")

  # global

  expect_that(test1,is_a("list"))
  expect_equal(length(test1),16)


  # by output object

  expect_true(is.character(test1[[1]]))
  expect_equal(test1[[2]],nrow(tab1[tab1[,1]=="A",]))
  expect_equal(test1[[3]],nrow(tab1[tab1[,1]=="B",]))
  expect_equal(dim(test1[[4]]),c(nrow(tab1),ncol(tab1)-3))
  expect_lte(nrow(test1[[5]]),nrow(test1[[4]]))

  expect_equal(as.vector(table(test1[[6]])),as.vector(table(tab1[,2])))
  expect_equal(as.vector(table(test1[[7]])),as.vector(table(tab1[,3])))

  expect_equal(is.matrix(test1[[8]]),TRUE)
  expect_equal(dim(test1[[8]]),c(test1[[2]],test1[[3]]))

  expect_length(test1[[9]],length(levels(tab1[,2])))
  expect_length(test1[[10]],length(levels(as.factor(tab1[,3]))))

  expect_that(test1[[11]],is_a("list"))
  expect_equal(sapply(test1[[11]],length),as.vector(table(tab1[,2])))

  expect_that(test1[[12]],is_a("list"))
  expect_equal(sapply(test1[[12]],length),as.vector(table(tab1[,3])))

  expect_that(test1[[13]],is_a("list"))
  expect_equal(length(test1[[13]]),nrow(test1[[5]]))

  expect_that(test1[[14]],is_a("list"))
  expect_equal(length(test1[[14]]),nrow(test1[[5]]))

  expect_equal(is.matrix(test1[[15]]),TRUE)
  expect_equal(dim(test1[[15]]),c(test1[[2]],test1[[2]]))

  expect_equal(is.matrix(test1[[16]]),TRUE)
  expect_equal(dim(test1[[16]]),c(test1[[3]],test1[[3]]))


  ## varying prox argument

  # test2 = proxim_dist(tab1, norm = "E", prox = 0.20)

  # expect_lt(sum(sapply(test2[[13]],length)),sum(sapply(test1[[13]],length)))
  # expect_lt(sum(sapply(test2[[14]],length)),sum(sapply(test1[[14]],length)))

  # test2b = proxim_dist(tab1, norm = "E", prox = 0)

  # expect_lt(sum(sapply(test2b[[13]],length)),sum(sapply(test1[[13]],length)))
  # expect_lt(sum(sapply(test2b[[14]],length)),sum(sapply(test1[[14]],length)))


  ## varying order for DB, Y and Z

  tab2       = tab1[,c(2,4,1,5:8,3,9)]
  test3      = proxim_dist(tab2,indx_DB_Y_Z = c(3,1,8),norm = "E")
  expect_identical(test3[[4]],test1[[4]])


  ## varying the remaining distances

  # Gower

  # tab3  = transfo_dist(simu_data,quanti = c(3,8), nominal = c(1,4:5,7),
  #                     ordinal = c(2,6), logic = NULL, prep_choice = "G")

  # test4 = proxim_dist(tab3, norm = "G")

  # expect_equal(dim(test4[[8]]),c(test1[[2]],test1[[3]]))
  # expect_equal(dim(test4[[15]]),c(test1[[2]],test1[[2]]))
  # expect_equal(dim(test4[[16]]),c(test1[[3]],test1[[3]]))

  # Hamming

  # sim_data = simu_data[c(1:100,301:450),1:7]

  # tab4  = transfo_dist(sim_data,quanti = 3, nominal = c(1,4:5,7),ordinal = c(2,6),prep_choice = "H")
  # test5 = proxim_dist(tab4,norm = "H")

  # expect_equal(dim(test5[[4]]),c(nrow(tab4),ncol(tab4)-3))
  # expect_lte(nrow(test5[[5]]),nrow(test5[[4]]))
  # expect_equal(dim(test5[[8]]),c(test5[[2]],test5[[3]]))
  # expect_equal(length(test5[[13]]),nrow(test5[[5]]))
  # expect_equal(length(test5[[14]]),nrow(test5[[5]]))
  # expect_equal(dim(test5[[15]]),c(test5[[2]],test5[[2]]))
  # expect_equal(dim(test5[[16]]),c(test5[[3]],test5[[3]]))

  # Only on covariate

  tab5  = tab1[,c(1:3,8)]
  tab6  = tab5[!is.na(tab5[,4]),]
  test6 = proxim_dist(tab6, norm = "G")

  expect_equal(dim(test6[[15]]),c(test6[[2]],test6[[2]]))
  expect_equal(dim(test6[[16]]),c(test6[[3]],test6[[3]]))


})
