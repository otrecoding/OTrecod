test_that("mice works",{

  data(simu_data)
  sim_data = simu_data[c(1:100,301:450),]

  ind1 = 4:8
  imp1 = imput_cov(sim_data, indcol = ind1, R_mice = 2, meth = c("logreg","polyreg","polr","logreg","pmm"))

  # global
  # expect_that(imp1,is_a("list"))
  expect_equal(length(imp1),4)

  # output 1
  expect_equal(dim(imp1[[1]]),dim(sim_data))

  # output 2
  expect_identical(imp1[[2]],"MICE")

  # output 3
  expect_equal(ncol(imp1[[3]]),length(ind1))
  expect_equal(nrow(imp1[[3]]),nrow(sim_data))
  expect_identical(colnames(imp1[[3]]),colnames(sim_data)[ind1])
  # expect_identical(lapply(imp1[[3]],typeof),lapply(simu_data[ind1],typeof))

  # output 4
  expect_equal(length(imp1[[4]]),2)
  expect_identical(colnames(imp1[[4]][[1]]),colnames(sim_data)[ind1])


  # ind3 = c(4,6:8)
  # imp3 = imput_cov(sim_data, indcol = ind3, R_mice = 2, meth = c("logreg","polr","logreg","pmm"))

  # expect_identical(colnames(imp3[[3]]),colnames(sim_data)[ind3])

  ind4 = c(6:8,4)
  imp4 = imput_cov(sim_data, indcol = ind4, R_mice = 2, meth = c("polr","logreg","pmm","logreg"))

  expect_identical(colnames(imp4[[3]]),colnames(sim_data)[ind4])


})


test_that("famd works",{

  data(simu_data)

  ind2 = 4:8
  imp2 = imput_cov(simu_data,indcol = ind2, meth = c("logreg","polyreg","polr","logreg","pmm"),
                   missMDA = TRUE, NB_COMP = 4)

  # global
  # expect_that(imp2,is_a("list"))
  expect_equal(length(imp2),3)

  # output 2
  expect_identical(imp2[[2]],"FAMD")

  # output 3
  expect_equal(ncol(imp2[[3]]),length(ind2))
  expect_equal(nrow(imp2[[3]]),nrow(simu_data))
  expect_identical(colnames(imp2[[3]]),colnames(simu_data)[ind2])

})

# test_that("there are no discrepancies between arguments",{
#
#  data(simu_data)
#
#  typ_mic = c("logreg","polyreg","polr")
#  expect_error(imput_cov(simu_data, indcol = 4:8, R_mice = 3, meth = typ_mic))
#
#})




