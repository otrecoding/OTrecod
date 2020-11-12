test_that("transfo_quali works", {

  treatment  = as.factor(c(rep("trt1",12),rep("trt2",15),rep("trt3",13)))
  treatment2 = c(rep("trt1",12),rep("trt2",15),rep("trt3",13))
  treatment3 = treatment
  treatment3[c(5,12,27)] = NA
  treat_bin  = transfo_quali(treatment,"trt")
  treat_bin3 = transfo_quali(treatment3,"trt")

  expect_equal(is.matrix(treat_bin), TRUE)
  expect_equal(dim(treat_bin),c(length(treatment),length(levels(treatment))-1))
  expect_error(transfo_quali(treatment2,"trt"))
  expect_error(transfo_quali(treatment,labx = trt))
  expect_gt(sum(is.na(treat_bin3)), 0)

  # Factor with only one level
  treatment4 = as.factor(rep("A",25))
  treat_bin2 = transfo_quali(treatment4,"trt")
  treat_bin4 = transfo_quali(treatment4)

  expect_equal(dim(treat_bin2),c(25,1))
  expect_identical(colnames(treat_bin2),"trt")

  expect_equal(dim(treat_bin4),c(25,1))
  expect_null(colnames(treat_bin4))

  treatment5 = as.factor(rep(0,25))
  treat_bin5 = transfo_quali(treatment5,"trt")
  expect_equal(sum(treat_bin5[,1]),0)

})
