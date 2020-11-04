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

})
